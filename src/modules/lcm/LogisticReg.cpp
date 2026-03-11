#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class LogisticRegModule : public Module {
public:
    QString name() const override { return "LOGISTIC_REG"; }
    QString description() const override {
        return "Logistic regression for land change transition potential modeling. "
               "Binary logistic regression: P(change) = 1/(1+exp(-XB)). Trains from "
               "change/no-change pixels and driver variables.";
    }
    QString category() const override { return "Land Change Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::raster("change_raster", "Change raster (binary 0/1)",
                "Binary raster where 1 = change, 0 = no change"),
            ParameterDef::file("drivers", "Driver variables (comma-separated)",
                "Raster layers of explanatory/driver variables"),
            ParameterDef::output("output_potential", "Output transition potential map"),
            ParameterDef::integer("max_iterations", "Maximum iterations", 20, 1, 500,
                "Maximum number of IRLS iterations for logistic regression"),
        };
    }

    bool execute() override {
        QString changePath = parameter("change_raster").toString();
        QStringList driverFiles = parameter("drivers").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : driverFiles) f = f.trimmed();
        int nDrivers = driverFiles.size();
        int maxIter = parameter("max_iterations").toInt();
        QString outPath = parameter("output_potential").toString();

        if (nDrivers < 1) {
            setError("At least one driver variable is required");
            return false;
        }

        reportProgress(0.0, "Reading change raster...");

        auto changeRaster = GdalIO::read(changePath);
        if (!changeRaster) { setError("Failed to read change raster"); return false; }

        int cols = changeRaster->cols(), rows = changeRaster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = changeRaster->noDataValue();
        bool hasND = changeRaster->hasNoData();

        reportProgress(0.05, "Reading driver variables...");

        std::vector<std::unique_ptr<Raster>> drivers;
        drivers.reserve(nDrivers);
        for (int d = 0; d < nDrivers; ++d) {
            auto r = GdalIO::read(driverFiles[d]);
            if (!r) { setError(QString("Failed to read: %1").arg(driverFiles[d])); return false; }
            if (r->cols() != cols || r->rows() != rows) {
                setError("All rasters must have the same dimensions");
                return false;
            }
            drivers.push_back(std::move(r));
        }

        reportProgress(0.1, "Extracting training samples...");

        // Collect valid pixels for training
        int nCoeffs = nDrivers + 1; // intercept + drivers
        std::vector<std::vector<double>> X; // design matrix rows
        std::vector<double> y; // response

        for (int64_t i = 0; i < total; ++i) {
            double changeVal = changeRaster->data(0)[i];
            if (hasND && changeVal == noData) continue;

            bool valid = true;
            std::vector<double> row(nCoeffs);
            row[0] = 1.0; // intercept
            for (int d = 0; d < nDrivers; ++d) {
                double val = drivers[d]->data(0)[i];
                if (hasND && val == noData) { valid = false; break; }
                row[d + 1] = val;
            }
            if (!valid) continue;

            X.push_back(row);
            y.push_back((changeVal != 0.0) ? 1.0 : 0.0);
        }

        int64_t nSamples = (int64_t)y.size();
        if (nSamples < nCoeffs + 1) {
            setError("Insufficient valid pixels for logistic regression");
            return false;
        }

        reportProgress(0.15, "Standardizing predictors...");

        // Standardize predictors for numerical stability
        std::vector<double> means(nCoeffs, 0.0), stds(nCoeffs, 1.0);
        for (int j = 1; j < nCoeffs; ++j) { // skip intercept
            double sum = 0.0;
            for (int64_t i = 0; i < nSamples; ++i) sum += X[i][j];
            means[j] = sum / nSamples;
            double ss = 0.0;
            for (int64_t i = 0; i < nSamples; ++i) {
                double d = X[i][j] - means[j];
                ss += d * d;
            }
            stds[j] = std::sqrt(ss / nSamples);
            if (stds[j] < 1e-15) stds[j] = 1.0;
            for (int64_t i = 0; i < nSamples; ++i)
                X[i][j] = (X[i][j] - means[j]) / stds[j];
        }

        reportProgress(0.2, "Fitting logistic regression (IRLS)...");

        // Iteratively Reweighted Least Squares (IRLS)
        std::vector<double> beta(nCoeffs, 0.0);

        for (int iter = 0; iter < maxIter; ++iter) {
            // Compute predictions: p = sigmoid(X * beta)
            std::vector<double> p(nSamples);
            for (int64_t i = 0; i < nSamples; ++i) {
                double eta = 0.0;
                for (int j = 0; j < nCoeffs; ++j) eta += X[i][j] * beta[j];
                eta = std::clamp(eta, -20.0, 20.0);
                p[i] = 1.0 / (1.0 + std::exp(-eta));
                p[i] = std::clamp(p[i], 1e-10, 1.0 - 1e-10);
            }

            // Compute gradient: X^T (y - p)
            std::vector<double> gradient(nCoeffs, 0.0);
            for (int j = 0; j < nCoeffs; ++j)
                for (int64_t i = 0; i < nSamples; ++i)
                    gradient[j] += X[i][j] * (y[i] - p[i]);

            // Compute Hessian: -X^T W X where W = diag(p * (1-p))
            std::vector<std::vector<double>> H(nCoeffs, std::vector<double>(nCoeffs, 0.0));
            for (int a = 0; a < nCoeffs; ++a)
                for (int b = a; b < nCoeffs; ++b) {
                    double sum = 0.0;
                    for (int64_t i = 0; i < nSamples; ++i)
                        sum += X[i][a] * p[i] * (1.0 - p[i]) * X[i][b];
                    H[a][b] = sum;
                    H[b][a] = sum;
                }

            // Solve H * delta = gradient
            auto delta = solveLinearSystem(H, gradient);

            // Update beta
            double maxDelta = 0.0;
            for (int j = 0; j < nCoeffs; ++j) {
                beta[j] += delta[j];
                maxDelta = std::max(maxDelta, std::abs(delta[j]));
            }

            reportProgress(0.2 + 0.5 * (iter + 1.0) / maxIter);

            // Check convergence
            if (maxDelta < 1e-6) break;
        }

        // Un-standardize coefficients for prediction on original scale
        // beta_orig[j] = beta[j] / stds[j] for j >= 1
        // beta_orig[0] = beta[0] - sum(beta[j] * means[j] / stds[j])
        std::vector<double> betaOrig(nCoeffs);
        betaOrig[0] = beta[0];
        for (int j = 1; j < nCoeffs; ++j) {
            betaOrig[j] = beta[j] / stds[j];
            betaOrig[0] -= beta[j] * means[j] / stds[j];
        }

        reportProgress(0.7, "Computing transition potential surface...");

        // Apply model to all pixels
        Raster potentialOut(cols, rows, 1, DataType::Float64);
        potentialOut.setGeoTransform(changeRaster->geoTransform());
        potentialOut.setProjection(changeRaster->projection());
        potentialOut.setNoDataValue(noData);
        auto& potData = potentialOut.data(0);

        for (int64_t i = 0; i < total; ++i) {
            bool valid = true;
            double eta = betaOrig[0];
            for (int d = 0; d < nDrivers; ++d) {
                double val = drivers[d]->data(0)[i];
                if (hasND && val == noData) { valid = false; break; }
                eta += betaOrig[d + 1] * val;
            }

            if (!valid) {
                potData[i] = noData;
                continue;
            }

            eta = std::clamp(eta, -20.0, 20.0);
            potData[i] = 1.0 / (1.0 + std::exp(-eta));

            if (i % 500000 == 0)
                reportProgress(0.7 + 0.25 * static_cast<double>(i) / total);
        }

        reportProgress(0.95, "Writing output...");
        if (!GdalIO::write(potentialOut, outPath)) {
            setError("Failed to write transition potential output");
            return false;
        }

        reportProgress(1.0);
        return true;
    }

private:
    static std::vector<double> solveLinearSystem(
        std::vector<std::vector<double>> A, std::vector<double> b) {
        int n = (int)b.size();
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i + 1; j < n; ++j)
                if (std::abs(A[j][i]) > std::abs(A[pivot][i])) pivot = j;
            std::swap(A[i], A[pivot]);
            std::swap(b[i], b[pivot]);
            double diag = A[i][i];
            if (std::abs(diag) < 1e-15) continue;
            for (int j = i + 1; j < n; ++j) {
                double factor = A[j][i] / diag;
                for (int k = i; k < n; ++k) A[j][k] -= factor * A[i][k];
                b[j] -= factor * b[i];
            }
        }
        std::vector<double> x(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            double sum = b[i];
            for (int j = i + 1; j < n; ++j) sum -= A[i][j] * x[j];
            x[i] = (std::abs(A[i][i]) > 1e-15) ? sum / A[i][i] : 0.0;
        }
        return x;
    }
};

REGISTER_MODULE(LogisticRegModule)

} // namespace aplaceholder
