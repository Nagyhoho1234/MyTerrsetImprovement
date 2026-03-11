#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class MultiRegModule : public Module {
public:
    QString name() const override { return "MULTI_REG"; }
    QString description() const override {
        return "Multiple linear regression per pixel across time. Regresses a dependent variable "
               "series on multiple independent variable series. Outputs R-squared, coefficients, "
               "and residuals.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dependent_bands", "Dependent variable bands (comma-separated)",
                "Raster bands for the dependent variable time series"),
            ParameterDef::file("independent_bands", "Independent variable bands (comma-separated)",
                "Raster bands for independent variable series. If multiple predictors, "
                "list all bands; the number of time steps is inferred from dependent_bands."),
            ParameterDef::output("output_prefix", "Output prefix"),
        };
    }

    bool execute() override {
        QStringList depFiles = parameter("dependent_bands").toString().split(",", Qt::SkipEmptyParts);
        QStringList indFiles = parameter("independent_bands").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : depFiles) f = f.trimmed();
        for (auto& f : indFiles) f = f.trimmed();

        int nTime = depFiles.size();
        QString outPrefix = parameter("output_prefix").toString();

        if (nTime < 2) {
            setError("At least 2 time steps required for regression");
            return false;
        }

        // Determine number of predictors: indFiles should be nPred * nTime
        if (indFiles.size() % nTime != 0) {
            setError(QString("Number of independent bands (%1) must be a multiple of "
                             "dependent bands count (%2)")
                     .arg(indFiles.size()).arg(nTime));
            return false;
        }
        int nPred = indFiles.size() / nTime;
        int nCoeffs = nPred + 1; // intercept + predictors

        if (nTime < nCoeffs) {
            setError("More time steps required than coefficients (predictors + 1)");
            return false;
        }

        reportProgress(0.0, "Reading dependent variable bands...");

        std::vector<std::unique_ptr<Raster>> depRasters;
        depRasters.reserve(nTime);
        for (int t = 0; t < nTime; ++t) {
            auto r = GdalIO::read(depFiles[t]);
            if (!r) { setError(QString("Failed to read: %1").arg(depFiles[t])); return false; }
            if (t > 0 && (r->cols() != depRasters[0]->cols() || r->rows() != depRasters[0]->rows())) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
            depRasters.push_back(std::move(r));
        }

        reportProgress(0.1, "Reading independent variable bands...");

        // Independent rasters: organized as [pred][time]
        std::vector<std::vector<std::unique_ptr<Raster>>> indRasters(nPred);
        for (int p = 0; p < nPred; ++p) {
            indRasters[p].reserve(nTime);
            for (int t = 0; t < nTime; ++t) {
                int fileIdx = p * nTime + t;
                auto r = GdalIO::read(indFiles[fileIdx]);
                if (!r) { setError(QString("Failed to read: %1").arg(indFiles[fileIdx])); return false; }
                if (r->cols() != depRasters[0]->cols() || r->rows() != depRasters[0]->rows()) {
                    setError("All input rasters must have the same dimensions");
                    return false;
                }
                indRasters[p].push_back(std::move(r));
            }
        }

        int cols = depRasters[0]->cols(), rows = depRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = depRasters[0]->noDataValue();
        bool hasND = depRasters[0]->hasNoData();

        // Outputs: R-squared, coefficients (nCoeffs bands), residuals (nTime bands)
        Raster r2Out(cols, rows, 1, DataType::Float64);
        r2Out.setGeoTransform(depRasters[0]->geoTransform());
        r2Out.setProjection(depRasters[0]->projection());
        r2Out.setNoDataValue(noData);

        Raster coeffOut(cols, rows, nCoeffs, DataType::Float64);
        coeffOut.setGeoTransform(depRasters[0]->geoTransform());
        coeffOut.setProjection(depRasters[0]->projection());
        coeffOut.setNoDataValue(noData);

        Raster residOut(cols, rows, nTime, DataType::Float64);
        residOut.setGeoTransform(depRasters[0]->geoTransform());
        residOut.setProjection(depRasters[0]->projection());
        residOut.setNoDataValue(noData);

        reportProgress(0.2, "Computing per-pixel multiple regression...");

        for (int64_t i = 0; i < total; ++i) {
            bool valid = true;

            // Gather dependent series
            std::vector<double> y(nTime);
            for (int t = 0; t < nTime && valid; ++t) {
                double val = depRasters[t]->data(0)[i];
                if (hasND && val == noData) valid = false;
                y[t] = val;
            }

            // Gather independent series
            // Design matrix X: nTime x nCoeffs (first col = 1 for intercept)
            std::vector<std::vector<double>> X(nTime, std::vector<double>(nCoeffs));
            for (int t = 0; t < nTime && valid; ++t) {
                X[t][0] = 1.0;
                for (int p = 0; p < nPred; ++p) {
                    double val = indRasters[p][t]->data(0)[i];
                    if (hasND && val == noData) valid = false;
                    X[t][p + 1] = val;
                }
            }

            if (!valid) {
                r2Out.data(0)[i] = noData;
                for (int c = 0; c < nCoeffs; ++c) coeffOut.data(c)[i] = noData;
                for (int t = 0; t < nTime; ++t) residOut.data(t)[i] = noData;
                continue;
            }

            // Compute X^T X and X^T y
            std::vector<std::vector<double>> XtX(nCoeffs, std::vector<double>(nCoeffs, 0.0));
            std::vector<double> XtY(nCoeffs, 0.0);
            for (int a = 0; a < nCoeffs; ++a) {
                for (int t = 0; t < nTime; ++t) XtY[a] += X[t][a] * y[t];
                for (int b = a; b < nCoeffs; ++b) {
                    double sum = 0.0;
                    for (int t = 0; t < nTime; ++t) sum += X[t][a] * X[t][b];
                    XtX[a][b] = sum;
                    XtX[b][a] = sum;
                }
            }

            // Solve normal equations via Gauss-Jordan
            auto beta = solveLinearSystem(XtX, XtY);

            // Compute R-squared
            double yMean = 0.0;
            for (int t = 0; t < nTime; ++t) yMean += y[t];
            yMean /= nTime;

            double ssTot = 0.0, ssRes = 0.0;
            for (int t = 0; t < nTime; ++t) {
                double predicted = 0.0;
                for (int c = 0; c < nCoeffs; ++c) predicted += X[t][c] * beta[c];
                double resid = y[t] - predicted;
                ssRes += resid * resid;
                ssTot += (y[t] - yMean) * (y[t] - yMean);
                residOut.data(t)[i] = resid;
            }

            r2Out.data(0)[i] = (ssTot > 1e-15) ? 1.0 - ssRes / ssTot : 0.0;
            for (int c = 0; c < nCoeffs; ++c) coeffOut.data(c)[i] = beta[c];

            if (i % 500000 == 0)
                reportProgress(0.2 + 0.7 * static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");
        if (!GdalIO::write(r2Out, QString("%1_r2.tif").arg(outPrefix))) {
            setError("Failed to write R-squared output");
            return false;
        }
        if (!GdalIO::write(coeffOut, QString("%1_coefficients.tif").arg(outPrefix))) {
            setError("Failed to write coefficients output");
            return false;
        }
        if (!GdalIO::write(residOut, QString("%1_residuals.tif").arg(outPrefix))) {
            setError("Failed to write residuals output");
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

REGISTER_MODULE(MultiRegModule)

} // namespace aplaceholder
