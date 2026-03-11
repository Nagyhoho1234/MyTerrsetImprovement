#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class MlrModule : public Module {
public:
    QString name() const override { return "MLR"; }
    QString description() const override {
        return "Multi-temporal Linear Regression. OLS regression per pixel over time with "
               "multiple predictors. Similar to TrendAnalysis but with covariates.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("response_bands", "Response variable bands (comma-separated)",
                "Time series of the response (dependent) variable"),
            ParameterDef::file("predictor_bands", "Predictor variable bands (comma-separated)",
                "Time series of predictor (independent) variables. Total count must be "
                "a multiple of the response band count."),
            ParameterDef::output("output_coefficients", "Output coefficients raster"),
            ParameterDef::output("output_rsquared", "Output R-squared raster"),
        };
    }

    bool execute() override {
        QStringList respFiles = parameter("response_bands").toString().split(",", Qt::SkipEmptyParts);
        QStringList predFiles = parameter("predictor_bands").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : respFiles) f = f.trimmed();
        for (auto& f : predFiles) f = f.trimmed();

        int nTime = respFiles.size();
        QString outCoeffs = parameter("output_coefficients").toString();
        QString outR2 = parameter("output_rsquared").toString();

        if (nTime < 2) {
            setError("At least 2 time steps required");
            return false;
        }
        if (predFiles.size() % nTime != 0) {
            setError(QString("Number of predictor bands (%1) must be a multiple of "
                             "response bands count (%2)")
                     .arg(predFiles.size()).arg(nTime));
            return false;
        }

        int nPred = predFiles.size() / nTime;
        // Model: y(t) = b0 + b1*t + b2*x1(t) + b3*x2(t) + ...
        // nCoeffs = 1 (intercept) + 1 (time) + nPred
        int nCoeffs = 2 + nPred;

        if (nTime < nCoeffs) {
            setError("More time steps required than model parameters");
            return false;
        }

        reportProgress(0.0, "Reading response bands...");

        std::vector<std::unique_ptr<Raster>> respRasters;
        respRasters.reserve(nTime);
        for (int t = 0; t < nTime; ++t) {
            auto r = GdalIO::read(respFiles[t]);
            if (!r) { setError(QString("Failed to read: %1").arg(respFiles[t])); return false; }
            if (t > 0 && (r->cols() != respRasters[0]->cols() || r->rows() != respRasters[0]->rows())) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
            respRasters.push_back(std::move(r));
        }

        reportProgress(0.1, "Reading predictor bands...");

        std::vector<std::vector<std::unique_ptr<Raster>>> predRasters(nPred);
        for (int p = 0; p < nPred; ++p) {
            predRasters[p].reserve(nTime);
            for (int t = 0; t < nTime; ++t) {
                int fileIdx = p * nTime + t;
                auto r = GdalIO::read(predFiles[fileIdx]);
                if (!r) { setError(QString("Failed to read: %1").arg(predFiles[fileIdx])); return false; }
                if (r->cols() != respRasters[0]->cols() || r->rows() != respRasters[0]->rows()) {
                    setError("All input rasters must have the same dimensions");
                    return false;
                }
                predRasters[p].push_back(std::move(r));
            }
        }

        int cols = respRasters[0]->cols(), rows = respRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = respRasters[0]->noDataValue();
        bool hasND = respRasters[0]->hasNoData();

        // Output: coefficients (nCoeffs bands) and R-squared (1 band)
        Raster coeffOut(cols, rows, nCoeffs, DataType::Float64);
        coeffOut.setGeoTransform(respRasters[0]->geoTransform());
        coeffOut.setProjection(respRasters[0]->projection());
        coeffOut.setNoDataValue(noData);

        Raster r2Out(cols, rows, 1, DataType::Float64);
        r2Out.setGeoTransform(respRasters[0]->geoTransform());
        r2Out.setProjection(respRasters[0]->projection());
        r2Out.setNoDataValue(noData);

        reportProgress(0.2, "Computing per-pixel multi-temporal regression...");

        for (int64_t i = 0; i < total; ++i) {
            bool valid = true;

            std::vector<double> y(nTime);
            for (int t = 0; t < nTime && valid; ++t) {
                double val = respRasters[t]->data(0)[i];
                if (hasND && val == noData) valid = false;
                y[t] = val;
            }

            // Design matrix: [1, t, x1(t), x2(t), ...]
            std::vector<std::vector<double>> X(nTime, std::vector<double>(nCoeffs));
            for (int t = 0; t < nTime && valid; ++t) {
                X[t][0] = 1.0;
                X[t][1] = t;
                for (int p = 0; p < nPred; ++p) {
                    double val = predRasters[p][t]->data(0)[i];
                    if (hasND && val == noData) valid = false;
                    X[t][2 + p] = val;
                }
            }

            if (!valid) {
                for (int c = 0; c < nCoeffs; ++c) coeffOut.data(c)[i] = noData;
                r2Out.data(0)[i] = noData;
                continue;
            }

            // Normal equations: X^T X beta = X^T y
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
            }

            r2Out.data(0)[i] = (ssTot > 1e-15) ? 1.0 - ssRes / ssTot : 0.0;
            for (int c = 0; c < nCoeffs; ++c) coeffOut.data(c)[i] = beta[c];

            if (i % 500000 == 0)
                reportProgress(0.2 + 0.7 * static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");
        if (!GdalIO::write(coeffOut, outCoeffs)) {
            setError("Failed to write coefficients output");
            return false;
        }
        if (!GdalIO::write(r2Out, outR2)) {
            setError("Failed to write R-squared output");
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

REGISTER_MODULE(MlrModule)

} // namespace aplaceholder
