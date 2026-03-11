#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class HarmonicRegModule : public Module {
public:
    QString name() const override { return "HARMONICREG"; }
    QString description() const override {
        return "Harmonic regression. Fits sin/cos harmonics to a time series using "
               "ordinary least squares. Outputs the intercept, amplitude and phase "
               "for each harmonic, and the fitted values. Useful for characterizing "
               "seasonal patterns and phenological cycles.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series", "Time series rasters (comma-separated)",
                "Comma-separated list of time series raster file paths in chronological order"),
            ParameterDef::output("output_prefix", "Output filename prefix",
                "Output: prefix_intercept, prefix_amp1, prefix_phase1, ..., prefix_fitted001, ..."),
            ParameterDef::integer("num_harmonics", "Number of harmonics", 2, 1, 10,
                "Number of sin/cos harmonic pairs to fit"),
        };
    }

    bool execute() override {
        QString tsParam = parameter("time_series").toString();
        QStringList tsPaths = tsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : tsPaths) p = p.trimmed();

        QString prefix = parameter("output_prefix").toString();
        int numHarmonics = parameter("num_harmonics").toInt();
        int N = tsPaths.size();

        int numCoeffs = 1 + 2 * numHarmonics;  // intercept + sin/cos pairs
        if (N < numCoeffs) {
            setError(QString("Need at least %1 time steps for %2 harmonics (got %3).")
                     .arg(numCoeffs).arg(numHarmonics).arg(N));
            return false;
        }

        // Read all time steps
        reportProgress(0.0, "Reading time series...");
        std::vector<std::unique_ptr<Raster>> rasters(N);
        int cols = 0, rows = 0;

        for (int t = 0; t < N; ++t) {
            rasters[t] = GdalIO::read(tsPaths[t]);
            if (!rasters[t]) {
                setError("Failed to read: " + tsPaths[t]);
                return false;
            }
            if (t == 0) {
                cols = rasters[0]->cols();
                rows = rasters[0]->rows();
            } else if (rasters[t]->cols() != cols || rasters[t]->rows() != rows) {
                setError("All time series rasters must have the same dimensions.");
                return false;
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = rasters[0]->hasNoData();
        double noData = rasters[0]->noDataValue();
        double outNoData = -9999.0;

        std::vector<const std::vector<double>*> data(N);
        for (int t = 0; t < N; ++t)
            data[t] = &rasters[t]->data(0);

        // Build design matrix X: [1, cos(2*pi*t/N), sin(2*pi*t/N), cos(4*pi*t/N), ...]
        // X is N x numCoeffs
        std::vector<std::vector<double>> X(N, std::vector<double>(numCoeffs));
        for (int t = 0; t < N; ++t) {
            X[t][0] = 1.0;
            for (int h = 0; h < numHarmonics; ++h) {
                double angle = 2.0 * M_PI * (h + 1) * t / N;
                X[t][1 + 2 * h] = std::cos(angle);
                X[t][2 + 2 * h] = std::sin(angle);
            }
        }

        // Compute X'X and its inverse (small matrix, numCoeffs x numCoeffs)
        reportProgress(0.1, "Computing design matrix...");
        std::vector<double> XtX(numCoeffs * numCoeffs, 0.0);
        for (int i = 0; i < numCoeffs; ++i) {
            for (int j = i; j < numCoeffs; ++j) {
                double sum = 0.0;
                for (int t = 0; t < N; ++t)
                    sum += X[t][i] * X[t][j];
                XtX[i * numCoeffs + j] = sum;
                XtX[j * numCoeffs + i] = sum;
            }
        }

        // Invert X'X using Gauss-Jordan elimination
        std::vector<double> inv(numCoeffs * numCoeffs, 0.0);
        for (int i = 0; i < numCoeffs; ++i)
            inv[i * numCoeffs + i] = 1.0;

        std::vector<double> aug = XtX;
        for (int i = 0; i < numCoeffs; ++i) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k < numCoeffs; ++k) {
                if (std::abs(aug[k * numCoeffs + i]) > std::abs(aug[maxRow * numCoeffs + i]))
                    maxRow = k;
            }
            // Swap rows
            if (maxRow != i) {
                for (int j = 0; j < numCoeffs; ++j) {
                    std::swap(aug[i * numCoeffs + j], aug[maxRow * numCoeffs + j]);
                    std::swap(inv[i * numCoeffs + j], inv[maxRow * numCoeffs + j]);
                }
            }

            double pivot = aug[i * numCoeffs + i];
            if (std::abs(pivot) < 1e-15) {
                setError("Design matrix is singular. Check input time series.");
                return false;
            }

            for (int j = 0; j < numCoeffs; ++j) {
                aug[i * numCoeffs + j] /= pivot;
                inv[i * numCoeffs + j] /= pivot;
            }

            for (int k = 0; k < numCoeffs; ++k) {
                if (k == i) continue;
                double factor = aug[k * numCoeffs + i];
                for (int j = 0; j < numCoeffs; ++j) {
                    aug[k * numCoeffs + j] -= factor * aug[i * numCoeffs + j];
                    inv[k * numCoeffs + j] -= factor * inv[i * numCoeffs + j];
                }
            }
        }

        // Precompute (X'X)^-1 X' = invXtX * X'
        // This is numCoeffs x N
        std::vector<std::vector<double>> invXtXXt(numCoeffs, std::vector<double>(N, 0.0));
        for (int i = 0; i < numCoeffs; ++i) {
            for (int t = 0; t < N; ++t) {
                double sum = 0.0;
                for (int j = 0; j < numCoeffs; ++j)
                    sum += inv[i * numCoeffs + j] * X[t][j];
                invXtXXt[i][t] = sum;
            }
        }

        // Create output rasters for coefficients
        // intercept, amp1, phase1, amp2, phase2, ...
        int numOutputCoeffs = 1 + 2 * numHarmonics;  // intercept + amplitude + phase per harmonic
        std::vector<Raster> coeffOutputs;
        coeffOutputs.reserve(numOutputCoeffs);
        for (int c = 0; c < numOutputCoeffs; ++c) {
            coeffOutputs.emplace_back(cols, rows, 1, DataType::Float32);
            coeffOutputs[c].setGeoTransform(rasters[0]->geoTransform());
            coeffOutputs[c].setProjection(rasters[0]->projection());
            coeffOutputs[c].setNoDataValue(outNoData);
        }

        reportProgress(0.2, "Fitting harmonics...");

        for (int64_t i = 0; i < total; ++i) {
            // Check nodata
            bool valid = true;
            if (hasND) {
                for (int t = 0; t < N && valid; ++t)
                    if ((*data[t])[i] == noData) valid = false;
            }

            if (!valid) {
                for (int c = 0; c < numOutputCoeffs; ++c)
                    coeffOutputs[c].data(0)[i] = outNoData;
                continue;
            }

            // Compute beta = (X'X)^-1 X' y
            std::vector<double> beta(numCoeffs, 0.0);
            for (int j = 0; j < numCoeffs; ++j) {
                double sum = 0.0;
                for (int t = 0; t < N; ++t)
                    sum += invXtXXt[j][t] * (*data[t])[i];
                beta[j] = sum;
            }

            // Store intercept
            coeffOutputs[0].data(0)[i] = beta[0];

            // Compute amplitude and phase for each harmonic
            for (int h = 0; h < numHarmonics; ++h) {
                double cosCoeff = beta[1 + 2 * h];
                double sinCoeff = beta[2 + 2 * h];
                double amplitude = std::sqrt(cosCoeff * cosCoeff + sinCoeff * sinCoeff);
                double phase = std::atan2(sinCoeff, cosCoeff);

                coeffOutputs[1 + 2 * h].data(0)[i] = amplitude;
                coeffOutputs[2 + 2 * h].data(0)[i] = phase;
            }

            if (i % 1000000 == 0)
                reportProgress(0.2 + 0.5 * static_cast<double>(i) / total);
        }

        // Write coefficient outputs
        reportProgress(0.7, "Writing coefficient rasters...");
        if (!GdalIO::write(coeffOutputs[0], prefix + "_intercept")) {
            setError("Failed to write intercept raster.");
            return false;
        }
        for (int h = 0; h < numHarmonics; ++h) {
            QString hStr = QString::number(h + 1);
            if (!GdalIO::write(coeffOutputs[1 + 2 * h], prefix + "_amp" + hStr)) {
                setError("Failed to write amplitude raster for harmonic " + hStr);
                return false;
            }
            if (!GdalIO::write(coeffOutputs[2 + 2 * h], prefix + "_phase" + hStr)) {
                setError("Failed to write phase raster for harmonic " + hStr);
                return false;
            }
        }

        // Compute and write fitted values
        reportProgress(0.8, "Computing fitted values...");
        for (int t = 0; t < N; ++t) {
            Raster fitted(cols, rows, 1, DataType::Float32);
            fitted.setGeoTransform(rasters[0]->geoTransform());
            fitted.setProjection(rasters[0]->projection());
            fitted.setNoDataValue(outNoData);
            auto& fData = fitted.data(0);

            for (int64_t i = 0; i < total; ++i) {
                if (coeffOutputs[0].data(0)[i] == outNoData) {
                    fData[i] = outNoData;
                    continue;
                }

                double val = coeffOutputs[0].data(0)[i];  // intercept
                for (int h = 0; h < numHarmonics; ++h) {
                    double amp = coeffOutputs[1 + 2 * h].data(0)[i];
                    double phase = coeffOutputs[2 + 2 * h].data(0)[i];
                    double angle = 2.0 * M_PI * (h + 1) * t / N;
                    val += amp * std::cos(angle - phase);
                }
                fData[i] = val;
            }

            QString outPath = prefix + "_fitted" + QString::number(t + 1).rightJustified(3, '0');
            if (!GdalIO::write(fitted, outPath)) {
                setError("Failed to write fitted raster: " + outPath);
                return false;
            }

            reportProgress(0.8 + 0.18 * static_cast<double>(t + 1) / N);
        }

        reportProgress(1.0, "Harmonic regression complete.");
        return true;
    }
};

REGISTER_MODULE(HarmonicRegModule)

} // namespace aplaceholder
