#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>

namespace aplaceholder {

class EotModule : public Module {
public:
    QString name() const override { return "EOT"; }
    QString description() const override {
        return "Empirical Orthogonal Teleconnections. Identifies spatial patterns in "
               "time series that explain variance at remote locations.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series_bands", "Time series bands (comma-separated)",
                "Input raster bands representing time steps"),
            ParameterDef::output("output_prefix", "Output prefix"),
            ParameterDef::integer("num_modes", "Number of EOT modes", 5, 1, 50,
                "Number of EOT modes to extract"),
        };
    }

    bool execute() override {
        QStringList files = parameter("time_series_bands").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : files) f = f.trimmed();
        int n = files.size();
        int numModes = parameter("num_modes").toInt();
        QString outPrefix = parameter("output_prefix").toString();

        if (n < 3) {
            setError("At least 3 time steps are required");
            return false;
        }

        reportProgress(0.0, "Reading time series...");

        std::vector<std::unique_ptr<Raster>> rasters;
        rasters.reserve(n);
        for (int t = 0; t < n; ++t) {
            auto r = GdalIO::read(files[t]);
            if (!r) { setError(QString("Failed to read: %1").arg(files[t])); return false; }
            if (t > 0 && (r->cols() != rasters[0]->cols() || r->rows() != rasters[0]->rows())) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
            rasters.push_back(std::move(r));
        }

        int cols = rasters[0]->cols(), rows = rasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = rasters[0]->noDataValue();
        bool hasND = rasters[0]->hasNoData();

        // Build matrix: each pixel has a time series of length n
        // Identify valid pixels
        std::vector<bool> validMask(total, true);
        for (int64_t i = 0; i < total; ++i) {
            for (int t = 0; t < n; ++t) {
                double val = rasters[t]->data(0)[i];
                if (hasND && val == noData) { validMask[i] = false; break; }
            }
        }

        // Collect valid pixel indices
        std::vector<int64_t> validIdx;
        validIdx.reserve(total);
        for (int64_t i = 0; i < total; ++i)
            if (validMask[i]) validIdx.push_back(i);

        int64_t nValid = validIdx.size();
        if (nValid < 2) {
            setError("Insufficient valid pixels");
            return false;
        }
        numModes = std::min(numModes, (int)std::min((int64_t)n, nValid));

        reportProgress(0.1, "Building pixel time series matrix...");

        // Build time series matrix: timeSeries[pixel_idx][time]
        // and compute per-pixel means for demeaning
        std::vector<std::vector<double>> ts(nValid, std::vector<double>(n));
        std::vector<double> pixelMeans(nValid, 0.0);
        for (int64_t pi = 0; pi < nValid; ++pi) {
            int64_t idx = validIdx[pi];
            double sum = 0.0;
            for (int t = 0; t < n; ++t) {
                ts[pi][t] = rasters[t]->data(0)[idx];
                sum += ts[pi][t];
            }
            pixelMeans[pi] = sum / n;
            for (int t = 0; t < n; ++t) ts[pi][t] -= pixelMeans[pi];
        }

        // Residual time series (modified during deflation)
        auto residual = ts;

        reportProgress(0.2, "Computing EOT modes...");

        for (int mode = 0; mode < numModes; ++mode) {
            reportProgress(0.2 + 0.7 * mode / numModes,
                           QString("Computing EOT mode %1...").arg(mode + 1));

            // For each candidate base point, compute total explained variance
            // at all other points via regression
            double bestTotalR2 = -1.0;
            int64_t bestBase = 0;

            // For efficiency, subsample if too many pixels
            int64_t step = std::max((int64_t)1, nValid / 2000);
            for (int64_t pi = 0; pi < nValid; pi += step) {
                // Compute sum of squared R^2 with all other points
                // Base point time series
                const auto& baseSeries = residual[pi];
                double baseSS = 0.0;
                for (int t = 0; t < n; ++t) baseSS += baseSeries[t] * baseSeries[t];
                if (baseSS < 1e-15) continue;

                double totalR2 = 0.0;
                for (int64_t pj = 0; pj < nValid; pj += step) {
                    double crossSum = 0.0, targetSS = 0.0;
                    for (int t = 0; t < n; ++t) {
                        crossSum += baseSeries[t] * residual[pj][t];
                        targetSS += residual[pj][t] * residual[pj][t];
                    }
                    if (targetSS < 1e-15) continue;
                    double r2 = (crossSum * crossSum) / (baseSS * targetSS);
                    totalR2 += r2;
                }

                if (totalR2 > bestTotalR2) {
                    bestTotalR2 = totalR2;
                    bestBase = pi;
                }
            }

            // Compute regression coefficients from best base point to all pixels
            const auto& baseSeries = residual[bestBase];
            double baseSS = 0.0;
            for (int t = 0; t < n; ++t) baseSS += baseSeries[t] * baseSeries[t];

            // Output: regression coefficient map (spatial pattern) and base point time series
            Raster patternOut(cols, rows, 1, DataType::Float64);
            patternOut.setGeoTransform(rasters[0]->geoTransform());
            patternOut.setProjection(rasters[0]->projection());
            patternOut.setNoDataValue(noData);
            auto& patData = patternOut.data(0);
            std::fill(patData.begin(), patData.end(), noData);

            Raster r2Out(cols, rows, 1, DataType::Float64);
            r2Out.setGeoTransform(rasters[0]->geoTransform());
            r2Out.setProjection(rasters[0]->projection());
            r2Out.setNoDataValue(noData);
            auto& r2Data = r2Out.data(0);
            std::fill(r2Data.begin(), r2Data.end(), noData);

            for (int64_t pi = 0; pi < nValid; ++pi) {
                int64_t idx = validIdx[pi];
                double crossSum = 0.0, targetSS = 0.0;
                for (int t = 0; t < n; ++t) {
                    crossSum += baseSeries[t] * residual[pi][t];
                    targetSS += residual[pi][t] * residual[pi][t];
                }
                double slope = (baseSS > 1e-15) ? crossSum / baseSS : 0.0;
                double r2 = (targetSS > 1e-15) ? (crossSum * crossSum) / (baseSS * targetSS) : 0.0;
                patData[idx] = slope;
                r2Data[idx] = r2;

                // Deflate: remove explained variance
                for (int t = 0; t < n; ++t)
                    residual[pi][t] -= slope * baseSeries[t];
            }

            // Write outputs
            QString patPath = QString("%1_pattern_%2.tif").arg(outPrefix).arg(mode + 1);
            QString r2Path = QString("%1_r2_%2.tif").arg(outPrefix).arg(mode + 1);
            if (!GdalIO::write(patternOut, patPath)) {
                setError(QString("Failed to write: %1").arg(patPath));
                return false;
            }
            if (!GdalIO::write(r2Out, r2Path)) {
                setError(QString("Failed to write: %1").arg(r2Path));
                return false;
            }
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(EotModule)

} // namespace aplaceholder
