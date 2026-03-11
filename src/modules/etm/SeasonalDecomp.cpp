#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class SeasonalDecompModule : public Module {
public:
    QString name() const override { return "SEASONAL_DECOMP"; }
    QString description() const override {
        return "Seasonal trend decomposition of a raster time series. "
               "Separates the signal into trend, seasonal, and residual components.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_series", "Input time series (multi-file)"),
            ParameterDef::integer("frequency", "Seasonal frequency", 12, 2, 365,
                "Number of observations per seasonal cycle (e.g. 12 for monthly)"),
            ParameterDef::output("output_trend", "Output trend component"),
            ParameterDef::output("output_seasonal", "Output seasonal component"),
            ParameterDef::output("output_residual", "Output residual component"),
        };
    }

    bool execute() override {
        // Parse comma-separated input file list
        QString inputStr = parameter("input_series").toString();
        QStringList files = inputStr.split(",", Qt::SkipEmptyParts);
        for (auto& f : files) f = f.trimmed();

        int n = files.size();
        int freq = parameter("frequency").toInt();

        if (n < freq * 2) {
            setError(QString("At least %1 time steps required (2 full cycles of frequency %2)")
                     .arg(freq * 2).arg(freq));
            return false;
        }

        reportProgress(0.0, "Reading input series...");

        // Read all rasters
        std::vector<std::unique_ptr<Raster>> rasters;
        rasters.reserve(n);
        for (int t = 0; t < n; ++t) {
            auto r = GdalIO::read(files[t]);
            if (!r) {
                setError(QString("Failed to read raster: %1").arg(files[t]));
                return false;
            }
            if (t > 0 && (r->cols() != rasters[0]->cols() || r->rows() != rasters[0]->rows())) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
            rasters.push_back(std::move(r));
            reportProgress(0.1 * (t + 1) / n);
        }

        int cols = rasters[0]->cols(), rows = rasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = rasters[0]->noDataValue();
        bool hasND = rasters[0]->hasNoData();

        // Output rasters: each is a multi-band raster with n bands (one per time step)
        Raster trendOut(cols, rows, n, DataType::Float64);
        trendOut.setGeoTransform(rasters[0]->geoTransform());
        trendOut.setProjection(rasters[0]->projection());
        trendOut.setNoDataValue(noData);

        Raster seasonalOut(cols, rows, n, DataType::Float64);
        seasonalOut.setGeoTransform(rasters[0]->geoTransform());
        seasonalOut.setProjection(rasters[0]->projection());
        seasonalOut.setNoDataValue(noData);

        Raster residualOut(cols, rows, n, DataType::Float64);
        residualOut.setGeoTransform(rasters[0]->geoTransform());
        residualOut.setProjection(rasters[0]->projection());
        residualOut.setNoDataValue(noData);

        reportProgress(0.1, "Decomposing time series...");

        // Half-window for centered moving average
        int halfWin = freq / 2;

        for (int64_t i = 0; i < total; ++i) {
            // Gather pixel time series
            std::vector<double> series(n);
            bool valid = true;

            for (int t = 0; t < n; ++t) {
                double val = rasters[t]->data(0)[i];
                if (hasND && val == noData) {
                    valid = false;
                    break;
                }
                series[t] = val;
            }

            if (!valid) {
                for (int t = 0; t < n; ++t) {
                    trendOut.data(t)[i] = noData;
                    seasonalOut.data(t)[i] = noData;
                    residualOut.data(t)[i] = noData;
                }
                continue;
            }

            // Step 1: Compute trend using centered moving average of period length
            std::vector<double> trend(n, 0.0);
            std::vector<bool> trendValid(n, false);

            for (int t = 0; t < n; ++t) {
                if (t - halfWin < 0 || t + halfWin >= n) continue;

                double sum = 0.0;
                int count = 0;

                if (freq % 2 == 0) {
                    // Even period: use weighted average (endpoints get half weight)
                    sum += 0.5 * series[t - halfWin];
                    for (int k = t - halfWin + 1; k < t + halfWin; ++k) {
                        sum += series[k];
                    }
                    sum += 0.5 * series[t + halfWin];
                    trend[t] = sum / freq;
                } else {
                    // Odd period: simple average
                    for (int k = t - halfWin; k <= t + halfWin; ++k) {
                        sum += series[k];
                        count++;
                    }
                    trend[t] = sum / count;
                }
                trendValid[t] = true;
            }

            // Step 2: Detrend
            std::vector<double> detrended(n, 0.0);
            for (int t = 0; t < n; ++t) {
                if (trendValid[t])
                    detrended[t] = series[t] - trend[t];
            }

            // Step 3: Compute seasonal component as mean of detrended values
            // for each position in the cycle
            std::vector<double> seasonalAvg(freq, 0.0);
            std::vector<int> seasonalCount(freq, 0);
            for (int t = 0; t < n; ++t) {
                if (!trendValid[t]) continue;
                int pos = t % freq;
                seasonalAvg[pos] += detrended[t];
                seasonalCount[pos]++;
            }
            for (int p = 0; p < freq; ++p) {
                if (seasonalCount[p] > 0)
                    seasonalAvg[p] /= seasonalCount[p];
            }

            // Center the seasonal component (subtract its mean so it sums to ~0)
            double seasonalMean = 0.0;
            for (int p = 0; p < freq; ++p)
                seasonalMean += seasonalAvg[p];
            seasonalMean /= freq;
            for (int p = 0; p < freq; ++p)
                seasonalAvg[p] -= seasonalMean;

            // Step 4: Assign outputs
            for (int t = 0; t < n; ++t) {
                double seasonal = seasonalAvg[t % freq];

                if (trendValid[t]) {
                    trendOut.data(t)[i] = trend[t];
                    seasonalOut.data(t)[i] = seasonal;
                    residualOut.data(t)[i] = series[t] - trend[t] - seasonal;
                } else {
                    // For endpoints where trend is not computable, use nearest valid trend
                    // or mark as nodata
                    trendOut.data(t)[i] = noData;
                    seasonalOut.data(t)[i] = seasonal;
                    residualOut.data(t)[i] = noData;
                }
            }

            if (i % 500000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");
        if (!GdalIO::write(trendOut, parameter("output_trend").toString())) {
            setError("Failed to write trend output");
            return false;
        }
        if (!GdalIO::write(seasonalOut, parameter("output_seasonal").toString())) {
            setError("Failed to write seasonal output");
            return false;
        }
        if (!GdalIO::write(residualOut, parameter("output_residual").toString())) {
            setError("Failed to write residual output");
            return false;
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(SeasonalDecompModule)

} // namespace aplaceholder
