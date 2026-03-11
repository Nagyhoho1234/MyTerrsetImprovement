#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>

namespace aplaceholder {

class TimeSeriesStatsModule : public Module {
public:
    QString name() const override { return "TIMESERIES_STATS"; }
    QString description() const override {
        return "Computes per-pixel temporal statistics across a time series of "
               "raster images. Supports mean, minimum, maximum, standard deviation, "
               "sum, and count operations.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_series", "Input time series (multi-file)"),
            ParameterDef::combo("statistic", "Statistic",
                {"Mean", "Min", "Max", "StdDev", "Sum", "Count"}, 0,
                "Temporal statistic to compute for each pixel"),
            ParameterDef::output("output", "Output statistics image"),
        };
    }

    bool execute() override {
        // Parse comma-separated input file list
        QString inputStr = parameter("input_series").toString();
        QStringList files = inputStr.split(",", Qt::SkipEmptyParts);
        for (auto& f : files) f = f.trimmed();

        int n = files.size();
        if (n < 1) {
            setError("At least 1 input raster is required");
            return false;
        }

        int statIndex = parameter("statistic").toInt();

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
            reportProgress(0.2 * (t + 1) / n);
        }

        int cols = rasters[0]->cols(), rows = rasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = rasters[0]->noDataValue();
        bool hasND = rasters[0]->hasNoData();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(rasters[0]->geoTransform());
        output.setProjection(rasters[0]->projection());
        output.setNoDataValue(noData);

        auto& dst = output.data(0);

        reportProgress(0.2, "Computing statistics...");

        for (int64_t i = 0; i < total; ++i) {
            // Gather valid values across time
            double sum = 0.0;
            double sumSq = 0.0;
            double minVal = std::numeric_limits<double>::max();
            double maxVal = -std::numeric_limits<double>::max();
            int validCount = 0;

            for (int t = 0; t < n; ++t) {
                double val = rasters[t]->data(0)[i];
                if (hasND && val == noData) continue;

                sum += val;
                sumSq += val * val;
                if (val < minVal) minVal = val;
                if (val > maxVal) maxVal = val;
                validCount++;
            }

            if (validCount == 0) {
                dst[i] = noData;
                continue;
            }

            switch (statIndex) {
                case 0: // Mean
                    dst[i] = sum / validCount;
                    break;
                case 1: // Min
                    dst[i] = minVal;
                    break;
                case 2: // Max
                    dst[i] = maxVal;
                    break;
                case 3: { // StdDev
                    if (validCount < 2) {
                        dst[i] = 0.0;
                    } else {
                        double mean = sum / validCount;
                        double variance = (sumSq - validCount * mean * mean) / (validCount - 1);
                        dst[i] = (variance > 0.0) ? std::sqrt(variance) : 0.0;
                    }
                    break;
                }
                case 4: // Sum
                    dst[i] = sum;
                    break;
                case 5: // Count of valid
                    dst[i] = static_cast<double>(validCount);
                    break;
            }

            if (i % 1000000 == 0)
                reportProgress(0.2 + 0.75 * static_cast<double>(i) / total);
        }

        reportProgress(0.95, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(TimeSeriesStatsModule)

} // namespace aplaceholder
