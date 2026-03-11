#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class DeseasonModule : public Module {
public:
    QString name() const override { return "DESEASON"; }
    QString description() const override {
        return "Remove the seasonal component from a time series. Computes the mean "
               "for each position within the seasonal cycle and subtracts it to "
               "produce anomaly values.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series", "Time series rasters (comma-separated)",
                "Comma-separated list of time series raster file paths in chronological order"),
            ParameterDef::output("output_prefix", "Output filename prefix",
                "Output anomaly rasters will be named prefix_anom001, prefix_anom002, etc."),
            ParameterDef::integer("period", "Seasonal period", 12, 2, 365,
                "Number of time steps per seasonal cycle (e.g., 12 for monthly data)"),
        };
    }

    bool execute() override {
        QString tsParam = parameter("time_series").toString();
        QStringList tsPaths = tsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : tsPaths) p = p.trimmed();

        QString prefix = parameter("output_prefix").toString();
        int period = parameter("period").toInt();
        int numSteps = tsPaths.size();

        if (numSteps < period) {
            setError("Time series must have at least one full seasonal period.");
            return false;
        }

        // Read all time steps
        reportProgress(0.0, "Reading time series...");
        std::vector<std::unique_ptr<Raster>> rasters(numSteps);
        int cols = 0, rows = 0;

        for (int t = 0; t < numSteps; ++t) {
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

        // Collect data pointers
        std::vector<const std::vector<double>*> data(numSteps);
        for (int t = 0; t < numSteps; ++t)
            data[t] = &rasters[t]->data(0);

        // Compute seasonal means per pixel
        reportProgress(0.2, "Computing seasonal means...");
        // seasonalMean[phase][pixel]
        std::vector<std::vector<double>> seasonalSum(period, std::vector<double>(total, 0.0));
        std::vector<std::vector<int>> seasonalCount(period, std::vector<int>(total, 0));

        for (int t = 0; t < numSteps; ++t) {
            int phase = t % period;
            for (int64_t i = 0; i < total; ++i) {
                double val = (*data[t])[i];
                if (hasND && val == noData) continue;
                seasonalSum[phase][i] += val;
                seasonalCount[phase][i] += 1;
            }
        }

        // Compute seasonal means
        std::vector<std::vector<double>> seasonalMean(period, std::vector<double>(total, 0.0));
        for (int phase = 0; phase < period; ++phase) {
            for (int64_t i = 0; i < total; ++i) {
                if (seasonalCount[phase][i] > 0)
                    seasonalMean[phase][i] = seasonalSum[phase][i] / seasonalCount[phase][i];
            }
        }

        // Subtract seasonal mean and write output
        reportProgress(0.5, "Computing anomalies...");
        for (int t = 0; t < numSteps; ++t) {
            int phase = t % period;

            Raster output(cols, rows, 1, DataType::Float32);
            output.setGeoTransform(rasters[0]->geoTransform());
            output.setProjection(rasters[0]->projection());
            output.setNoDataValue(outNoData);
            auto& outData = output.data(0);

            for (int64_t i = 0; i < total; ++i) {
                double val = (*data[t])[i];
                if (hasND && val == noData) {
                    outData[i] = outNoData;
                } else {
                    outData[i] = val - seasonalMean[phase][i];
                }
            }

            QString outPath = prefix + "_anom" + QString::number(t + 1).rightJustified(3, '0');
            if (!GdalIO::write(output, outPath)) {
                setError("Failed to write: " + outPath);
                return false;
            }

            reportProgress(0.5 + 0.45 * static_cast<double>(t + 1) / numSteps);
        }

        reportProgress(1.0, "Deseasoning complete.");
        return true;
    }
};

REGISTER_MODULE(DeseasonModule)

} // namespace aplaceholder
