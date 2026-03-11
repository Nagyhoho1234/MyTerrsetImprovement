#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace aplaceholder {

class TAggregateModule : public Module {
public:
    QString name() const override { return "TAGGREGATE"; }
    QString description() const override {
        return "Temporal aggregation. Aggregates a time series to coarser temporal "
               "resolution (e.g., monthly to annual) using mean, sum, max, or min "
               "over groups of consecutive time steps.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series", "Time series rasters (comma-separated)",
                "Comma-separated list of time series raster file paths in chronological order"),
            ParameterDef::output("output_prefix", "Output filename prefix",
                "Output rasters named prefix_agg001, prefix_agg002, etc."),
            ParameterDef::integer("group_size", "Group size", 12, 2, 365,
                "Number of time steps per aggregation group"),
            ParameterDef::combo("method", "Aggregation method",
                {"mean", "sum", "max", "min"}, 0,
                "Statistical method for aggregation"),
        };
    }

    bool execute() override {
        QString tsParam = parameter("time_series").toString();
        QStringList tsPaths = tsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : tsPaths) p = p.trimmed();

        QString prefix = parameter("output_prefix").toString();
        int groupSize = parameter("group_size").toInt();
        int methodIdx = parameter("method").toInt();
        int numSteps = tsPaths.size();

        if (numSteps < groupSize) {
            setError("Time series must have at least one full group.");
            return false;
        }

        int numGroups = numSteps / groupSize;

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

        reportProgress(0.2, "Aggregating...");

        for (int g = 0; g < numGroups; ++g) {
            int tStart = g * groupSize;
            int tEnd = tStart + groupSize;

            Raster output(cols, rows, 1, DataType::Float32);
            output.setGeoTransform(rasters[0]->geoTransform());
            output.setProjection(rasters[0]->projection());
            output.setNoDataValue(outNoData);
            auto& outData = output.data(0);

            for (int64_t i = 0; i < total; ++i) {
                double result = 0.0;
                int count = 0;
                bool first = true;

                for (int t = tStart; t < tEnd; ++t) {
                    double val = (*data[t])[i];
                    if (hasND && val == noData) continue;

                    if (first) {
                        result = val;
                        first = false;
                    } else {
                        switch (methodIdx) {
                        case 0: // mean
                        case 1: // sum
                            result += val;
                            break;
                        case 2: // max
                            if (val > result) result = val;
                            break;
                        case 3: // min
                            if (val < result) result = val;
                            break;
                        }
                    }
                    ++count;
                }

                if (count == 0) {
                    outData[i] = outNoData;
                } else {
                    if (methodIdx == 0) // mean
                        result /= count;
                    outData[i] = result;
                }
            }

            QString outPath = prefix + "_agg" + QString::number(g + 1).rightJustified(3, '0');
            if (!GdalIO::write(output, outPath)) {
                setError("Failed to write: " + outPath);
                return false;
            }

            reportProgress(0.2 + 0.75 * static_cast<double>(g + 1) / numGroups);
        }

        reportProgress(1.0, "Temporal aggregation complete.");
        return true;
    }
};

REGISTER_MODULE(TAggregateModule)

} // namespace aplaceholder
