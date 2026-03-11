#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class TDiffModule : public Module {
public:
    QString name() const override { return "TDIFF"; }
    QString description() const override {
        return "Temporal differencing. Computes the difference between consecutive "
               "time steps in a time series, producing N-1 difference images from "
               "N input images.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series", "Time series rasters (comma-separated)",
                "Comma-separated list of time series raster file paths in chronological order"),
            ParameterDef::output("output_prefix", "Output filename prefix",
                "Output rasters named prefix_diff001, prefix_diff002, etc."),
        };
    }

    bool execute() override {
        QString tsParam = parameter("time_series").toString();
        QStringList tsPaths = tsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : tsPaths) p = p.trimmed();

        QString prefix = parameter("output_prefix").toString();
        int numSteps = tsPaths.size();

        if (numSteps < 2) {
            setError("At least 2 time steps are required for differencing.");
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
        int numDiffs = numSteps - 1;

        reportProgress(0.2, "Computing temporal differences...");

        for (int t = 0; t < numDiffs; ++t) {
            const auto& data1 = rasters[t]->data(0);
            const auto& data2 = rasters[t + 1]->data(0);

            Raster output(cols, rows, 1, DataType::Float32);
            output.setGeoTransform(rasters[0]->geoTransform());
            output.setProjection(rasters[0]->projection());
            output.setNoDataValue(outNoData);
            auto& outData = output.data(0);

            for (int64_t i = 0; i < total; ++i) {
                double v1 = data1[i];
                double v2 = data2[i];
                if (hasND && (v1 == noData || v2 == noData)) {
                    outData[i] = outNoData;
                } else {
                    outData[i] = v2 - v1;
                }
            }

            QString outPath = prefix + "_diff" + QString::number(t + 1).rightJustified(3, '0');
            if (!GdalIO::write(output, outPath)) {
                setError("Failed to write: " + outPath);
                return false;
            }

            reportProgress(0.2 + 0.75 * static_cast<double>(t + 1) / numDiffs);
        }

        reportProgress(1.0, "Temporal differencing complete.");
        return true;
    }
};

REGISTER_MODULE(TDiffModule)

} // namespace aplaceholder
