#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <limits>

namespace aplaceholder {

class NormalizeModule : public Module {
public:
    QString name() const override { return "NORMALIZE"; }
    QString description() const override {
        return "Normalize raster values to the [0,1] range using min-max scaling.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output image"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input").toString());
        if (!r1) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        const auto& d1 = r1->data(0);
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Pass 1: find min and max
        double minVal = std::numeric_limits<double>::max();
        double maxVal = std::numeric_limits<double>::lowest();
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && d1[i] == noData) continue;
            if (d1[i] < minVal) minVal = d1[i];
            if (d1[i] > maxVal) maxVal = d1[i];
        }

        double range = maxVal - minVal;
        if (range == 0.0) range = 1.0; // avoid division by zero

        // Pass 2: scale to [0,1]
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(noData);

        auto& out = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && d1[i] == noData) {
                out[i] = noData;
                continue;
            }
            out[i] = (d1[i] - minVal) / range;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(NormalizeModule)

} // namespace aplaceholder
