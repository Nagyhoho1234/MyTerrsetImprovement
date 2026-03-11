#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class HighlightModule : public Module {
public:
    QString name() const override { return "HIGHLIGHT"; }
    QString description() const override {
        return "Highlight specific value ranges. Outputs binary raster: "
               "1 where input is within [min_val, max_val], 0 elsewhere.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::real("min_val", "Minimum value", 0.0, -1e15, 1e15,
                "Lower bound of highlight range (inclusive)"),
            ParameterDef::real("max_val", "Maximum value", 1.0, -1e15, 1e15,
                "Upper bound of highlight range (inclusive)"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        double minVal = parameter("min_val").toDouble();
        double maxVal = parameter("max_val").toDouble();

        int cols = raster->cols(), rows = raster->rows();
        Raster output(cols, rows, 1, DataType::Byte);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());

        const auto& src = raster->data(0);
        auto& dst = output.data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && src[i] == noData) {
                dst[i] = 0.0;
            } else {
                dst[i] = (src[i] >= minVal && src[i] <= maxVal) ? 1.0 : 0.0;
            }

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(HighlightModule)

} // namespace aplaceholder
