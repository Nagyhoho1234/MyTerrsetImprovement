#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class LogicNotModule : public Module {
public:
    QString name() const override { return "LOGICNOT"; }
    QString description() const override {
        return "Logical NOT. Output is 1 where input is zero, 0 where input is non-zero.";
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
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(r1->noDataValue());

        const auto& d1 = r1->data(0);
        auto& out = output.data(0);
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && d1[i] == noData) {
                out[i] = noData;
                continue;
            }
            out[i] = (d1[i] == 0.0) ? 1.0 : 0.0;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(LogicNotModule)

} // namespace aplaceholder
