#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class LogicXorModule : public Module {
public:
    QString name() const override { return "LOGICXOR"; }
    QString description() const override {
        return "Logical XOR overlay. Output is 1 where exactly one input is non-zero, else 0.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First input image"),
            ParameterDef::file("input2", "Second input image"),
            ParameterDef::output("output", "Output image"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input1").toString());
        auto r2 = GdalIO::read(parameter("input2").toString());
        if (!r1 || !r2) {
            setError("Failed to read input rasters");
            return false;
        }
        if (r1->cols() != r2->cols() || r1->rows() != r2->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(r1->noDataValue());

        const auto& d1 = r1->data(0);
        const auto& d2 = r2->data(0);
        auto& out = output.data(0);
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (d1[i] == noData || d2[i] == noData)) {
                out[i] = noData;
                continue;
            }
            bool a = (d1[i] != 0.0);
            bool b = (d2[i] != 0.0);
            out[i] = (a != b) ? 1.0 : 0.0;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(LogicXorModule)

} // namespace aplaceholder
