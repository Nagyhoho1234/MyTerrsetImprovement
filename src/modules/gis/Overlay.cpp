#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class OverlayModule : public Module {
public:
    QString name() const override { return "OVERLAY"; }
    QString description() const override {
        return "Mathematical overlay of two raster images. "
               "Performs pixel-by-pixel arithmetic operations.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First image"),
            ParameterDef::file("input2", "Second image"),
            ParameterDef::output("output", "Output image"),
            ParameterDef::combo("operation", "Operation",
                {"Add", "Subtract", "Multiply", "Divide", "Maximum", "Minimum",
                 "Average", "Cover"}, 0,
                "Mathematical operation to perform"),
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

        int op = parameter("operation").toInt();
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

            switch (op) {
                case 0: out[i] = d1[i] + d2[i]; break;
                case 1: out[i] = d1[i] - d2[i]; break;
                case 2: out[i] = d1[i] * d2[i]; break;
                case 3: out[i] = (d2[i] != 0) ? d1[i] / d2[i] : noData; break;
                case 4: out[i] = std::max(d1[i], d2[i]); break;
                case 5: out[i] = std::min(d1[i], d2[i]); break;
                case 6: out[i] = (d1[i] + d2[i]) / 2.0; break;
                case 7: out[i] = (d2[i] != noData && d2[i] != 0) ? d2[i] : d1[i]; break;
            }

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(OverlayModule)

} // namespace aplaceholder
