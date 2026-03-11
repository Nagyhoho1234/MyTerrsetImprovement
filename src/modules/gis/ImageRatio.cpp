#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ImageRatioModule : public Module {
public:
    QString name() const override { return "IMAGERATIO"; }
    QString description() const override {
        return "Image ratioing. Divides one raster by another for change detection "
               "or spectral band ratioing. Pixels where the divisor is zero are "
               "assigned NoData.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "Numerator raster"),
            ParameterDef::file("input2", "Denominator raster"),
            ParameterDef::output("output", "Output ratio raster"),
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
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& d1 = r1->data(0);
        const auto& d2 = r2->data(0);
        double noData1 = r1->noDataValue();
        double noData2 = r2->noDataValue();
        bool hasND1 = r1->hasNoData();
        bool hasND2 = r2->hasNoData();

        reportProgress(0.0, "Computing image ratio...");

        double outNoData = -9999.0;
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(outNoData);
        auto& out = output.data(0);

        int64_t divByZeroCount = 0;
        for (int64_t i = 0; i < total; ++i) {
            if ((hasND1 && d1[i] == noData1) || (hasND2 && d2[i] == noData2)) {
                out[i] = outNoData;
            } else if (std::abs(d2[i]) < 1e-15) {
                out[i] = outNoData;
                ++divByZeroCount;
            } else {
                out[i] = d1[i] / d2[i];
            }

            if (i % 1000000 == 0)
                reportProgress(0.8 * static_cast<double>(i) / total);
        }

        reportProgress(0.85, "Writing output raster...");

        if (!GdalIO::write(output, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0,
            QString("Image ratio complete: %1 pixels, %2 division-by-zero pixels set to NoData")
                .arg(total)
                .arg(divByZeroCount));

        return true;
    }
};

REGISTER_MODULE(ImageRatioModule)

} // namespace aplaceholder
