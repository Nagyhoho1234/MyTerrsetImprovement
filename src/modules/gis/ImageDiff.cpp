#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ImageDiffModule : public Module {
public:
    QString name() const override { return "IMAGEDIFF"; }
    QString description() const override {
        return "Image differencing for change detection. Performs pixel-by-pixel "
               "subtraction of two rasters with optional thresholding. "
               "When threshold > 0, output is binary (1 = change, 0 = no change).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First input raster (earlier time)"),
            ParameterDef::file("input2", "Second input raster (later time)"),
            ParameterDef::output("output", "Output difference raster"),
            ParameterDef::real("threshold", "Change threshold", 0.0, 0.0, 999999.0,
                "If > 0, output is binary: 1 where |diff| > threshold, 0 otherwise"),
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
        double threshold = parameter("threshold").toDouble();

        reportProgress(0.0, "Computing image difference...");

        double outNoData = -9999.0;
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(outNoData);
        auto& out = output.data(0);

        int64_t changeCount = 0;
        for (int64_t i = 0; i < total; ++i) {
            if ((hasND1 && d1[i] == noData1) || (hasND2 && d2[i] == noData2)) {
                out[i] = outNoData;
            } else {
                double diff = d2[i] - d1[i];
                if (threshold > 0.0) {
                    out[i] = (std::abs(diff) > threshold) ? 1.0 : 0.0;
                    if (out[i] == 1.0) ++changeCount;
                } else {
                    out[i] = diff;
                }
            }

            if (i % 1000000 == 0)
                reportProgress(0.8 * static_cast<double>(i) / total);
        }

        reportProgress(0.85, "Writing output raster...");

        if (!GdalIO::write(output, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        QString msg;
        if (threshold > 0.0) {
            msg = QString("Image difference (threshold=%1): %2 change pixels of %3 total")
                      .arg(threshold).arg(changeCount).arg(total);
        } else {
            msg = QString("Image difference computed for %1 pixels").arg(total);
        }
        reportProgress(1.0, msg);

        return true;
    }
};

REGISTER_MODULE(ImageDiffModule)

} // namespace aplaceholder
