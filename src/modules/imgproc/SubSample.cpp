#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class SubSampleModule : public Module {
public:
    QString name() const override { return "SUBSAMPLE"; }
    QString description() const override {
        return "Spatial subsampling. Extracts every Nth pixel from the input image, "
               "reducing spatial resolution by the given factor.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output subsampled image"),
            ParameterDef::integer("factor", "Subsampling factor", 2, 2, 100,
                "Extract every Nth pixel (e.g., 2 means half resolution)"),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString outputPath = parameter("output").toString();
        int factor = parameter("factor").toInt();

        auto input = GdalIO::read(inputPath);
        if (!input) {
            setError("Failed to read input image: " + inputPath);
            return false;
        }

        int srcCols = input->cols();
        int srcRows = input->rows();
        int numBands = input->bands();

        int dstCols = (srcCols + factor - 1) / factor;
        int dstRows = (srcRows + factor - 1) / factor;

        bool hasND = input->hasNoData();
        double inputND = input->noDataValue();

        // Update geo-transform
        GeoTransform gt = input->geoTransform();
        gt.pixelWidth *= factor;
        gt.pixelHeight *= factor;

        Raster output(dstCols, dstRows, numBands, DataType::Float32);
        output.setGeoTransform(gt);
        output.setProjection(input->projection());
        if (hasND)
            output.setNoDataValue(inputND);

        reportProgress(0.0, "Subsampling...");

        for (int b = 0; b < numBands; ++b) {
            const auto& srcData = input->data(b);
            auto& dstData = output.data(b);

            for (int dr = 0; dr < dstRows; ++dr) {
                int sr = dr * factor;
                for (int dc = 0; dc < dstCols; ++dc) {
                    int sc = dc * factor;
                    int64_t srcIdx = static_cast<int64_t>(sr) * srcCols + sc;
                    int64_t dstIdx = static_cast<int64_t>(dr) * dstCols + dc;
                    dstData[dstIdx] = srcData[srcIdx];
                }
            }

            reportProgress(static_cast<double>(b + 1) / numBands);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
    }
};

REGISTER_MODULE(SubSampleModule)

} // namespace aplaceholder
