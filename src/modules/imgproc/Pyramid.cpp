#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class PyramidModule : public Module {
public:
    QString name() const override { return "PYRAMID"; }
    QString description() const override {
        return "Image pyramid (multi-resolution). Creates a series of reduced-resolution "
               "images, each at half the size of the previous level, using pixel averaging.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output_prefix", "Output filename prefix",
                "Output rasters will be named prefix_level1, prefix_level2, etc."),
            ParameterDef::integer("num_levels", "Number of pyramid levels", 4, 1, 16,
                "Number of reduced-resolution levels to generate"),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString prefix = parameter("output_prefix").toString();
        int numLevels = parameter("num_levels").toInt();

        auto current = GdalIO::read(inputPath);
        if (!current) {
            setError("Failed to read input image: " + inputPath);
            return false;
        }

        bool hasND = current->hasNoData();
        double noData = current->noDataValue();
        double outNoData = -9999.0;

        for (int level = 1; level <= numLevels; ++level) {
            reportProgress(static_cast<double>(level - 1) / numLevels,
                           QString("Building pyramid level %1...").arg(level));

            int srcCols = current->cols();
            int srcRows = current->rows();
            int numBands = current->bands();

            int dstCols = srcCols / 2;
            int dstRows = srcRows / 2;

            if (dstCols < 1 || dstRows < 1) {
                reportProgress(1.0, QString("Stopped at level %1: image too small.").arg(level));
                break;
            }

            // Update geo-transform for the reduced resolution
            GeoTransform gt = current->geoTransform();
            gt.pixelWidth *= 2.0;
            gt.pixelHeight *= 2.0;

            auto reduced = std::make_unique<Raster>(dstCols, dstRows, numBands, DataType::Float32);
            reduced->setGeoTransform(gt);
            reduced->setProjection(current->projection());
            reduced->setNoDataValue(outNoData);

            for (int b = 0; b < numBands; ++b) {
                const auto& srcData = current->data(b);
                auto& dstData = reduced->data(b);

                for (int dr = 0; dr < dstRows; ++dr) {
                    for (int dc = 0; dc < dstCols; ++dc) {
                        int sr = dr * 2;
                        int sc = dc * 2;

                        // Average 2x2 block, skipping nodata
                        double sum = 0.0;
                        int count = 0;

                        for (int yr = 0; yr < 2; ++yr) {
                            for (int xc = 0; xc < 2; ++xc) {
                                int64_t si = static_cast<int64_t>(sr + yr) * srcCols + (sc + xc);
                                double v = srcData[si];
                                if (hasND && v == noData) continue;
                                sum += v;
                                ++count;
                            }
                        }

                        int64_t di = static_cast<int64_t>(dr) * dstCols + dc;
                        dstData[di] = (count > 0) ? sum / count : outNoData;
                    }
                }
            }

            QString outPath = prefix + "_level" + QString::number(level);
            if (!GdalIO::write(*reduced, outPath)) {
                setError("Failed to write pyramid level: " + outPath);
                return false;
            }

            current = std::move(reduced);
        }

        reportProgress(1.0, "Image pyramid complete.");
        return true;
    }
};

REGISTER_MODULE(PyramidModule)

} // namespace aplaceholder
