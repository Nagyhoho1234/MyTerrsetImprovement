#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class LandsatUpsampleModule : public Module {
public:
    QString name() const override { return "LANDSAT_UPSAMPLE"; }
    QString description() const override {
        return "Upsample Landsat thermal/SWIR bands to match higher resolution bands. "
               "Uses bilinear interpolation to resample the low-resolution input band "
               "to the dimensions of a high-resolution reference band.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_low", "Input low-resolution band"),
            ParameterDef::file("reference_high", "High-resolution reference band"),
            ParameterDef::output("output", "Output upsampled image"),
        };
    }

    bool execute() override {
        auto lowRes = GdalIO::read(parameter("input_low").toString());
        if (!lowRes) {
            setError("Failed to read input low-resolution raster");
            return false;
        }

        auto highRes = GdalIO::read(parameter("reference_high").toString());
        if (!highRes) {
            setError("Failed to read high-resolution reference raster");
            return false;
        }

        int srcCols = lowRes->cols(), srcRows = lowRes->rows();
        int dstCols = highRes->cols(), dstRows = highRes->rows();
        int64_t dstTotal = static_cast<int64_t>(dstCols) * dstRows;

        Raster output(dstCols, dstRows, 1, DataType::Float64);
        output.setGeoTransform(highRes->geoTransform());
        output.setProjection(highRes->projection());
        output.setNoDataValue(lowRes->noDataValue());

        const auto& src = lowRes->data(0);
        auto& dst = output.data(0);
        double noData = lowRes->noDataValue();
        bool hasND = lowRes->hasNoData();

        // Compute scale factors from destination to source coordinates
        double scaleX = static_cast<double>(srcCols) / dstCols;
        double scaleY = static_cast<double>(srcRows) / dstRows;

        reportProgress(0.0, "Upsampling via bilinear interpolation...");

        for (int dstR = 0; dstR < dstRows; ++dstR) {
            for (int dstC = 0; dstC < dstCols; ++dstC) {
                int64_t dstIdx = static_cast<int64_t>(dstR) * dstCols + dstC;

                // Map destination pixel to source coordinates
                double srcX = (dstC + 0.5) * scaleX - 0.5;
                double srcY = (dstR + 0.5) * scaleY - 0.5;

                int x0 = static_cast<int>(std::floor(srcX));
                int y0 = static_cast<int>(std::floor(srcY));
                int x1 = x0 + 1;
                int y1 = y0 + 1;

                // Clamp to source raster bounds
                x0 = std::clamp(x0, 0, srcCols - 1);
                x1 = std::clamp(x1, 0, srcCols - 1);
                y0 = std::clamp(y0, 0, srcRows - 1);
                y1 = std::clamp(y1, 0, srcRows - 1);

                double fx = srcX - std::floor(srcX);
                double fy = srcY - std::floor(srcY);

                double v00 = src[static_cast<int64_t>(y0) * srcCols + x0];
                double v10 = src[static_cast<int64_t>(y0) * srcCols + x1];
                double v01 = src[static_cast<int64_t>(y1) * srcCols + x0];
                double v11 = src[static_cast<int64_t>(y1) * srcCols + x1];

                // Check for nodata in any of the four neighbors
                if (hasND && (v00 == noData || v10 == noData ||
                              v01 == noData || v11 == noData)) {
                    dst[dstIdx] = noData;
                    continue;
                }

                // Bilinear interpolation
                double top = v00 + fx * (v10 - v00);
                double bot = v01 + fx * (v11 - v01);
                dst[dstIdx] = top + fy * (bot - top);
            }

            if (dstR % 100 == 0)
                reportProgress(static_cast<double>(dstR) / dstRows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(LandsatUpsampleModule)

} // namespace aplaceholder
