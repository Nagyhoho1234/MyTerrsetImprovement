#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class ZeroPadModule : public Module {
public:
    QString name() const override { return "ZEROPAD"; }
    QString description() const override {
        return "Zero-pads a raster image to the next power-of-2 dimensions in both columns "
               "and rows, as required for FFT-based processing (FOURIER module). Padding is "
               "added symmetrically to mitigate edge effects.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::output("output", "Output zero-padded raster image"),
        };
    }

    bool execute() override {
        // --------------------------------------------------------------------
        // 1. Read input raster
        // --------------------------------------------------------------------
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input image");
            return false;
        }

        int origCols = raster->cols();
        int origRows = raster->rows();
        int numBands = raster->bands();
        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        reportProgress(0.05, "Computing padded dimensions...");

        // --------------------------------------------------------------------
        // 2. Compute next power-of-2 dimensions
        // --------------------------------------------------------------------
        auto nextPow2 = [](int n) -> int {
            if (n <= 0) return 1;
            int p = 1;
            while (p < n) p <<= 1;
            return p;
        };

        int padCols = nextPow2(origCols);
        int padRows = nextPow2(origRows);

        // If already power of 2, no padding needed (but still produce output)
        // Compute offsets for symmetric padding
        int offsetCol = (padCols - origCols) / 2;
        int offsetRow = (padRows - origRows) / 2;

        reportProgress(0.1, QString("Padding from %1x%2 to %3x%4...")
                       .arg(origCols).arg(origRows).arg(padCols).arg(padRows));

        // --------------------------------------------------------------------
        // 3. Create output raster and copy data with zero-padding
        // --------------------------------------------------------------------
        Raster output(padCols, padRows, numBands, raster->dataType());

        // Adjust geo-transform to account for the offset
        GeoTransform gt = raster->geoTransform();
        gt.originX -= offsetCol * gt.pixelWidth;
        gt.originY -= offsetRow * gt.pixelHeight;
        output.setGeoTransform(gt);
        output.setProjection(raster->projection());
        if (hasND)
            output.setNoDataValue(noData);

        int64_t padTotal = static_cast<int64_t>(padCols) * padRows;

        for (int band = 0; band < numBands; ++band) {
            const auto& srcData = raster->data(band);
            auto& dstData = output.data(band);

            // Initialize all to zero
            std::fill(dstData.begin(), dstData.end(), 0.0);

            // Copy original data at the offset position
            for (int row = 0; row < origRows; ++row) {
                int dstRow = row + offsetRow;
                for (int col = 0; col < origCols; ++col) {
                    int dstCol = col + offsetCol;
                    int64_t srcIdx = static_cast<int64_t>(row) * origCols + col;
                    int64_t dstIdx = static_cast<int64_t>(dstRow) * padCols + dstCol;
                    dstData[dstIdx] = srcData[srcIdx];
                }
            }

            reportProgress(0.1 + 0.8 * static_cast<double>(band + 1) / numBands);
        }

        // --------------------------------------------------------------------
        // 4. Write output
        // --------------------------------------------------------------------
        reportProgress(0.9, "Writing output...");
        QString outputPath = parameter("output").toString();
        if (!GdalIO::write(output, outputPath)) {
            setError("Failed to write output: " + outputPath);
            return false;
        }

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(ZeroPadModule)

} // namespace aplaceholder
