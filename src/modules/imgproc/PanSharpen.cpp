#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class PanSharpenModule : public Module {
public:
    QString name() const override { return "PANSHARPEN"; }
    QString description() const override {
        return "IHS pan-sharpening. Converts RGB bands to IHS color space, replaces "
               "the Intensity component with a high-resolution panchromatic band, "
               "and converts back to RGB to produce a sharpened color image.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("red", "Red band image"),
            ParameterDef::file("green", "Green band image"),
            ParameterDef::file("blue", "Blue band image"),
            ParameterDef::file("pan_band", "Panchromatic (high-res) band"),
            ParameterDef::output("output_prefix", "Output filename prefix",
                "Output rasters will be named prefix_red, prefix_green, prefix_blue"),
        };
    }

    bool execute() override {
        QString redPath = parameter("red").toString();
        QString greenPath = parameter("green").toString();
        QString bluePath = parameter("blue").toString();
        QString panPath = parameter("pan_band").toString();
        QString prefix = parameter("output_prefix").toString();

        // Read inputs
        reportProgress(0.0, "Reading input bands...");
        auto redRst = GdalIO::read(redPath);
        auto greenRst = GdalIO::read(greenPath);
        auto blueRst = GdalIO::read(bluePath);
        auto panRst = GdalIO::read(panPath);

        if (!redRst || !greenRst || !blueRst || !panRst) {
            setError("Failed to read one or more input images.");
            return false;
        }

        // Use pan dimensions as output dimensions
        int cols = panRst->cols();
        int rows = panRst->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // The multispectral bands must match each other
        if (redRst->cols() != greenRst->cols() || redRst->cols() != blueRst->cols() ||
            redRst->rows() != greenRst->rows() || redRst->rows() != blueRst->rows()) {
            setError("Red, green, and blue band images must have the same dimensions.");
            return false;
        }

        int msCols = redRst->cols();
        int msRows = redRst->rows();

        // Compute scale factors for resampling multispectral to pan resolution
        double scaleX = static_cast<double>(msCols) / cols;
        double scaleY = static_cast<double>(msRows) / rows;

        bool hasNoData = panRst->hasNoData();
        double noData = panRst->noDataValue();
        double outNoData = -9999.0;

        const auto& panData = panRst->data(0);
        const auto& rData = redRst->data(0);
        const auto& gData = greenRst->data(0);
        const auto& bData = blueRst->data(0);

        // Create output rasters
        Raster outRed(cols, rows, 1, DataType::Float32);
        Raster outGreen(cols, rows, 1, DataType::Float32);
        Raster outBlue(cols, rows, 1, DataType::Float32);

        for (auto* r : {&outRed, &outGreen, &outBlue}) {
            r->setGeoTransform(panRst->geoTransform());
            r->setProjection(panRst->projection());
            r->setNoDataValue(outNoData);
        }

        auto& oR = outRed.data(0);
        auto& oG = outGreen.data(0);
        auto& oB = outBlue.data(0);

        reportProgress(0.1, "Pan-sharpening via IHS transform...");

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                int64_t idx = static_cast<int64_t>(row) * cols + col;
                double pan = panData[idx];

                if (hasNoData && pan == noData) {
                    oR[idx] = outNoData;
                    oG[idx] = outNoData;
                    oB[idx] = outNoData;
                    continue;
                }

                // Resample multispectral to pan grid using nearest neighbor
                int msCol = std::min(static_cast<int>(col * scaleX), msCols - 1);
                int msRow = std::min(static_cast<int>(row * scaleY), msRows - 1);
                int64_t msIdx = static_cast<int64_t>(msRow) * msCols + msCol;

                double r = rData[msIdx];
                double g = gData[msIdx];
                double b = bData[msIdx];

                // RGB to IHS
                double I = (r + g + b) / 3.0;
                double v1 = (-r - g + 2.0 * b) / (2.0 * std::sqrt(2.0));  // related to H,S
                double v2 = (r - g) / std::sqrt(2.0);

                // Replace intensity with pan band
                double newI = pan;

                // IHS to RGB: add the chromaticity difference
                double diff = newI - I;
                double newR = r + diff;
                double newG = g + diff;
                double newB = b + diff;

                oR[idx] = newR;
                oG[idx] = newG;
                oB[idx] = newB;
            }

            if (row % 500 == 0)
                reportProgress(0.1 + 0.7 * static_cast<double>(row) / rows);
        }

        // Write outputs
        reportProgress(0.8, "Writing output bands...");
        if (!GdalIO::write(outRed, prefix + "_red") ||
            !GdalIO::write(outGreen, prefix + "_green") ||
            !GdalIO::write(outBlue, prefix + "_blue")) {
            setError("Failed to write output rasters.");
            return false;
        }

        reportProgress(1.0, "Pan-sharpening complete.");
        return true;
    }
};

REGISTER_MODULE(PanSharpenModule)

} // namespace aplaceholder
