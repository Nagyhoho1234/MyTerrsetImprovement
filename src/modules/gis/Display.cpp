#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class DisplayModule : public Module {
public:
    QString name() const override { return "DISPLAY"; }
    QString description() const override {
        return "Display launcher. Opens a raster image for visualization by generating "
               "an auto-stretched RGB or grayscale preview. Supports single-band "
               "grayscale, pseudo-color, and RGB composite display modes.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster to display"),
            ParameterDef::output("output", "Output display image (8-bit stretched)"),
            ParameterDef::combo("mode", "Display mode",
                {"Grayscale", "Pseudo-color", "RGB Composite"}, 0,
                "Visualization mode"),
            ParameterDef::integer("red_band", "Red band (for RGB)", 1, 1, 256,
                "Band number for red channel in RGB composite"),
            ParameterDef::integer("green_band", "Green band (for RGB)", 2, 1, 256,
                "Band number for green channel in RGB composite"),
            ParameterDef::integer("blue_band", "Blue band (for RGB)", 3, 1, 256,
                "Band number for blue channel in RGB composite"),
            ParameterDef::real("stretch_min_pct", "Stretch min percentile", 2.0, 0.0, 50.0,
                "Lower percentile for contrast stretch"),
            ParameterDef::real("stretch_max_pct", "Stretch max percentile", 98.0, 50.0, 100.0,
                "Upper percentile for contrast stretch"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int mode = parameter("mode").toInt();

        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        reportProgress(0.05, "Preparing display...");

        // Helper: linear stretch a band to 0-255 using percentile clipping
        auto stretchBand = [&](int bandIdx) -> std::vector<double> {
            const auto& src = raster->data(bandIdx);
            std::vector<double> valid;
            valid.reserve(total);
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && src[i] == noData) continue;
                valid.push_back(src[i]);
            }
            std::sort(valid.begin(), valid.end());

            double lo = 0.0, hi = 255.0;
            if (!valid.empty()) {
                double minPct = parameter("stretch_min_pct").toDouble() / 100.0;
                double maxPct = parameter("stretch_max_pct").toDouble() / 100.0;
                int64_t nv = static_cast<int64_t>(valid.size());
                lo = valid[std::min(static_cast<int64_t>(minPct * nv), nv - 1)];
                hi = valid[std::min(static_cast<int64_t>(maxPct * nv), nv - 1)];
            }

            std::vector<double> result(total);
            double range = hi - lo;
            if (range < 1e-10) range = 1.0;
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && src[i] == noData) {
                    result[i] = 0.0;
                } else {
                    result[i] = std::clamp((src[i] - lo) / range * 255.0, 0.0, 255.0);
                }
            }
            return result;
        };

        int outBands = 1;
        if (mode == 2) outBands = 3;  // RGB

        Raster output(cols, rows, outBands, DataType::Byte);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());

        if (mode == 0) {
            // Grayscale: stretch band 0
            auto stretched = stretchBand(0);
            auto& dst = output.data(0);
            for (int64_t i = 0; i < total; ++i)
                dst[i] = stretched[i];
        } else if (mode == 1) {
            // Pseudo-color: stretch band 0, output as single band
            // (actual colormap rendering is a GUI concern; we produce the stretched values)
            auto stretched = stretchBand(0);
            auto& dst = output.data(0);
            for (int64_t i = 0; i < total; ++i)
                dst[i] = stretched[i];
        } else if (mode == 2) {
            // RGB composite
            int rBand = parameter("red_band").toInt() - 1;
            int gBand = parameter("green_band").toInt() - 1;
            int bBand = parameter("blue_band").toInt() - 1;

            rBand = std::clamp(rBand, 0, raster->bands() - 1);
            gBand = std::clamp(gBand, 0, raster->bands() - 1);
            bBand = std::clamp(bBand, 0, raster->bands() - 1);

            reportProgress(0.20, "Stretching red band...");
            auto rStr = stretchBand(rBand);
            reportProgress(0.40, "Stretching green band...");
            auto gStr = stretchBand(gBand);
            reportProgress(0.60, "Stretching blue band...");
            auto bStr = stretchBand(bBand);

            auto& dstR = output.data(0);
            auto& dstG = output.data(1);
            auto& dstB = output.data(2);
            for (int64_t i = 0; i < total; ++i) {
                dstR[i] = rStr[i];
                dstG[i] = gStr[i];
                dstB[i] = bStr[i];
            }
        }

        reportProgress(1.0, "Writing display output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(DisplayModule)

} // namespace aplaceholder
