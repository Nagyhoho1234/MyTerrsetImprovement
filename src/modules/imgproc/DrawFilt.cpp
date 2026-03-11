#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class DrawFiltModule : public Module {
public:
    QString name() const override { return "DRAWFILT"; }
    QString description() const override {
        return "Design frequency domain filter mask. Creates a filter raster with "
               "specified geometry (circle, ring, rectangle, wedge) for use with "
               "frequency domain filtering operations.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("reference_raster", "Reference raster (for dimensions)",
                "A raster whose dimensions define the filter size (e.g., FFT magnitude output)"),
            ParameterDef::output("output", "Output filter mask image"),
            ParameterDef::combo("shape", "Filter shape",
                {"circle", "ring", "rectangle", "wedge"}, 0,
                "Geometric shape of the filter region"),
            ParameterDef::real("center_x", "Center X (fraction 0-1)", 0.5, 0.0, 1.0,
                "Horizontal center of filter as fraction of image width"),
            ParameterDef::real("center_y", "Center Y (fraction 0-1)", 0.5, 0.0, 1.0,
                "Vertical center of filter as fraction of image height"),
            ParameterDef::real("radius1", "Radius 1 (fraction)", 0.25, 0.0, 1.0,
                "Primary radius as fraction of image diagonal. For ring: inner radius. "
                "For rectangle: half-width. For wedge: radius."),
            ParameterDef::real("radius2", "Radius 2 (fraction)", 0.5, 0.0, 1.0,
                "Secondary radius as fraction of image diagonal. For ring: outer radius. "
                "For rectangle: half-height. For wedge: angular half-width in degrees (0-180)."),
            ParameterDef::combo("pass_type", "Pass type",
                {"pass", "reject"}, 0,
                "Pass: keep frequencies inside shape. Reject: remove frequencies inside shape."),
        };
    }

    bool execute() override {
        // Read reference raster for dimensions only
        auto refRaster = GdalIO::readMetadata(parameter("reference_raster").toString());
        if (!refRaster) {
            setError("Failed to read reference raster");
            return false;
        }

        int cols = refRaster->cols();
        int rows = refRaster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        int shapeIdx = parameter("shape").toInt();
        double cx = parameter("center_x").toDouble() * cols;
        double cy = parameter("center_y").toDouble() * rows;
        double r1 = parameter("radius1").toDouble();
        double r2 = parameter("radius2").toDouble();
        int passIdx = parameter("pass_type").toInt();
        bool passMode = (passIdx == 0); // pass = inside is 1; reject = inside is 0

        double diag = std::sqrt(static_cast<double>(cols * cols + rows * rows));
        double absR1 = r1 * diag;
        double absR2 = r2 * diag;

        reportProgress(0.1, "Generating filter mask...");

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(refRaster->geoTransform());
        output.setProjection(refRaster->projection());
        auto& outData = output.data(0);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                double dx = c - cx;
                double dy = r - cy;
                double dist = std::sqrt(dx * dx + dy * dy);

                bool inside = false;

                switch (shapeIdx) {
                case 0: // circle
                    inside = (dist <= absR1);
                    break;

                case 1: // ring
                    inside = (dist >= absR1 && dist <= absR2);
                    break;

                case 2: { // rectangle
                    double halfW = r1 * cols;
                    double halfH = r2 * rows;
                    inside = (std::abs(dx) <= halfW && std::abs(dy) <= halfH);
                    break;
                }

                case 3: { // wedge
                    // r1 = radius, r2 = angular half-width in degrees (stored as fraction,
                    // so we interpret r2 * 180 as degrees)
                    double angle = std::atan2(dy, dx) * 180.0 / M_PI; // -180 to 180
                    double halfAngle = r2 * 180.0; // convert fraction to degrees
                    // Wedge centered on the positive x-axis from center
                    inside = (dist <= absR1) &&
                             (std::abs(angle) <= halfAngle);
                    break;
                }
                }

                if (passMode) {
                    outData[idx] = inside ? 1.0 : 0.0;
                } else {
                    outData[idx] = inside ? 0.0 : 1.0;
                }
            }

            if (r % 100 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Filter mask complete.");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(DrawFiltModule)

} // namespace aplaceholder
