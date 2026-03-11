#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class SurfaceModule : public Module {
public:
    QString name() const override { return "SURFACE"; }
    QString description() const override {
        return "Calculates slope, aspect, or analytical hillshading from a DEM.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input DEM"),
            ParameterDef::output("output", "Output image"),
            ParameterDef::combo("mode", "Output type",
                {"Slope (degrees)", "Slope (percent)", "Aspect (degrees)",
                 "Hillshade", "Profile Curvature", "Plan Curvature"}, 0),
            ParameterDef::real("azimuth", "Sun azimuth (hillshade)", 315, 0, 360,
                "Azimuth of light source in degrees"),
            ParameterDef::real("altitude", "Sun altitude (hillshade)", 45, 0, 90,
                "Altitude of light source in degrees"),
            ParameterDef::real("z_factor", "Z factor", 1.0, 0.0001, 100000,
                "Vertical exaggeration factor"),
        };
    }

    bool execute() override {
        auto dem = GdalIO::read(parameter("input").toString());
        if (!dem) { setError("Failed to read DEM"); return false; }

        int cols = dem->cols(), rows = dem->rows();
        int mode = parameter("mode").toInt();
        double zFactor = parameter("z_factor").toDouble();
        double azimuth = parameter("azimuth").toDouble() * M_PI / 180.0;
        double altitude = parameter("altitude").toDouble() * M_PI / 180.0;

        double dx = std::abs(dem->geoTransform().pixelWidth);
        double dy = std::abs(dem->geoTransform().pixelHeight);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(dem->geoTransform());
        output.setProjection(dem->projection());
        output.setNoDataValue(-9999);

        const auto& elev = dem->data(0);
        auto& out = output.data(0);

        for (int r = 1; r < rows - 1; ++r) {
            for (int c = 1; c < cols - 1; ++c) {
                // 3x3 neighborhood (Horn's method)
                double z1 = elev[(r-1)*cols+(c-1)] * zFactor;
                double z2 = elev[(r-1)*cols+c]     * zFactor;
                double z3 = elev[(r-1)*cols+(c+1)] * zFactor;
                double z4 = elev[r*cols+(c-1)]     * zFactor;
                double z6 = elev[r*cols+(c+1)]     * zFactor;
                double z7 = elev[(r+1)*cols+(c-1)] * zFactor;
                double z8 = elev[(r+1)*cols+c]     * zFactor;
                double z9 = elev[(r+1)*cols+(c+1)] * zFactor;

                double dzdx = ((z3 + 2*z6 + z9) - (z1 + 2*z4 + z7)) / (8.0 * dx);
                double dzdy = ((z7 + 2*z8 + z9) - (z1 + 2*z2 + z3)) / (8.0 * dy);

                size_t idx = static_cast<size_t>(r) * cols + c;

                switch (mode) {
                case 0: { // Slope degrees
                    double slope = std::atan(std::sqrt(dzdx*dzdx + dzdy*dzdy));
                    out[idx] = slope * 180.0 / M_PI;
                    break;
                }
                case 1: { // Slope percent
                    out[idx] = std::sqrt(dzdx*dzdx + dzdy*dzdy) * 100.0;
                    break;
                }
                case 2: { // Aspect
                    double aspect = std::atan2(dzdy, -dzdx) * 180.0 / M_PI;
                    if (aspect < 0) aspect += 360.0;
                    out[idx] = aspect;
                    break;
                }
                case 3: { // Hillshade
                    double slope = std::atan(std::sqrt(dzdx*dzdx + dzdy*dzdy));
                    double aspect = std::atan2(dzdy, -dzdx);
                    double hs = 255.0 * (std::cos(altitude) * std::cos(slope) +
                                std::sin(altitude) * std::sin(slope) *
                                std::cos(azimuth - aspect));
                    out[idx] = std::clamp(hs, 0.0, 255.0);
                    break;
                }
                case 4: { // Profile curvature
                    double p = dzdx * dzdx + dzdy * dzdy;
                    if (p < 1e-10) { out[idx] = 0; break; }
                    double dzdxx = ((z3 - 2*elev[r*cols+c]*zFactor + z4)) / (dx*dx);
                    double dzdyy = ((z2 - 2*elev[r*cols+c]*zFactor + z8)) / (dy*dy);
                    double dzdxy = ((z3 - z1 - z9 + z7) * zFactor) / (4*dx*dy);
                    out[idx] = -(dzdx*dzdx*dzdxx + 2*dzdx*dzdy*dzdxy + dzdy*dzdy*dzdyy) /
                               (p * std::pow(1 + p, 1.5));
                    break;
                }
                case 5: { // Plan curvature
                    double p = dzdx * dzdx + dzdy * dzdy;
                    if (p < 1e-10) { out[idx] = 0; break; }
                    double dzdxx = ((z6 - 2*elev[r*cols+c]*zFactor + z4)) / (dx*dx);
                    double dzdyy = ((z2 - 2*elev[r*cols+c]*zFactor + z8)) / (dy*dy);
                    double dzdxy = ((z3 - z1 - z9 + z7) * zFactor) / (4*dx*dy);
                    out[idx] = (dzdy*dzdy*dzdxx - 2*dzdx*dzdy*dzdxy + dzdx*dzdx*dzdyy) /
                               (p * std::pow(1 + p, 0.5));
                    break;
                }
                }
            }
            if (r % 100 == 0)
                reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SurfaceModule)

} // namespace aplaceholder
