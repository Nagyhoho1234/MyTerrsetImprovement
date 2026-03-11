#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class HillshadeModule : public Module {
public:
    QString name() const override { return "HILLSHADE"; }
    QString description() const override {
        return "Standalone Analytical Hillshading. Creates a shaded relief image from "
               "a DEM using Horn's method for slope and aspect calculation. Output is "
               "scaled 0-255.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dem", "Input DEM"),
            ParameterDef::output("output", "Output hillshade image"),
            ParameterDef::real("sun_azimuth", "Sun azimuth (degrees)", 315.0, 0.0, 360.0,
                "Direction of illumination in degrees clockwise from north"),
            ParameterDef::real("sun_elevation", "Sun elevation (degrees)", 45.0, 0.0, 90.0,
                "Angle of illumination above the horizon in degrees"),
            ParameterDef::real("z_factor", "Z factor (vertical exaggeration)", 1.0, 0.0001, 100000.0,
                "Vertical exaggeration factor for the elevation values"),
        };
    }

    bool execute() override {
        auto dem = GdalIO::read(parameter("dem").toString());
        if (!dem) {
            setError("Failed to read DEM");
            return false;
        }

        int cols = dem->cols(), rows = dem->rows();
        double azimuthDeg = parameter("sun_azimuth").toDouble();
        double elevationDeg = parameter("sun_elevation").toDouble();
        double zFactor = parameter("z_factor").toDouble();

        // Convert to radians
        double azimuth = azimuthDeg * M_PI / 180.0;
        double zenith = (90.0 - elevationDeg) * M_PI / 180.0;

        double dx = std::abs(dem->geoTransform().pixelWidth);
        double dy = std::abs(dem->geoTransform().pixelHeight);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(dem->geoTransform());
        output.setProjection(dem->projection());
        output.setNoDataValue(-9999.0);

        const auto& elev = dem->data(0);
        auto& out = output.data(0);
        bool hasND = dem->hasNoData();
        double noData = dem->noDataValue();

        reportProgress(0.0, "Computing hillshade...");

        // Set border pixels to nodata
        for (int c = 0; c < cols; ++c) {
            out[c] = -9999.0;
            out[static_cast<size_t>(rows - 1) * cols + c] = -9999.0;
        }
        for (int r = 0; r < rows; ++r) {
            out[static_cast<size_t>(r) * cols] = -9999.0;
            out[static_cast<size_t>(r) * cols + (cols - 1)] = -9999.0;
        }

        for (int r = 1; r < rows - 1; ++r) {
            for (int c = 1; c < cols - 1; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;

                // Check for nodata in 3x3 neighborhood
                if (hasND) {
                    bool hasNoData = false;
                    for (int dr = -1; dr <= 1 && !hasNoData; ++dr)
                        for (int dc = -1; dc <= 1 && !hasNoData; ++dc)
                            if (elev[static_cast<size_t>(r + dr) * cols + (c + dc)] == noData)
                                hasNoData = true;
                    if (hasNoData) {
                        out[idx] = -9999.0;
                        continue;
                    }
                }

                // 3x3 neighborhood with Horn's method
                double z1 = elev[(r-1)*cols+(c-1)] * zFactor;
                double z2 = elev[(r-1)*cols+c]     * zFactor;
                double z3 = elev[(r-1)*cols+(c+1)] * zFactor;
                double z4 = elev[r*cols+(c-1)]     * zFactor;
                double z6 = elev[r*cols+(c+1)]     * zFactor;
                double z7 = elev[(r+1)*cols+(c-1)] * zFactor;
                double z8 = elev[(r+1)*cols+c]     * zFactor;
                double z9 = elev[(r+1)*cols+(c+1)] * zFactor;

                double dzdx = ((z3 + 2.0*z6 + z9) - (z1 + 2.0*z4 + z7)) / (8.0 * dx);
                double dzdy = ((z7 + 2.0*z8 + z9) - (z1 + 2.0*z2 + z3)) / (8.0 * dy);

                // Slope and aspect
                double slope = std::atan(std::sqrt(dzdx*dzdx + dzdy*dzdy));
                double aspect = std::atan2(dzdy, -dzdx);

                // Hillshade formula
                double shade = std::cos(zenith) * std::cos(slope) +
                               std::sin(zenith) * std::sin(slope) *
                               std::cos(azimuth - aspect);

                // Scale to 0-255
                out[idx] = std::clamp(shade * 255.0, 0.0, 255.0);
            }

            if (r % 100 == 0)
                reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(HillshadeModule)

} // namespace aplaceholder
