#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class FacetModule : public Module {
public:
    QString name() const override { return "FACET"; }
    QString description() const override {
        return "Terrain facet classification. Classifies terrain into flat, N, NE, E, SE, "
               "S, SW, W, NW facing based on slope and aspect from a DEM.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dem", "Input DEM"),
            ParameterDef::output("output", "Output facet classification raster"),
            ParameterDef::real("slope_threshold", "Slope threshold for flat (degrees)",
                5.0, 0.0, 90.0, "Cells with slope below this are classified as flat"),
        };
    }

    bool execute() override {
        auto dem = GdalIO::read(parameter("dem").toString());
        if (!dem) { setError("Failed to read DEM"); return false; }

        int cols = dem->cols(), rows = dem->rows();
        const auto& elev = dem->data(0);
        double noData = dem->noDataValue();
        bool hasND = dem->hasNoData();

        double dx = std::abs(dem->geoTransform().pixelWidth);
        double dy = std::abs(dem->geoTransform().pixelHeight);
        double slopeThreshold = parameter("slope_threshold").toDouble();

        // Classification codes:
        // 0 = NoData, 1 = Flat, 2 = N, 3 = NE, 4 = E, 5 = SE, 6 = S, 7 = SW, 8 = W, 9 = NW
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(dem->geoTransform());
        output.setProjection(dem->projection());
        output.setNoDataValue(0);
        auto& out = output.data(0);

        reportProgress(0.0, "Classifying terrain facets...");
        for (int r = 1; r < rows - 1; ++r) {
            for (int c = 1; c < cols - 1; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;

                if (hasND && elev[idx] == noData) {
                    out[idx] = 0;
                    continue;
                }

                // Horn's method
                double z1 = elev[(r-1)*cols+(c-1)];
                double z2 = elev[(r-1)*cols+c];
                double z3 = elev[(r-1)*cols+(c+1)];
                double z4 = elev[r*cols+(c-1)];
                double z6 = elev[r*cols+(c+1)];
                double z7 = elev[(r+1)*cols+(c-1)];
                double z8 = elev[(r+1)*cols+c];
                double z9 = elev[(r+1)*cols+(c+1)];

                double dzdx = ((z3 + 2*z6 + z9) - (z1 + 2*z4 + z7)) / (8.0 * dx);
                double dzdy = ((z7 + 2*z8 + z9) - (z1 + 2*z2 + z3)) / (8.0 * dy);

                double slopeDeg = std::atan(std::sqrt(dzdx*dzdx + dzdy*dzdy)) * 180.0 / M_PI;

                if (slopeDeg < slopeThreshold) {
                    out[idx] = 1; // Flat
                } else {
                    double aspect = std::atan2(dzdy, -dzdx) * 180.0 / M_PI;
                    if (aspect < 0) aspect += 360.0;

                    // Classify into 8 cardinal/ordinal directions
                    // N=0/360, NE=45, E=90, SE=135, S=180, SW=225, W=270, NW=315
                    if (aspect >= 337.5 || aspect < 22.5)
                        out[idx] = 2; // N
                    else if (aspect < 67.5)
                        out[idx] = 3; // NE
                    else if (aspect < 112.5)
                        out[idx] = 4; // E
                    else if (aspect < 157.5)
                        out[idx] = 5; // SE
                    else if (aspect < 202.5)
                        out[idx] = 6; // S
                    else if (aspect < 247.5)
                        out[idx] = 7; // SW
                    else if (aspect < 292.5)
                        out[idx] = 8; // W
                    else
                        out[idx] = 9; // NW
                }
            }
            if (r % 100 == 0)
                reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(FacetModule)

} // namespace aplaceholder
