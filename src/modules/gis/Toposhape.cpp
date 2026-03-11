#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ToposhapeModule : public Module {
public:
    QString name() const override { return "TOPOSHAPE"; }
    QString description() const override {
        return "Classify topographic shape from a DEM using profile and plan curvature. "
               "Output categories: 1=Peak, 2=Ridge, 3=Saddle, 4=Flat, "
               "5=Ravine, 6=Valley, 7=Pit, 8=Convex slope, 9=Concave slope, 10=Inflection.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input DEM raster"),
            ParameterDef::output("output", "Output topographic shape classification"),
            ParameterDef::real("tolerance", "Curvature tolerance", 0.0001, 0.0, 1.0,
                "Threshold below which curvature is considered zero"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input DEM"); return false; }

        int cols = input->cols(), rows = input->rows();
        double tol = parameter("tolerance").toDouble();
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        double cs = std::abs(input->geoTransform().pixelWidth);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(0.0);
        auto& out = output.data(0);

        reportProgress(0.0, "Classifying topographic shapes...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;
                if (hasND && data[idx] == noData) {
                    out[idx] = 0.0;
                    continue;
                }

                auto val = [&](int rr, int cc) -> double {
                    rr = std::max(0, std::min(rows - 1, rr));
                    cc = std::max(0, std::min(cols - 1, cc));
                    size_t i = static_cast<size_t>(rr) * cols + cc;
                    if (hasND && data[i] == noData) return data[idx];
                    return data[i];
                };

                double z1 = val(r-1,c-1), z2 = val(r-1,c), z3 = val(r-1,c+1);
                double z4 = val(r,c-1),    z5 = data[idx],   z6 = val(r,c+1);
                double z7 = val(r+1,c-1), z8 = val(r+1,c), z9 = val(r+1,c+1);

                // Partial derivatives
                double D = ((z4 + z6) / 2.0 - z5) / (cs * cs);
                double E = ((z2 + z8) / 2.0 - z5) / (cs * cs);
                double F = (-z1 + z3 + z7 - z9) / (4.0 * cs * cs);
                double G = (-z4 + z6) / (2.0 * cs);
                double H = (z2 - z8) / (2.0 * cs);

                double slope2 = G * G + H * H;
                double profileCurv = 0.0, planCurv = 0.0;

                if (slope2 > 1e-15) {
                    profileCurv = -2.0 * (D*G*G + E*H*H + F*G*H) / (slope2 * std::sqrt(1.0 + slope2));
                    planCurv = -2.0 * (D*H*H - F*G*H + E*G*G) / std::pow(slope2, 1.5);
                }

                double slopeMag = std::sqrt(slope2);
                bool flatSlope = (slopeMag < tol * 10);
                int profSign = (std::abs(profileCurv) < tol) ? 0 : (profileCurv > 0 ? 1 : -1);
                int planSign = (std::abs(planCurv) < tol) ? 0 : (planCurv > 0 ? 1 : -1);

                // Classification based on slope, profile and plan curvature
                if (flatSlope) {
                    if (profSign > 0 && planSign > 0) out[idx] = 1.0;      // Peak
                    else if (profSign < 0 && planSign < 0) out[idx] = 7.0;  // Pit
                    else if (profSign > 0 && planSign < 0) out[idx] = 3.0;  // Saddle
                    else if (profSign < 0 && planSign > 0) out[idx] = 3.0;  // Saddle
                    else out[idx] = 4.0;                                     // Flat
                } else {
                    if (profSign > 0 && planSign > 0) out[idx] = 2.0;       // Ridge
                    else if (profSign < 0 && planSign < 0) out[idx] = 6.0;  // Valley
                    else if (profSign > 0 && planSign == 0) out[idx] = 8.0;  // Convex slope
                    else if (profSign < 0 && planSign == 0) out[idx] = 9.0;  // Concave slope
                    else if (profSign > 0 && planSign < 0) out[idx] = 5.0;  // Ravine
                    else if (profSign < 0 && planSign > 0) out[idx] = 5.0;  // Ravine
                    else out[idx] = 10.0;                                     // Inflection
                }
            }
            if (r % 100 == 0) reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ToposhapeModule)

} // namespace aplaceholder
