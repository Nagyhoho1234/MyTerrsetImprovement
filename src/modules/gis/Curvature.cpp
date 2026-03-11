#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class CurvatureModule : public Module {
public:
    QString name() const override { return "CURVATURE"; }
    QString description() const override {
        return "Compute surface curvature from a DEM. "
               "Supports profile curvature (in the direction of maximum slope), "
               "plan curvature (perpendicular to slope direction), and total curvature.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input DEM raster"),
            ParameterDef::output("output", "Output curvature image"),
            ParameterDef::combo("type", "Curvature type",
                {"Profile", "Plan", "Total"}, 2,
                "Profile: along steepest slope; Plan: across slope; Total: combined"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input DEM"); return false; }

        int cols = input->cols(), rows = input->rows();
        int curvType = parameter("type").toInt();
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        const auto& gt = input->geoTransform();
        double cs = std::abs(gt.pixelWidth); // cell size

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing curvature...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;
                if (hasND && data[idx] == noData) {
                    out[idx] = noData;
                    continue;
                }

                auto val = [&](int rr, int cc) -> double {
                    rr = std::max(0, std::min(rows - 1, rr));
                    cc = std::max(0, std::min(cols - 1, cc));
                    size_t i = static_cast<size_t>(rr) * cols + cc;
                    if (hasND && data[i] == noData) return data[idx];
                    return data[i];
                };

                // Zevenbergen-Thorne method
                double z1 = val(r-1, c-1), z2 = val(r-1, c), z3 = val(r-1, c+1);
                double z4 = val(r, c-1),   z5 = data[idx],    z6 = val(r, c+1);
                double z7 = val(r+1, c-1), z8 = val(r+1, c), z9 = val(r+1, c+1);

                double D = ((z4 + z6) / 2.0 - z5) / (cs * cs);
                double E = ((z2 + z8) / 2.0 - z5) / (cs * cs);
                double F = (-z1 + z3 + z7 - z9) / (4.0 * cs * cs);
                double G = (-z4 + z6) / (2.0 * cs);
                double H = (z2 - z8) / (2.0 * cs);

                double result = 0.0;
                if (curvType == 2) {
                    // Total curvature
                    result = -2.0 * (D + E);
                } else {
                    double p = G * G + H * H;
                    if (p < 1e-15) {
                        result = 0.0;
                    } else if (curvType == 0) {
                        // Profile curvature
                        result = -2.0 * (D * G * G + E * H * H + F * G * H) / (p * std::sqrt(1.0 + p));
                    } else {
                        // Plan curvature
                        result = -2.0 * (D * H * H - F * G * H + E * G * G) / std::pow(p, 1.5);
                    }
                }
                out[idx] = result;
            }
            if (r % 100 == 0) reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(CurvatureModule)

} // namespace aplaceholder
