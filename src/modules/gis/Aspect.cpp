#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class AspectModule : public Module {
public:
    QString name() const override { return "ASPECT"; }
    QString description() const override {
        return "Compute slope aspect (direction of steepest descent) from a DEM. "
               "Output is in degrees clockwise from north (0-360), with flat areas "
               "assigned a value of -1.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input DEM raster"),
            ParameterDef::output("output", "Output aspect image (degrees)"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input DEM"); return false; }

        int cols = input->cols(), rows = input->rows();
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        const auto& gt = input->geoTransform();
        double cellX = std::abs(gt.pixelWidth);
        double cellY = std::abs(gt.pixelHeight);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(-1.0);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing aspect...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;

                if (hasND && data[idx] == noData) {
                    out[idx] = -1.0;
                    continue;
                }

                // 3x3 neighborhood with boundary replication
                auto val = [&](int rr, int cc) -> double {
                    rr = std::max(0, std::min(rows - 1, rr));
                    cc = std::max(0, std::min(cols - 1, cc));
                    size_t i = static_cast<size_t>(rr) * cols + cc;
                    if (hasND && data[i] == noData) return data[idx];
                    return data[i];
                };

                // Horn's method for gradient
                double dzdx = ((val(r-1,c+1) + 2*val(r,c+1) + val(r+1,c+1)) -
                               (val(r-1,c-1) + 2*val(r,c-1) + val(r+1,c-1))) / (8.0 * cellX);
                double dzdy = ((val(r+1,c-1) + 2*val(r+1,c) + val(r+1,c+1)) -
                               (val(r-1,c-1) + 2*val(r-1,c) + val(r-1,c+1))) / (8.0 * cellY);

                if (std::abs(dzdx) < 1e-10 && std::abs(dzdy) < 1e-10) {
                    out[idx] = -1.0; // Flat
                } else {
                    // atan2 gives angle from east, convert to degrees CW from north
                    double aspect = std::atan2(dzdy, -dzdx) * 180.0 / M_PI;
                    // Convert: mathematical angle to compass bearing
                    aspect = 90.0 - aspect;
                    if (aspect < 0.0) aspect += 360.0;
                    if (aspect >= 360.0) aspect -= 360.0;
                    out[idx] = aspect;
                }
            }
            if (r % 100 == 0) reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(AspectModule)

} // namespace aplaceholder
