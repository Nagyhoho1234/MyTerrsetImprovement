#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class SedimentationModule : public Module {
public:
    QString name() const override { return "SEDIMENTATION"; }
    QString description() const override {
        return "Sediment transport modeling based on slope and contributing area. "
               "Computes sediment transport capacity and net erosion/deposition "
               "using a stream power or transport-limited approach.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dem", "Input DEM raster"),
            ParameterDef::file("flow_accum", "Flow accumulation raster"),
            ParameterDef::output("output", "Output sediment transport capacity image"),
            ParameterDef::real("k_factor", "Erodibility factor (K)", 0.03, 0.0, 1.0,
                "Soil erodibility factor"),
            ParameterDef::real("m_exp", "Area exponent (m)", 1.6, 0.0, 5.0,
                "Exponent for contributing area term"),
            ParameterDef::real("n_exp", "Slope exponent (n)", 1.3, 0.0, 5.0,
                "Exponent for slope term"),
        };
    }

    bool execute() override {
        auto dem = GdalIO::read(parameter("dem").toString());
        if (!dem) { setError("Failed to read DEM"); return false; }

        auto flowAcc = GdalIO::read(parameter("flow_accum").toString());
        if (!flowAcc) { setError("Failed to read flow accumulation raster"); return false; }

        int cols = dem->cols(), rows = dem->rows();
        if (flowAcc->cols() != cols || flowAcc->rows() != rows) {
            setError("DEM and flow accumulation rasters must have the same dimensions");
            return false;
        }

        double K = parameter("k_factor").toDouble();
        double mExp = parameter("m_exp").toDouble();
        double nExp = parameter("n_exp").toDouble();

        const auto& demData = dem->data(0);
        const auto& accData = flowAcc->data(0);
        double noData = dem->noDataValue();
        bool hasND = dem->hasNoData();
        double cs = std::abs(dem->geoTransform().pixelWidth);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(dem->geoTransform());
        output.setProjection(dem->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing sediment transport...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;

                if (hasND && (demData[idx] == noData || accData[idx] == noData)) {
                    out[idx] = noData;
                    continue;
                }

                // Compute local slope using Horn's method
                auto val = [&](int rr, int cc) -> double {
                    rr = std::max(0, std::min(rows - 1, rr));
                    cc = std::max(0, std::min(cols - 1, cc));
                    size_t i = static_cast<size_t>(rr) * cols + cc;
                    if (hasND && demData[i] == noData) return demData[idx];
                    return demData[i];
                };

                double dzdx = ((val(r-1,c+1) + 2*val(r,c+1) + val(r+1,c+1)) -
                               (val(r-1,c-1) + 2*val(r,c-1) + val(r+1,c-1))) / (8.0 * cs);
                double dzdy = ((val(r+1,c-1) + 2*val(r+1,c) + val(r+1,c+1)) -
                               (val(r-1,c-1) + 2*val(r-1,c) + val(r-1,c+1))) / (8.0 * cs);

                double slope = std::sqrt(dzdx * dzdx + dzdy * dzdy);
                double area = accData[idx] * cs; // specific contributing area

                // Transport capacity: T = K * A^m * S^n
                double T = K * std::pow(std::max(area, 0.0), mExp) * std::pow(slope, nExp);
                out[idx] = T;
            }
            if (r % 100 == 0) reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SedimentationModule)

} // namespace aplaceholder
