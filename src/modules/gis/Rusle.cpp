#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class RusleModule : public Module {
public:
    QString name() const override { return "RUSLE"; }
    QString description() const override {
        return "RUSLE soil erosion model: A = R * K * LS * C * P. "
               "Multiplies five factor rasters to estimate annual soil loss.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("r_factor", "R factor (rainfall erosivity)"),
            ParameterDef::file("k_factor", "K factor (soil erodibility)"),
            ParameterDef::file("ls_factor", "LS factor (slope length/steepness)"),
            ParameterDef::file("c_factor", "C factor (cover management)"),
            ParameterDef::file("p_factor", "P factor (support practice)"),
            ParameterDef::output("output", "Output soil loss raster"),
        };
    }

    bool execute() override {
        auto rFac = GdalIO::read(parameter("r_factor").toString());
        auto kFac = GdalIO::read(parameter("k_factor").toString());
        auto lsFac = GdalIO::read(parameter("ls_factor").toString());
        auto cFac = GdalIO::read(parameter("c_factor").toString());
        auto pFac = GdalIO::read(parameter("p_factor").toString());

        if (!rFac || !kFac || !lsFac || !cFac || !pFac) {
            setError("Failed to read one or more input rasters");
            return false;
        }

        int cols = rFac->cols(), rows = rFac->rows();
        if (kFac->cols() != cols || kFac->rows() != rows ||
            lsFac->cols() != cols || lsFac->rows() != rows ||
            cFac->cols() != cols || cFac->rows() != rows ||
            pFac->cols() != cols || pFac->rows() != rows) {
            setError("All input rasters must have the same dimensions");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& R = rFac->data(0);
        const auto& K = kFac->data(0);
        const auto& LS = lsFac->data(0);
        const auto& C = cFac->data(0);
        const auto& P = pFac->data(0);

        double nd = -9999.0;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(rFac->geoTransform());
        output.setProjection(rFac->projection());
        output.setNoDataValue(nd);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing RUSLE...");
        for (int64_t i = 0; i < total; ++i) {
            // Check nodata for any input
            bool isND = false;
            if (rFac->hasNoData() && R[i] == rFac->noDataValue()) isND = true;
            if (kFac->hasNoData() && K[i] == kFac->noDataValue()) isND = true;
            if (lsFac->hasNoData() && LS[i] == lsFac->noDataValue()) isND = true;
            if (cFac->hasNoData() && C[i] == cFac->noDataValue()) isND = true;
            if (pFac->hasNoData() && P[i] == pFac->noDataValue()) isND = true;

            if (isND) {
                out[i] = nd;
            } else {
                out[i] = R[i] * K[i] * LS[i] * C[i] * P[i];
            }

            if (i % 500000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(RusleModule)

} // namespace aplaceholder
