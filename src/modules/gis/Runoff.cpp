#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class RunoffModule : public Module {
public:
    QString name() const override { return "RUNOFF"; }
    QString description() const override {
        return "Simple runoff modeling using the SCS Curve Number method. "
               "Q = (P - 0.2*S)^2 / (P + 0.8*S) where S = (1000/CN - 10) * 25.4.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("precipitation_raster", "Precipitation raster (mm)"),
            ParameterDef::file("curve_number_raster", "Curve number raster"),
            ParameterDef::output("output", "Output runoff raster (mm)"),
        };
    }

    bool execute() override {
        auto precip = GdalIO::read(parameter("precipitation_raster").toString());
        auto cn = GdalIO::read(parameter("curve_number_raster").toString());
        if (!precip || !cn) { setError("Failed to read input rasters"); return false; }

        if (precip->cols() != cn->cols() || precip->rows() != cn->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = precip->cols(), rows = precip->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& pData = precip->data(0);
        const auto& cnData = cn->data(0);

        double pND = precip->noDataValue();
        bool pHasND = precip->hasNoData();
        double cnND = cn->noDataValue();
        bool cnHasND = cn->hasNoData();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(precip->geoTransform());
        output.setProjection(precip->projection());
        output.setNoDataValue(-9999);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing runoff...");
        for (int64_t i = 0; i < total; ++i) {
            if ((pHasND && pData[i] == pND) || (cnHasND && cnData[i] == cnND)) {
                out[i] = -9999;
                continue;
            }

            double P = pData[i];
            double CN = cnData[i];

            if (CN <= 0 || CN > 100 || P <= 0) {
                out[i] = 0.0;
                continue;
            }

            double S = (1000.0 / CN - 10.0) * 25.4;
            double Ia = 0.2 * S; // initial abstraction

            if (P <= Ia) {
                out[i] = 0.0;
            } else {
                out[i] = (P - Ia) * (P - Ia) / (P + 0.8 * S);
            }

            if (i % 500000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(RunoffModule)

} // namespace aplaceholder
