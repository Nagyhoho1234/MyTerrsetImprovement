#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class WaterIndexModule : public Module {
public:
    QString name() const override { return "WATERINDEX"; }
    QString description() const override {
        return "NDWI (Normalized Difference Water Index): (Green - NIR) / (Green + NIR). "
               "Values near +1 indicate water, values near -1 indicate vegetation.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("green_band", "Green band raster"),
            ParameterDef::file("nir_band", "NIR band raster"),
            ParameterDef::output("output", "Output NDWI raster"),
        };
    }

    bool execute() override {
        auto green = GdalIO::read(parameter("green_band").toString());
        auto nir = GdalIO::read(parameter("nir_band").toString());
        if (!green || !nir) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = green->cols(), rows = green->rows();
        if (nir->cols() != cols || nir->rows() != rows) {
            setError("Green and NIR rasters must have the same dimensions");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& gData = green->data(0);
        const auto& nData = nir->data(0);
        double gND = green->noDataValue();
        bool gHasND = green->hasNoData();
        double nND = nir->noDataValue();
        bool nHasND = nir->hasNoData();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(green->geoTransform());
        output.setProjection(green->projection());
        output.setNoDataValue(-9999);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing NDWI...");
        for (int64_t i = 0; i < total; ++i) {
            if ((gHasND && gData[i] == gND) || (nHasND && nData[i] == nND)) {
                out[i] = -9999;
                continue;
            }

            double sum = gData[i] + nData[i];
            if (std::abs(sum) < 1e-10) {
                out[i] = 0.0;
            } else {
                out[i] = (gData[i] - nData[i]) / sum;
            }

            if (i % 500000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(WaterIndexModule)

} // namespace aplaceholder
