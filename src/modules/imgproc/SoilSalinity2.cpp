#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class SoilSalinityModule : public Module {
public:
    QString name() const override { return "SOILSALINITY"; }
    QString description() const override {
        return "Soil salinity index computation from multispectral bands (SI = sqrt(Red * NIR)).";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("red_band", "Red band input raster"),
            ParameterDef::file("nir_band", "NIR band input raster"),
            ParameterDef::output("output", "Output salinity index raster"),
        };
    }

    bool execute() override {
        auto red = GdalIO::read(parameter("red_band").toString());
        if (!red) { setError("Failed to read Red band raster"); return false; }

        auto nir = GdalIO::read(parameter("nir_band").toString());
        if (!nir) { setError("Failed to read NIR band raster"); return false; }

        int cols = red->cols(), rows = red->rows();
        if (cols != nir->cols() || rows != nir->rows()) {
            setError("Red and NIR bands must have the same dimensions");
            return false;
        }

        Raster output(cols, rows, 1, DataType::Float32);
        output.setGeoTransform(red->geoTransform());
        output.setProjection(red->projection());
        if (red->hasNoData())
            output.setNoDataValue(red->noDataValue());

        const auto& srcRed = red->data(0);
        const auto& srcNir = nir->data(0);
        auto& dst = output.data(0);
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int64_t i = 0; i < total; ++i) {
            double r = srcRed[i];
            double n = srcNir[i];
            double product = r * n;
            dst[i] = (product >= 0.0) ? std::sqrt(product) : 0.0;

            if (i % 100000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SoilSalinityModule)

} // namespace aplaceholder
