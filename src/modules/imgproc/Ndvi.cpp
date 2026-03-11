#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class NdviModule : public Module {
public:
    QString name() const override { return "NDVI"; }
    QString description() const override {
        return "Normalized Difference Vegetation Index calculation. Computes "
               "NDVI = (NIR - Red) / (NIR + Red) from two input bands.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("nir_band", "NIR band file"),
            ParameterDef::file("red_band", "Red band file"),
            ParameterDef::output("output", "Output NDVI image"),
        };
    }

    bool execute() override {
        auto nir = GdalIO::read(parameter("nir_band").toString());
        auto red = GdalIO::read(parameter("red_band").toString());
        if (!nir || !red) {
            setError("Failed to read input band(s)");
            return false;
        }

        if (nir->cols() != red->cols() || nir->rows() != red->rows()) {
            setError("NIR and Red bands must have the same dimensions");
            return false;
        }

        int cols = nir->cols(), rows = nir->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(nir->geoTransform());
        output.setProjection(nir->projection());
        double noData = -9999.0;
        output.setNoDataValue(noData);

        const auto& nirData = nir->data(0);
        const auto& redData = red->data(0);
        auto& out = output.data(0);

        bool nirHasND = nir->hasNoData();
        bool redHasND = red->hasNoData();
        double nirND = nir->noDataValue();
        double redND = red->noDataValue();

        for (int64_t i = 0; i < total; ++i) {
            double n = nirData[i];
            double r = redData[i];

            // Check for nodata in either input
            if ((nirHasND && n == nirND) || (redHasND && r == redND)) {
                out[i] = noData;
                continue;
            }

            double sum = n + r;
            if (std::abs(sum) < 1e-10) {
                // Division by zero: both bands are zero or cancel out
                out[i] = noData;
            } else {
                double ndvi = (n - r) / sum;
                // Clamp to valid NDVI range
                out[i] = std::clamp(ndvi, -1.0, 1.0);
            }

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(NdviModule)

} // namespace aplaceholder
