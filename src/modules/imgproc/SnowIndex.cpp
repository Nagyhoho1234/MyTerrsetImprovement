#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class SnowIndexModule : public Module {
public:
    QString name() const override { return "SNOWINDEX"; }
    QString description() const override {
        return "Normalized Difference Snow Index calculation. Computes "
               "NDSI = (Green - SWIR) / (Green + SWIR) from two input bands.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("green_band", "Green band file"),
            ParameterDef::file("swir_band", "SWIR band file"),
            ParameterDef::output("output", "Output NDSI image"),
        };
    }

    bool execute() override {
        auto green = GdalIO::read(parameter("green_band").toString());
        auto swir  = GdalIO::read(parameter("swir_band").toString());
        if (!green || !swir) {
            setError("Failed to read input band(s)");
            return false;
        }

        if (green->cols() != swir->cols() || green->rows() != swir->rows()) {
            setError("Green and SWIR bands must have the same dimensions");
            return false;
        }

        int cols = green->cols(), rows = green->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(green->geoTransform());
        output.setProjection(green->projection());
        double noData = -9999.0;
        output.setNoDataValue(noData);

        const auto& greenData = green->data(0);
        const auto& swirData  = swir->data(0);
        auto& out = output.data(0);

        bool greenHasND = green->hasNoData();
        bool swirHasND  = swir->hasNoData();
        double greenND  = green->noDataValue();
        double swirND   = swir->noDataValue();

        for (int64_t i = 0; i < total; ++i) {
            double g = greenData[i];
            double s = swirData[i];

            // Check for nodata in either input
            if ((greenHasND && g == greenND) || (swirHasND && s == swirND)) {
                out[i] = noData;
                continue;
            }

            double sum = g + s;
            if (std::abs(sum) < 1e-10) {
                // Division by zero: both bands are zero or cancel out
                out[i] = noData;
            } else {
                double ndsi = (g - s) / sum;
                // Clamp to valid NDSI range
                out[i] = std::clamp(ndsi, -1.0, 1.0);
            }

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SnowIndexModule)

} // namespace aplaceholder
