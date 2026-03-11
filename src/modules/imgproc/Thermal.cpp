#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ThermalModule : public Module {
public:
    QString name() const override { return "THERMAL"; }
    QString description() const override {
        return "Thermal band processing. Converts raw DN values to brightness "
               "temperature using calibration coefficients. Uses the formula "
               "T = K2 / ln(K1 / radiance + 1) where radiance = gain * DN + offset.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_thermal", "Input thermal band file"),
            ParameterDef::real("gain", "Radiometric rescaling gain (multiplicative factor)", 0.0),
            ParameterDef::real("offset", "Radiometric rescaling offset (additive factor)", 0.0),
            ParameterDef::real("k1", "Calibration constant K1 (in W/(m^2 sr um))", 0.0),
            ParameterDef::real("k2", "Calibration constant K2 (in Kelvin)", 0.0),
            ParameterDef::output("output", "Output brightness temperature image"),
        };
    }

    bool execute() override {
        auto thermal = GdalIO::read(parameter("input_thermal").toString());
        if (!thermal) {
            setError("Failed to read input thermal band");
            return false;
        }

        double gain   = parameter("gain").toDouble();
        double offset = parameter("offset").toDouble();
        double k1     = parameter("k1").toDouble();
        double k2     = parameter("k2").toDouble();

        if (k1 <= 0.0 || k2 <= 0.0) {
            setError("K1 and K2 constants must be positive values");
            return false;
        }

        int cols = thermal->cols(), rows = thermal->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(thermal->geoTransform());
        output.setProjection(thermal->projection());
        double noData = -9999.0;
        output.setNoDataValue(noData);

        const auto& thermalData = thermal->data(0);
        auto& out = output.data(0);

        bool hasND = thermal->hasNoData();
        double nd = thermal->noDataValue();

        for (int64_t i = 0; i < total; ++i) {
            double dn = thermalData[i];

            // Check for nodata
            if (hasND && dn == nd) {
                out[i] = noData;
                continue;
            }

            // Convert DN to at-sensor spectral radiance
            double radiance = gain * dn + offset;

            if (radiance <= 0.0) {
                // Invalid radiance value
                out[i] = noData;
                continue;
            }

            // Convert radiance to brightness temperature
            // T = K2 / ln(K1 / radiance + 1)
            double temp = k2 / std::log(k1 / radiance + 1.0);
            out[i] = temp;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ThermalModule)

} // namespace aplaceholder
