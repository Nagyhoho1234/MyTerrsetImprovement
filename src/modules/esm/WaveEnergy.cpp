#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class WaveEnergyModule : public Module {
public:
    QString name() const override { return "WAVE_ENERGY"; }
    QString description() const override {
        return "Wave energy potential estimation. Computes wave power density using the "
               "deep-water wave energy formula: E = (rho * g^2 * H^2 * T) / (64 * pi), "
               "where H is significant wave height (m) and T is wave period (s). "
               "Output is in watts per metre of wave crest (W/m).";
    }
    QString category() const override { return "Ecosystem Services"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("wave_height_raster", "Wave height raster (m)",
                "Raster of significant wave height Hs in metres"),
            ParameterDef::file("wave_period_raster", "Wave period raster (s)",
                "Raster of wave energy period Te in seconds"),
            ParameterDef::output("output", "Output wave energy raster (W/m)",
                "Continuous raster of wave power density in watts per metre"),
            ParameterDef::real("water_density", "Water density (kg/m^3)", 1025.0, 900.0, 1100.0,
                "Density of seawater. Default 1025 kg/m^3."),
        };
    }

    bool execute() override {
        // 1. Load input rasters
        auto waveHeight = GdalIO::read(parameter("wave_height_raster").toString());
        if (!waveHeight) {
            setError("Failed to read wave height raster.");
            return false;
        }

        auto wavePeriod = GdalIO::read(parameter("wave_period_raster").toString());
        if (!wavePeriod) {
            setError("Failed to read wave period raster.");
            return false;
        }

        reportProgress(0.10, "Input rasters loaded.");

        // 2. Validate dimensions
        int cols = waveHeight->cols();
        int rows = waveHeight->rows();
        int64_t total = waveHeight->cellCount();

        if (wavePeriod->cols() != cols || wavePeriod->rows() != rows) {
            setError(QString("Dimension mismatch: wave height is %1x%2 but wave period is %3x%4")
                         .arg(cols).arg(rows)
                         .arg(wavePeriod->cols()).arg(wavePeriod->rows()));
            return false;
        }

        // 3. Constants
        double rho = parameter("water_density").toDouble();
        const double g = 9.80665;           // gravitational acceleration (m/s^2)
        const double pi = 3.14159265358979;
        double coeff = (rho * g * g) / (64.0 * pi);

        // 4. NoData handling
        double noData = waveHeight->noDataValue();
        bool hasND_H = waveHeight->hasNoData();
        bool hasND_T = wavePeriod->hasNoData();
        double ndH = waveHeight->noDataValue();
        double ndT = wavePeriod->noDataValue();

        // 5. Compute wave energy: E = coeff * H^2 * T
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(waveHeight->geoTransform());
        output.setProjection(waveHeight->projection());
        output.setNoDataValue(noData);

        const auto& hData = waveHeight->data(0);
        const auto& tData = wavePeriod->data(0);
        auto& outData = output.data(0);

        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            if (hasND_H && hData[i] == ndH) skip = true;
            if (hasND_T && tData[i] == ndT) skip = true;

            if (skip) {
                outData[i] = noData;
                continue;
            }

            double H = hData[i];
            double T = tData[i];

            outData[i] = coeff * H * H * T;
            ++validCount;

            if (i % 500000 == 0)
                reportProgress(0.10 + 0.80 * static_cast<double>(i) / total);
        }

        reportProgress(0.90, "Writing output...");

        // 6. Write output
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        auto stats = output.computeStats(0);
        reportProgress(1.0,
            QString("Done. Wave energy — min: %1 W/m, max: %2 W/m, mean: %3 W/m")
                .arg(stats.min, 0, 'f', 2)
                .arg(stats.max, 0, 'f', 2)
                .arg(stats.mean, 0, 'f', 2));

        return true;
    }
};

REGISTER_MODULE(WaveEnergyModule)

} // namespace aplaceholder
