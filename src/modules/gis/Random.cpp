#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <random>

namespace aplaceholder {

class RandomModule : public Module {
public:
    QString name() const override { return "RANDOM"; }
    QString description() const override {
        return "Generates a raster with random values from a specified distribution. "
               "Supports uniform (rectilinear) and normal distributions. Used for "
               "Monte Carlo simulation of measurement error propagation.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("reference", "Reference raster (for extent/dimensions)"),
            ParameterDef::output("output", "Output random image"),
            ParameterDef::combo("distribution", "Distribution model",
                {"Uniform", "Normal"}, 0,
                "Statistical distribution for random values"),
            ParameterDef::separator("Uniform distribution parameters"),
            ParameterDef::real("min_value", "Minimum value",
                0.0, -999999, 999999,
                "Minimum value for uniform distribution"),
            ParameterDef::real("max_value", "Maximum value",
                1.0, -999999, 999999,
                "Maximum value for uniform distribution"),
            ParameterDef::separator("Normal distribution parameters"),
            ParameterDef::real("mean", "Mean",
                0.0, -999999, 999999,
                "Mean for normal distribution"),
            ParameterDef::real("stddev", "Standard deviation",
                1.0, 0.0, 999999,
                "Standard deviation (RMS) for normal distribution"),
        };
    }

    bool execute() override {
        auto ref = GdalIO::read(parameter("reference").toString());
        if (!ref) {
            setError("Failed to read reference raster");
            return false;
        }

        int cols = ref->cols(), rows = ref->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int distType = parameter("distribution").toInt();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(ref->geoTransform());
        output.setProjection(ref->projection());

        // Respect nodata from reference raster
        const auto& refData = ref->data(0);
        double refNoData = ref->noDataValue();
        bool refHasND = ref->hasNoData();
        if (refHasND)
            output.setNoDataValue(refNoData);

        auto& out = output.data(0);

        std::mt19937_64 rng(std::random_device{}());

        reportProgress(0.1, "Generating random values...");

        if (distType == 0) {
            // Uniform distribution
            double minVal = parameter("min_value").toDouble();
            double maxVal = parameter("max_value").toDouble();

            if (minVal > maxVal) {
                setError("Minimum value must be less than or equal to maximum value");
                return false;
            }

            std::uniform_real_distribution<double> dist(minVal, maxVal);

            for (int64_t i = 0; i < total; ++i) {
                if (refHasND && refData[i] == refNoData) {
                    out[i] = refNoData;
                } else {
                    out[i] = dist(rng);
                }

                if (i % 1000000 == 0)
                    reportProgress(0.1 + 0.85 * static_cast<double>(i) / total);
            }
        } else {
            // Normal distribution
            double mean = parameter("mean").toDouble();
            double stddev = parameter("stddev").toDouble();

            if (stddev < 0) {
                setError("Standard deviation must be non-negative");
                return false;
            }

            std::normal_distribution<double> dist(mean, stddev);

            for (int64_t i = 0; i < total; ++i) {
                if (refHasND && refData[i] == refNoData) {
                    out[i] = refNoData;
                } else {
                    out[i] = dist(rng);
                }

                if (i % 1000000 == 0)
                    reportProgress(0.1 + 0.85 * static_cast<double>(i) / total);
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(RandomModule)

} // namespace aplaceholder
