#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class PClassModule : public Module {
public:
    QString name() const override { return "PCLASS"; }
    QString description() const override {
        return "Probability classification. Evaluates the probability that each cell's "
               "true value exceeds a specified threshold, given measurement error (standard deviation).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image",
                "Quantitative raster image"),
            ParameterDef::file("error_raster", "Error (std deviation) image",
                "RMS error / standard deviation raster"),
            ParameterDef::real("threshold", "Threshold value", 0, -999999, 999999,
                "Critical value to evaluate against"),
            ParameterDef::output("output", "Output probability image"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        auto errorRst = GdalIO::read(parameter("error_raster").toString());
        if (!input || !errorRst) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = input->cols(), rows = input->rows();
        if (errorRst->cols() != cols || errorRst->rows() != rows) {
            setError("Input and error rasters must have the same dimensions");
            return false;
        }

        double threshold = parameter("threshold").toDouble();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(input->noDataValue());

        const auto& dIn = input->data(0);
        const auto& dErr = errorRst->data(0);
        auto& out = output.data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (dIn[i] == noData || dErr[i] == noData)) {
                out[i] = noData;
                continue;
            }

            double sd = dErr[i];
            if (sd <= 0.0) {
                // Zero error: deterministic result
                out[i] = (dIn[i] > threshold) ? 1.0 : 0.0;
            } else {
                // P(x > threshold) = 1 - Phi((threshold - value) / sd)
                // Using complementary error function: Phi(z) = 0.5 * erfc(-z / sqrt(2))
                double z = (threshold - dIn[i]) / sd;
                out[i] = 0.5 * std::erfc(z / std::sqrt(2.0));
            }

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PClassModule)

} // namespace aplaceholder
