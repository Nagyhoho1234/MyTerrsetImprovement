#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class FuzzyModule : public Module {
public:
    QString name() const override { return "FUZZY"; }
    QString description() const override {
        return "Standardizes factor images into fuzzy set membership functions, "
               "transforming raw criterion values into a 0-1 scale representing "
               "degree of membership in a decision set.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output membership image"),
            ParameterDef::combo("function", "Membership function type",
                {"sigmoidal-inc", "sigmoidal-dec", "j-shaped-inc", "j-shaped-dec",
                 "linear-inc", "linear-dec"}, 0,
                "Type and direction of the fuzzy membership function"),
            ParameterDef::real("control_a", "Control point A (inflection low)",
                0.0, -999999, 999999,
                "Lower inflection point where membership begins to rise (or ends falling)"),
            ParameterDef::real("control_b", "Control point B (inflection high)",
                1.0, -999999, 999999,
                "Upper inflection point where membership reaches 1 (or begins falling)"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = input->cols(), rows = input->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int funcType = parameter("function").toInt();
        double a = parameter("control_a").toDouble();
        double b = parameter("control_b").toDouble();

        if (a >= b) {
            setError("Control point A must be less than control point B");
            return false;
        }

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(input->noDataValue());

        const auto& inData = input->data(0);
        auto& out = output.data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();

        reportProgress(0.0, "Computing fuzzy membership...");

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && inData[i] == noData) {
                out[i] = noData;
                continue;
            }

            double x = inData[i];
            double membership = 0.0;

            switch (funcType) {
                case 0: // sigmoidal-inc
                    membership = sigmoidalInc(x, a, b);
                    break;
                case 1: // sigmoidal-dec
                    membership = sigmoidalDec(x, a, b);
                    break;
                case 2: // j-shaped-inc
                    membership = jShapedInc(x, a, b);
                    break;
                case 3: // j-shaped-dec
                    membership = jShapedDec(x, a, b);
                    break;
                case 4: // linear-inc
                    membership = linearInc(x, a, b);
                    break;
                case 5: // linear-dec
                    membership = linearDec(x, a, b);
                    break;
            }

            out[i] = membership;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }

private:
    // Sigmoidal monotonically increasing: cos^2 based
    // Below a -> 0, above b -> 1, between a and b -> cos^2 transition
    double sigmoidalInc(double x, double a, double b) const {
        if (x <= a) return 0.0;
        if (x >= b) return 1.0;
        double t = (b - x) / (b - a);
        double cosVal = std::cos(t * M_PI / 2.0);
        return 1.0 - cosVal * cosVal;
    }

    // Sigmoidal monotonically decreasing
    double sigmoidalDec(double x, double a, double b) const {
        if (x <= a) return 1.0;
        if (x >= b) return 0.0;
        double t = (x - a) / (b - a);
        double cosVal = std::cos(t * M_PI / 2.0);
        return cosVal * cosVal;
    }

    // J-shaped monotonically increasing: 1/(1+((x-b)/(b-a))^2)
    // Below a -> low membership approaching 0, at b -> 1
    double jShapedInc(double x, double a, double b) const {
        if (x >= b) return 1.0;
        double ratio = (x - b) / (b - a);
        return 1.0 / (1.0 + ratio * ratio);
    }

    // J-shaped monotonically decreasing
    double jShapedDec(double x, double a, double b) const {
        if (x <= a) return 1.0;
        double ratio = (x - a) / (b - a);
        return 1.0 / (1.0 + ratio * ratio);
    }

    // Linear monotonically increasing
    double linearInc(double x, double a, double b) const {
        if (x <= a) return 0.0;
        if (x >= b) return 1.0;
        return (x - a) / (b - a);
    }

    // Linear monotonically decreasing
    double linearDec(double x, double a, double b) const {
        if (x <= a) return 1.0;
        if (x >= b) return 0.0;
        return (b - x) / (b - a);
    }
};

REGISTER_MODULE(FuzzyModule)

} // namespace aplaceholder
