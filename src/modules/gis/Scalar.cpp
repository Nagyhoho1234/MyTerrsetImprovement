#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ScalarModule : public Module {
public:
    QString name() const override { return "SCALAR"; }
    QString description() const override {
        return "Applies a scalar (constant) mathematical operation to every pixel "
               "in a raster image. Supports add, subtract, multiply, divide, power, "
               "log, exp, sqrt, and abs operations.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output image"),
            ParameterDef::combo("operation", "Operation",
                {"Add", "Subtract", "Multiply", "Divide", "Power",
                 "Log", "Exp", "Sqrt", "Abs"}, 0,
                "Mathematical operation to apply"),
            ParameterDef::real("scalar", "Scalar value",
                1.0, -999999, 999999,
                "Constant value for the operation"),
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
        int op = parameter("operation").toInt();
        double scalar = parameter("scalar").toDouble();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(input->noDataValue());

        const auto& inData = input->data(0);
        auto& out = output.data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();

        reportProgress(0.0, "Applying scalar operation...");

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && inData[i] == noData) {
                out[i] = noData;
                continue;
            }

            double x = inData[i];
            double result = noData;

            switch (op) {
                case 0: // Add
                    result = x + scalar;
                    break;
                case 1: // Subtract
                    result = x - scalar;
                    break;
                case 2: // Multiply
                    result = x * scalar;
                    break;
                case 3: // Divide
                    if (scalar != 0.0)
                        result = x / scalar;
                    else
                        result = noData;
                    break;
                case 4: // Power
                    result = std::pow(x, scalar);
                    break;
                case 5: // Log
                    if (x > 0.0)
                        result = std::log(x) / std::log(scalar > 0.0 && scalar != 1.0 ? scalar : M_E);
                    else
                        result = noData;
                    break;
                case 6: // Exp
                    result = std::exp(x * scalar);
                    break;
                case 7: // Sqrt
                    if (x >= 0.0)
                        result = std::sqrt(x);
                    else
                        result = noData;
                    break;
                case 8: // Abs
                    result = std::abs(x);
                    break;
            }

            out[i] = result;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ScalarModule)

} // namespace aplaceholder
