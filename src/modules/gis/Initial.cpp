#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class InitialModule : public Module {
public:
    QString name() const override { return "INITIAL"; }
    QString description() const override {
        return "Initialize a new raster filled with a constant value.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::integer("rows", "Number of rows", 100, 1, 999999,
                "Number of rows in the output raster"),
            ParameterDef::integer("cols", "Number of columns", 100, 1, 999999,
                "Number of columns in the output raster"),
            ParameterDef::real("value", "Fill value", 0.0, -999999.0, 999999.0,
                "Constant value to fill the raster with"),
            ParameterDef::combo("data_type", "Data type",
                {"Byte", "Int16", "Int32", "Float32", "Float64"}, 3,
                "Data type for the output raster"),
            ParameterDef::output("output", "Output raster"),
        };
    }

    bool execute() override {
        int rows = parameter("rows").toInt();
        int cols = parameter("cols").toInt();
        double fillValue = parameter("value").toDouble();
        int typeIdx = parameter("data_type").toInt();

        DataType targetType;
        switch (typeIdx) {
            case 0: targetType = DataType::Byte;    break;
            case 1: targetType = DataType::Int16;   break;
            case 2: targetType = DataType::Int32;   break;
            case 3: targetType = DataType::Float32; break;
            case 4: targetType = DataType::Float64; break;
            default: targetType = DataType::Float64; break;
        }

        Raster output(cols, rows, 1, targetType);

        auto& dst = output.data(0);
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int64_t i = 0; i < total; ++i) {
            dst[i] = fillValue;

            if (i % 100000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(InitialModule)

} // namespace aplaceholder
