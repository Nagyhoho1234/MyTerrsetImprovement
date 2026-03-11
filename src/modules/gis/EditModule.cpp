#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class EditModule : public Module {
public:
    QString name() const override { return "EDIT"; }
    QString description() const override {
        return "Simple raster value editor. Modifies a single cell value.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::integer("row", "Row index (0-based)", 0),
            ParameterDef::integer("col", "Column index (0-based)", 0),
            ParameterDef::real("new_value", "New cell value", 0.0),
            ParameterDef::output("output", "Output raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int row = parameter("row").toInt();
        int col = parameter("col").toInt();
        double newVal = parameter("new_value").toDouble();

        int cols = raster->cols(), rows = raster->rows();
        if (row < 0 || row >= rows || col < 0 || col >= cols) {
            setError(QString("Row/col (%1,%2) out of bounds (%3 x %4)")
                         .arg(row).arg(col).arg(rows).arg(cols));
            return false;
        }

        int64_t idx = static_cast<int64_t>(row) * cols + col;
        auto& data = raster->data(0);
        data[idx] = newVal;

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString());
    }
};

REGISTER_MODULE(EditModule)

} // namespace aplaceholder
