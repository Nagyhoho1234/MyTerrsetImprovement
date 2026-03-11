#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class SstIdrisModule : public Module {
public:
    QString name() const override { return "SSTIDRIS"; }
    QString description() const override {
        return "Import ASCII grid or spreadsheet data to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input ASCII grid/spreadsheet file"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::combo("format", "Input format",
                {"ASCII Grid", "Spreadsheet (space-delimited)",
                 "Spreadsheet (tab-delimited)", "Spreadsheet (comma-delimited)"}, 0),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input file"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(SstIdrisModule)

} // namespace aplaceholder
