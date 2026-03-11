#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class StatIdrisModule : public Module {
public:
    QString name() const override { return "STATIDRIS"; }
    QString description() const override {
        return "Convert between Statistica format and TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input file"),
            ParameterDef::output("output", "Output file"),
            ParameterDef::combo("direction", "Direction",
                {"Statistica to TerrSet", "TerrSet to Statistica"}, 0),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input file"); return false; }

        int dir = parameter("direction").toInt();
        QString driver = (dir == 0) ? "RST" : "AAIGrid";

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), driver);
    }
};

REGISTER_MODULE(StatIdrisModule)

} // namespace aplaceholder
