#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class SpotModule : public Module {
public:
    QString name() const override { return "SPOT"; }
    QString description() const override {
        return "Import SPOT satellite imagery to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input SPOT file"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::boolean("all_bands", "Import all bands", true),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read SPOT input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(SpotModule)

} // namespace aplaceholder
