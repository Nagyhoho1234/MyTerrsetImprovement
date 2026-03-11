#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class SacIdrisModule : public Module {
public:
    QString name() const override { return "SACIDRIS"; }
    QString description() const override {
        return "Import SAC (Argentine Space Agency) satellite data to TerrSet format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input SAC file"),
            ParameterDef::output("output", "Output TerrSet raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read SAC input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(SacIdrisModule)

} // namespace aplaceholder
