#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class GacpIdrisiModule : public Module {
public:
    QString name() const override { return "GACPIDRISI"; }
    QString description() const override {
        return "Import GACP (Global Aerosol Climatology Project) data to TerrSet format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input GACP file"),
            ParameterDef::output("output", "Output TerrSet raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read GACP input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(GacpIdrisiModule)

} // namespace aplaceholder
