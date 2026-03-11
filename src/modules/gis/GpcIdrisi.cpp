#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class GpcIdrisiModule : public Module {
public:
    QString name() const override { return "GPCIDRISI"; }
    QString description() const override {
        return "Import GPC ground control point data to TerrSet format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input GPC file"),
            ParameterDef::output("output", "Output TerrSet file"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read GPC input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(GpcIdrisiModule)

} // namespace aplaceholder
