#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class DigitalGlobeModule : public Module {
public:
    QString name() const override { return "DIGITALGLOBE"; }
    QString description() const override {
        return "Import DigitalGlobe satellite imagery to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input DigitalGlobe file"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::boolean("apply_rpc", "Apply RPC correction", false),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read DigitalGlobe input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(DigitalGlobeModule)

} // namespace aplaceholder
