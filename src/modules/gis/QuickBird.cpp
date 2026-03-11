#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class QuickBirdModule : public Module {
public:
    QString name() const override { return "QUICKBIRD"; }
    QString description() const override {
        return "Import QuickBird satellite imagery to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input QuickBird file"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::boolean("apply_rpc", "Apply RPC correction", false),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read QuickBird input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(QuickBirdModule)

} // namespace aplaceholder
