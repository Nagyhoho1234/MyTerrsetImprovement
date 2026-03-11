#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class Ndvi3gModule : public Module {
public:
    QString name() const override { return "NDVI3G"; }
    QString description() const override {
        return "Import GIMMS NDVI3g data to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input NDVI3g file"),
            ParameterDef::output("output", "Output TerrSet raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read NDVI3g input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(Ndvi3gModule)

} // namespace aplaceholder
