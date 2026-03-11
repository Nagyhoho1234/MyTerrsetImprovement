#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class SentinelModule : public Module {
public:
    QString name() const override { return "SENTINEL"; }
    QString description() const override {
        return "Import Sentinel satellite imagery to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input Sentinel file or SAFE directory"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::combo("sensor", "Sensor",
                {"Sentinel-1 (SAR)", "Sentinel-2 (MSI)", "Sentinel-3 (OLCI)"}, 1),
            ParameterDef::boolean("all_bands", "Import all bands", true),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read Sentinel input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(SentinelModule)

} // namespace aplaceholder
