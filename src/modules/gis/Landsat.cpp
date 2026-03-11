#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class LandsatModule : public Module {
public:
    QString name() const override { return "LANDSAT"; }
    QString description() const override {
        return "Import Landsat data archive imagery to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input Landsat file or metadata"),
            ParameterDef::output("output", "Output TerrSet raster (base name)"),
            ParameterDef::combo("sensor", "Sensor type",
                {"Landsat 5 TM", "Landsat 7 ETM+", "Landsat 8 OLI/TIRS",
                 "Landsat 9 OLI-2/TIRS-2"}, 2),
            ParameterDef::boolean("all_bands", "Import all bands", true),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read Landsat input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(LandsatModule)

} // namespace aplaceholder
