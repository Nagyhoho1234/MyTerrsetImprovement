#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class PsdIdrisiModule : public Module {
public:
    QString name() const override { return "PSDIDRISI"; }
    QString description() const override {
        return "Import PSD (Power Spectral Density) data to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input PSD file"),
            ParameterDef::output("output", "Output TerrSet raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read PSD input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(PsdIdrisiModule)

} // namespace aplaceholder
