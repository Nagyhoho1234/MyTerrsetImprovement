#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class ModisConvModule : public Module {
public:
    QString name() const override { return "MODISCONV"; }
    QString description() const override {
        return "Convert MODIS HDF data to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input MODIS HDF file"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::integer("subdataset", "Subdataset index", 1, 1, 999),
            ParameterDef::boolean("apply_scale", "Apply scale factor", true),
        };
    }

    bool execute() override {
        int band = parameter("subdataset").toInt();
        auto raster = GdalIO::readBand(parameter("input").toString(), band);
        if (!raster) { setError("Failed to read MODIS input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(ModisConvModule)

} // namespace aplaceholder
