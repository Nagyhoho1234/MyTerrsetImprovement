#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class NetCdfModule : public Module {
public:
    QString name() const override { return "NETCDF"; }
    QString description() const override {
        return "Import NetCDF data to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input NetCDF file"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::integer("band", "Band/variable to import", 1, 1, 9999,
                "Band or variable index to extract"),
        };
    }

    bool execute() override {
        int band = parameter("band").toInt();
        auto raster = GdalIO::readBand(parameter("input").toString(), band);
        if (!raster) { setError("Failed to read NetCDF input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(NetCdfModule)

} // namespace aplaceholder
