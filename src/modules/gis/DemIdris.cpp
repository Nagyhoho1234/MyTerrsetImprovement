#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class DemIdrisModule : public Module {
public:
    QString name() const override { return "DEMIDRIS"; }
    QString description() const override {
        return "Import USGS DEM format to IDRISI raster via GDAL.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input USGS DEM file"),
            ParameterDef::output("output", "Output raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input USGS DEM file");
            return false;
        }

        reportProgress(0.5, "Writing output raster...");
        if (!GdalIO::write(*raster, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(DemIdrisModule)

} // namespace aplaceholder
