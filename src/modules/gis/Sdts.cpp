#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class SdtsModule : public Module {
public:
    QString name() const override { return "SDTS"; }
    QString description() const override {
        return "Import SDTS (Spatial Data Transfer Standard) format via GDAL.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input SDTS file"),
            ParameterDef::output("output", "Output raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input SDTS file");
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

REGISTER_MODULE(SdtsModule)

} // namespace aplaceholder
