#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class DlgModule : public Module {
public:
    QString name() const override { return "DLG"; }
    QString description() const override {
        return "Import USGS DLG (Digital Line Graph) format via GDAL.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input USGS DLG file"),
            ParameterDef::output("output", "Output raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input DLG file");
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

REGISTER_MODULE(DlgModule)

} // namespace aplaceholder
