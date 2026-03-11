#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class VxpModule : public Module {
public:
    QString name() const override { return "VXP"; }
    QString description() const override {
        return "IDRISI Vector Export - export vector data to VXP format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input vector data"),
            ParameterDef::output("output", "Output VXP file"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input vector data");
            return false;
        }

        reportProgress(0.5, "Exporting to VXP format...");
        if (!GdalIO::write(*raster, parameter("output").toString())) {
            setError("Failed to write VXP output");
            return false;
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(VxpModule)

} // namespace aplaceholder
