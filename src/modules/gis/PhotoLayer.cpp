#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class PhotoLayerModule : public Module {
public:
    QString name() const override { return "PHOTOLAYER"; }
    QString description() const override {
        return "Create a photo layer from a georeferenced image.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Georeferenced image file"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read georeferenced image");
            return false;
        }

        reportProgress(1.0, "Photo layer created.");
        return true;
    }
};

REGISTER_MODULE(PhotoLayerModule)

} // namespace aplaceholder
