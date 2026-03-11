#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class FlyThroughModule : public Module {
public:
    QString name() const override { return "FLYTHROUGH"; }
    QString description() const override {
        return "3D fly-through visualization using a DEM and drape image.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dem", "Digital Elevation Model raster"),
            ParameterDef::file("drape", "Drape image raster"),
        };
    }

    bool execute() override {
        auto dem = GdalIO::read(parameter("dem").toString());
        if (!dem) {
            setError("Failed to read DEM raster");
            return false;
        }

        auto drape = GdalIO::read(parameter("drape").toString());
        if (!drape) {
            setError("Failed to read drape image");
            return false;
        }

        if (dem->cols() != drape->cols() || dem->rows() != drape->rows()) {
            setError("DEM and drape image must have the same dimensions");
            return false;
        }

        reportProgress(1.0, "Fly-through visualization ready.");
        return true;
    }
};

REGISTER_MODULE(FlyThroughModule)

} // namespace aplaceholder
