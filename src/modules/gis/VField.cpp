#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class VFieldModule : public Module {
public:
    QString name() const override { return "VFIELD"; }
    QString description() const override {
        return "Vector field display from U and V component rasters.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("u_component", "U (east-west) component raster"),
            ParameterDef::file("v_component", "V (north-south) component raster"),
        };
    }

    bool execute() override {
        auto uRaster = GdalIO::read(parameter("u_component").toString());
        if (!uRaster) {
            setError("Failed to read U component raster");
            return false;
        }

        auto vRaster = GdalIO::read(parameter("v_component").toString());
        if (!vRaster) {
            setError("Failed to read V component raster");
            return false;
        }

        if (uRaster->cols() != vRaster->cols() || uRaster->rows() != vRaster->rows()) {
            setError("U and V component rasters must have the same dimensions");
            return false;
        }

        reportProgress(1.0, "Vector field display ready.");
        return true;
    }
};

REGISTER_MODULE(VFieldModule)

} // namespace aplaceholder
