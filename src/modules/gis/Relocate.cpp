#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class RelocateModule : public Module {
public:
    QString name() const override { return "RELOCATE"; }
    QString description() const override {
        return "Shift/relocate a raster by changing its geotransform origin.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::real("shift_x", "Shift X", 0.0, -1e15, 1e15,
                "Horizontal shift in map units"),
            ParameterDef::real("shift_y", "Shift Y", 0.0, -1e15, 1e15,
                "Vertical shift in map units"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        double shiftX = parameter("shift_x").toDouble();
        double shiftY = parameter("shift_y").toDouble();

        GeoTransform gt = raster->geoTransform();
        gt.originX += shiftX;
        gt.originY += shiftY;
        raster->setGeoTransform(gt);

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString());
    }
};

REGISTER_MODULE(RelocateModule)

} // namespace aplaceholder
