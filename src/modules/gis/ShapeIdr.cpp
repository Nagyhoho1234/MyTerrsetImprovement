#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class ShapeIdrModule : public Module {
public:
    QString name() const override { return "SHAPEIDR"; }
    QString description() const override {
        return "Convert between ESRI Shapefile and TerrSet vector format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input file"),
            ParameterDef::output("output", "Output file"),
            ParameterDef::combo("direction", "Direction",
                {"Shapefile to TerrSet", "TerrSet to Shapefile"}, 0),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input file"); return false; }

        int dir = parameter("direction").toInt();
        QString driver = (dir == 0) ? "RST" : "ESRI Shapefile";

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), driver);
    }
};

REGISTER_MODULE(ShapeIdrModule)

} // namespace aplaceholder
