#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class XyzIdrisModule : public Module {
public:
    QString name() const override { return "XYZIDRIS"; }
    QString description() const override {
        return "Import or export ASCII XYZ point data to/from TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input file"),
            ParameterDef::output("output", "Output file"),
            ParameterDef::combo("direction", "Direction",
                {"XYZ to TerrSet", "TerrSet to XYZ"}, 0),
        };
    }

    bool execute() override {
        // GDAL XYZ driver handles ASCII XYZ grids
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input file"); return false; }

        int dir = parameter("direction").toInt();
        QString driver = (dir == 0) ? "RST" : "XYZ";

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), driver);
    }
};

REGISTER_MODULE(XyzIdrisModule)

} // namespace aplaceholder
