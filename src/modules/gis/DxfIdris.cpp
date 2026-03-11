#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class DxfIdrisModule : public Module {
public:
    QString name() const override { return "DXFIDRIS"; }
    QString description() const override {
        return "Convert between DXF vector format and TerrSet vector format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input file"),
            ParameterDef::output("output", "Output file"),
            ParameterDef::combo("direction", "Direction",
                {"DXF to TerrSet", "TerrSet to DXF"}, 0),
        };
    }

    bool execute() override {
        // DXF is a vector format; use GDAL/OGR vector path
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input file"); return false; }

        int dir = parameter("direction").toInt();
        QString driver = (dir == 0) ? "RST" : "DXF";

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), driver);
    }
};

REGISTER_MODULE(DxfIdrisModule)

} // namespace aplaceholder
