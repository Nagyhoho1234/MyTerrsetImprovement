#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class GenericRasterModule : public Module {
public:
    QString name() const override { return "GENERICRASTER"; }
    QString description() const override {
        return "Import generic binary raster (BIL, BIP, BSQ) to TerrSet format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster file"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::combo("interleave", "Interleave type",
                {"BIL (Band Interleaved by Line)",
                 "BIP (Band Interleaved by Pixel)",
                 "BSQ (Band Sequential)"}, 0),
            ParameterDef::integer("cols", "Number of columns", 0, 1, 999999),
            ParameterDef::integer("rows", "Number of rows", 0, 1, 999999),
            ParameterDef::integer("bands", "Number of bands", 1, 1, 9999),
            ParameterDef::combo("data_type", "Data type",
                {"Byte", "Int16", "Int32", "Float32", "Float64"}, 0),
            ParameterDef::integer("header_bytes", "Header bytes to skip", 0, 0, 999999),
        };
    }

    bool execute() override {
        // GDAL can read generic binary rasters via the EHdr driver
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input raster"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(GenericRasterModule)

} // namespace aplaceholder
