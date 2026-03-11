#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class HdfEosModule : public Module {
public:
    QString name() const override { return "HDFEOS"; }
    QString description() const override {
        return "Import HDF or HDF-EOS data to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input HDF/HDF-EOS file"),
            ParameterDef::output("output", "Output TerrSet raster"),
            ParameterDef::integer("band", "Band/dataset to import", 1, 1, 9999,
                "Band or subdataset index to extract"),
        };
    }

    bool execute() override {
        int band = parameter("band").toInt();
        auto raster = GdalIO::readBand(parameter("input").toString(), band);
        if (!raster) { setError("Failed to read HDF input"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(HdfEosModule)

} // namespace aplaceholder
