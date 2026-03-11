#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class ModisQcModule : public Module {
public:
    QString name() const override { return "MODISQC"; }
    QString description() const override {
        return "Extract and decode MODIS quality control bit fields.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input MODIS QC raster"),
            ParameterDef::output("output", "Output decoded QC raster"),
            ParameterDef::integer("start_bit", "Start bit position", 0, 0, 31),
            ParameterDef::integer("num_bits", "Number of bits to extract", 2, 1, 16),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read MODIS QC input"); return false; }

        int startBit = parameter("start_bit").toInt();
        int numBits = parameter("num_bits").toInt();
        int mask = (1 << numBits) - 1;

        int64_t total = raster->cellCount();
        for (int b = 0; b < raster->bands(); ++b) {
            auto& data = raster->data(b);
            for (int64_t i = 0; i < total; ++i) {
                int val = static_cast<int>(data[i]);
                data[i] = static_cast<double>((val >> startBit) & mask);
            }
            reportProgress(static_cast<double>(b + 1) / raster->bands());
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(ModisQcModule)

} // namespace aplaceholder
