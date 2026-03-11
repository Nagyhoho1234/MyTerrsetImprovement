#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class NanFixModule : public Module {
public:
    QString name() const override { return "NANFIX"; }
    QString description() const override {
        return "Replace NaN and Inf values in a raster with a specified value.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::real("replace_value", "Replacement value", 0.0,
                -1e38, 1e38, "Value to replace NaN/Inf with"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input raster"); return false; }

        double replaceVal = parameter("replace_value").toDouble();
        int bands = raster->bands();
        int64_t total = raster->cellCount();

        for (int b = 0; b < bands; ++b) {
            auto& data = raster->data(b);
            for (int64_t i = 0; i < total; ++i) {
                if (std::isnan(data[i]) || std::isinf(data[i]))
                    data[i] = replaceVal;
            }
            reportProgress(static_cast<double>(b + 1) / bands);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(NanFixModule)

} // namespace aplaceholder
