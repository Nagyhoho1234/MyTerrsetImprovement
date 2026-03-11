#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class XyzMonthlyModule : public Module {
public:
    QString name() const override { return "XYZMONTHLY"; }
    QString description() const override {
        return "Import XYZ monthly time-series data to TerrSet raster format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input XYZ monthly data file"),
            ParameterDef::output("output", "Output TerrSet raster (base name)"),
            ParameterDef::integer("start_year", "Start year", 2000, 1900, 2100),
            ParameterDef::integer("start_month", "Start month", 1, 1, 12),
            ParameterDef::integer("num_months", "Number of months", 12, 1, 9999),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input file"); return false; }

        reportProgress(0.5, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString(), "RST");
    }
};

REGISTER_MODULE(XyzMonthlyModule)

} // namespace aplaceholder
