#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <set>
#include <string>

namespace aplaceholder {

class BreakoutModule : public Module {
public:
    QString name() const override { return "BREAKOUT"; }
    QString description() const override {
        return "Break out individual categories from a categorical raster into "
               "separate Boolean images. For each unique value in the input, "
               "creates a Boolean raster (1 where value matches, 0 elsewhere).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input categorical raster image"),
            ParameterDef::output("output_prefix", "Output filename prefix"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Collect unique values
        reportProgress(0.0, "Scanning for unique categories...");
        std::set<int> uniqueValues;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            uniqueValues.insert(static_cast<int>(data[i]));
        }

        if (uniqueValues.empty()) {
            setError("No valid pixels found in input raster");
            return false;
        }

        QString prefix = parameter("output_prefix").toString();
        int catIndex = 0;
        int numCats = static_cast<int>(uniqueValues.size());

        for (int val : uniqueValues) {
            reportProgress(static_cast<double>(catIndex) / numCats,
                "Creating Boolean image for category " + QString::number(val) + "...");

            Raster output(cols, rows, 1, DataType::Float64);
            output.setGeoTransform(input->geoTransform());
            output.setProjection(input->projection());
            output.setNoDataValue(noData);
            auto& out = output.data(0);

            for (int64_t i = 0; i < total; ++i) {
                if (hasND && data[i] == noData) {
                    out[i] = noData;
                    continue;
                }
                out[i] = (static_cast<int>(data[i]) == val) ? 1.0 : 0.0;
            }

            QString outPath = prefix + "_" + QString::number(val);
            if (!GdalIO::write(output, outPath)) {
                setError("Failed to write output for category " + QString::number(val));
                return false;
            }

            ++catIndex;
        }

        reportProgress(1.0, "Complete.");
        return true;
    }
};

REGISTER_MODULE(BreakoutModule)

} // namespace aplaceholder
