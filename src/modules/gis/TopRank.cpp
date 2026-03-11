#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class TopRankModule : public Module {
public:
    QString name() const override { return "TOPRANK"; }
    QString description() const override {
        return "Selects the top-ranked cells from a suitability surface. "
               "Outputs a binary raster with 1 for the best cells, used as "
               "a simple choice heuristic after MCE evaluation.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input suitability image"),
            ParameterDef::output("output", "Output selection image"),
            ParameterDef::combo("mode", "Selection mode",
                {"Number of cells", "Percentage"}, 0,
                "Select top N cells or top percentage"),
            ParameterDef::integer("num_cells", "Number of cells to select",
                100, 1, 999999999,
                "Number of top-ranked cells to extract (used in Number mode)"),
            ParameterDef::real("percentage", "Percentage to select",
                5.0, 0.0001, 100.0,
                "Percentage of valid cells to extract (used in Percentage mode)"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = input->cols(), rows = input->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();

        // Collect valid cell indices
        reportProgress(0.0, "Collecting valid cells...");
        std::vector<int64_t> validIndices;
        validIndices.reserve(total);
        for (int64_t i = 0; i < total; ++i) {
            if (!hasND || data[i] != noData)
                validIndices.push_back(i);
        }

        if (validIndices.empty()) {
            setError("No valid cells found in input");
            return false;
        }

        // Determine how many cells to select
        int mode = parameter("mode").toInt();
        int64_t numSelect;
        if (mode == 0) {
            numSelect = parameter("num_cells").toLongLong();
        } else {
            double pct = parameter("percentage").toDouble();
            numSelect = static_cast<int64_t>(
                std::ceil(static_cast<double>(validIndices.size()) * pct / 100.0));
        }

        if (numSelect <= 0) numSelect = 1;
        if (numSelect > static_cast<int64_t>(validIndices.size()))
            numSelect = static_cast<int64_t>(validIndices.size());

        // Sort valid indices by suitability value (descending = highest first)
        reportProgress(0.2, "Ranking cells...");
        std::sort(validIndices.begin(), validIndices.end(),
            [&data](int64_t a, int64_t b) {
                return data[a] > data[b];
            });

        // Create output: 1 for top cells, 0 for others
        reportProgress(0.7, "Creating output...");
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(input->noDataValue());
        auto& out = output.data(0);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData)
                out[i] = noData;
            else
                out[i] = 0.0;
        }

        for (int64_t i = 0; i < numSelect; ++i)
            out[validIndices[i]] = 1.0;

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(TopRankModule)

} // namespace aplaceholder
