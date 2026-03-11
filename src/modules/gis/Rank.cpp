#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class RankModule : public Module {
public:
    QString name() const override { return "RANK"; }
    QString description() const override {
        return "Rank-orders pixel values in a suitability surface. Output pixel "
               "values are rank positions. Divide by maximum rank to produce a "
               "relative risk image (0-1).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input suitability image"),
            ParameterDef::output("output", "Output ranked image"),
            ParameterDef::boolean("ascending", "Ascending order",
                false,
                "If false (default), highest suitability gets rank 1. "
                "If true, lowest value gets rank 1."),
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
        bool ascending = parameter("ascending").toBool();

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

        // Sort valid indices by value
        reportProgress(0.1, "Sorting cells...");
        if (ascending) {
            // Ascending: lowest value gets rank 1
            std::sort(validIndices.begin(), validIndices.end(),
                [&data](int64_t a, int64_t b) {
                    return data[a] < data[b];
                });
        } else {
            // Descending (default): highest value gets rank 1
            std::sort(validIndices.begin(), validIndices.end(),
                [&data](int64_t a, int64_t b) {
                    return data[a] > data[b];
                });
        }

        // Assign ranks
        reportProgress(0.6, "Assigning ranks...");
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(input->noDataValue());
        auto& out = output.data(0);

        // Initialize all to nodata
        for (int64_t i = 0; i < total; ++i)
            out[i] = noData;

        // Assign ranks with tied values getting the same rank
        int64_t rank = 1;
        for (int64_t i = 0; i < static_cast<int64_t>(validIndices.size()); ++i) {
            if (i > 0 && data[validIndices[i]] != data[validIndices[i - 1]])
                rank = i + 1;
            out[validIndices[i]] = static_cast<double>(rank);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(RankModule)

} // namespace aplaceholder
