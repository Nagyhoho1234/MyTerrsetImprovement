#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <map>
#include <string>

namespace aplaceholder {

class CellAtomModule : public Module {
public:
    QString name() const override { return "CELLATOM"; }
    QString description() const override {
        return "Cellular automaton simulation. Each iteration examines a 3x3 "
               "neighborhood majority for every cell and applies transition "
               "rules to produce the next state. Outputs the final state raster.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("initial_state", "Input initial state raster"),
            ParameterDef::file("transition_rules", "Transition rules raster"),
            ParameterDef::integer("iterations", "Number of iterations", 10, 1, 10000,
                "Number of simulation iterations to run"),
            ParameterDef::output("output", "Output final state raster"),
        };
    }

    bool execute() override {
        auto initialInput = GdalIO::read(parameter("initial_state").toString());
        if (!initialInput) { setError("Failed to read initial state raster"); return false; }

        auto rulesInput = GdalIO::read(parameter("transition_rules").toString());
        if (!rulesInput) { setError("Failed to read transition rules raster"); return false; }

        int cols = initialInput->cols(), rows = initialInput->rows();

        if (rulesInput->cols() != cols || rulesInput->rows() != rows) {
            setError("Initial state and transition rules rasters must have the same dimensions");
            return false;
        }

        double noData = initialInput->noDataValue();
        bool hasND = initialInput->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int iterations = parameter("iterations").toInt();

        const auto& rulesData = rulesInput->data(0);

        // Working buffers: current state and next state
        std::vector<double> current(total);
        std::vector<double> next(total);

        // Copy initial state
        const auto& initData = initialInput->data(0);
        for (int64_t i = 0; i < total; ++i) {
            current[i] = initData[i];
        }

        for (int iter = 0; iter < iterations; ++iter) {
            reportProgress(static_cast<double>(iter) / iterations,
                "Iteration " + QString::number(iter + 1) + " of " + QString::number(iterations) + "...");

            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    int64_t idx = static_cast<int64_t>(r) * cols + c;

                    if (hasND && current[idx] == noData) {
                        next[idx] = noData;
                        continue;
                    }

                    // Count neighborhood majority in 3x3 window
                    std::map<int, int> counts;
                    for (int dr = -1; dr <= 1; ++dr) {
                        for (int dc = -1; dc <= 1; ++dc) {
                            int nr = r + dr;
                            int nc = c + dc;
                            if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                            int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                            if (hasND && current[nIdx] == noData) continue;
                            int val = static_cast<int>(current[nIdx]);
                            counts[val]++;
                        }
                    }

                    // Find majority value
                    int majorityVal = static_cast<int>(current[idx]);
                    int maxCount = 0;
                    for (const auto& pair : counts) {
                        if (pair.second > maxCount) {
                            maxCount = pair.second;
                            majorityVal = pair.first;
                        }
                    }

                    // Apply transition: use the transition rules raster value
                    // as a modifier. If the rules cell is non-zero and majority
                    // differs from current, transition to majority; otherwise keep.
                    double ruleVal = rulesData[idx];
                    if (hasND && ruleVal == noData) {
                        next[idx] = current[idx];
                    } else if (static_cast<int>(ruleVal) != 0 &&
                               majorityVal != static_cast<int>(current[idx])) {
                        next[idx] = static_cast<double>(majorityVal);
                    } else {
                        next[idx] = current[idx];
                    }
                }
            }

            std::swap(current, next);
        }

        // Write final state
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(initialInput->geoTransform());
        output.setProjection(initialInput->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            out[i] = current[i];
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(CellAtomModule)

} // namespace aplaceholder
