#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <map>
#include <vector>

namespace aplaceholder {

class PatternModule : public Module {
public:
    QString name() const override { return "PATTERN"; }
    QString description() const override {
        return "Compute spatial pattern metrics over a moving neighborhood. "
               "Calculates entropy, dominance, and diversity indices from "
               "categorical raster data using a user-defined kernel.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input categorical raster image"),
            ParameterDef::output("output", "Output pattern metric image"),
            ParameterDef::combo("metric", "Pattern metric",
                {"Entropy", "Dominance", "Diversity"}, 0,
                "Entropy: Shannon entropy; Dominance: 1 - evenness; Diversity: number of classes"),
            ParameterDef::integer("kernel_size", "Kernel size", 3, 3, 25,
                "Size of the moving window (must be odd)"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        int metric = parameter("metric").toInt();
        int kernelSize = parameter("kernel_size").toInt();
        if (kernelSize % 2 == 0) kernelSize++;
        int half = kernelSize / 2;

        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing pattern metrics...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;

                if (hasND && data[idx] == noData) {
                    out[idx] = noData;
                    continue;
                }

                // Collect category frequencies in the neighborhood
                std::map<int, int> freq;
                int total = 0;
                for (int kr = -half; kr <= half; ++kr) {
                    for (int kc = -half; kc <= half; ++kc) {
                        int nr = r + kr, nc = c + kc;
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                        size_t nIdx = static_cast<size_t>(nr) * cols + nc;
                        if (hasND && data[nIdx] == noData) continue;
                        int cat = static_cast<int>(std::round(data[nIdx]));
                        freq[cat]++;
                        total++;
                    }
                }

                if (total == 0) {
                    out[idx] = noData;
                    continue;
                }

                if (metric == 2) {
                    // Diversity: number of distinct classes
                    out[idx] = static_cast<double>(freq.size());
                } else {
                    // Shannon entropy: H = -sum(p_i * ln(p_i))
                    double entropy = 0.0;
                    for (const auto& kv : freq) {
                        double p = static_cast<double>(kv.second) / total;
                        if (p > 0.0) entropy -= p * std::log(p);
                    }

                    if (metric == 0) {
                        // Entropy
                        out[idx] = entropy;
                    } else {
                        // Dominance: ln(numClasses) - entropy
                        double maxEntropy = std::log(static_cast<double>(freq.size()));
                        out[idx] = (maxEntropy > 0.0) ? (maxEntropy - entropy) : 0.0;
                    }
                }
            }
            if (r % 100 == 0) reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PatternModule)

} // namespace aplaceholder
