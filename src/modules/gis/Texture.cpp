#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <map>

namespace aplaceholder {

class TextureModule : public Module {
public:
    QString name() const override { return "TEXTURE"; }
    QString description() const override {
        return "Compute texture metrics from raster images using a moving window. "
               "Supports variance, coefficient of variation, range, contrast, "
               "and angular second moment (homogeneity) texture measures.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::output("output", "Output texture metric image"),
            ParameterDef::combo("metric", "Texture metric",
                {"Variance", "Coefficient of Variation", "Range",
                 "Contrast", "Angular Second Moment"}, 0,
                "Texture measure to compute"),
            ParameterDef::integer("kernel_size", "Window size", 3, 3, 25,
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

        reportProgress(0.0, "Computing texture...");

        std::vector<double> neighborhood;
        neighborhood.reserve(static_cast<size_t>(kernelSize) * kernelSize);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;
                if (hasND && data[idx] == noData) {
                    out[idx] = noData;
                    continue;
                }

                neighborhood.clear();
                for (int kr = -half; kr <= half; ++kr) {
                    for (int kc = -half; kc <= half; ++kc) {
                        int nr = r + kr, nc = c + kc;
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                        size_t nIdx = static_cast<size_t>(nr) * cols + nc;
                        if (hasND && data[nIdx] == noData) continue;
                        neighborhood.push_back(data[nIdx]);
                    }
                }

                if (neighborhood.empty()) {
                    out[idx] = noData;
                    continue;
                }

                size_t n = neighborhood.size();

                switch (metric) {
                case 0: { // Variance
                    double sum = 0, sumSq = 0;
                    for (double v : neighborhood) { sum += v; sumSq += v * v; }
                    double mean = sum / n;
                    out[idx] = sumSq / n - mean * mean;
                    break;
                }
                case 1: { // Coefficient of variation
                    double sum = 0, sumSq = 0;
                    for (double v : neighborhood) { sum += v; sumSq += v * v; }
                    double mean = sum / n;
                    double var = sumSq / n - mean * mean;
                    out[idx] = (std::abs(mean) > 1e-10) ? std::sqrt(std::max(var, 0.0)) / std::abs(mean) : 0.0;
                    break;
                }
                case 2: { // Range
                    double minV = neighborhood[0], maxV = neighborhood[0];
                    for (double v : neighborhood) {
                        if (v < minV) minV = v;
                        if (v > maxV) maxV = v;
                    }
                    out[idx] = maxV - minV;
                    break;
                }
                case 3: { // Contrast (sum of squared differences between neighbors)
                    double contrast = 0.0;
                    for (size_t i = 0; i < n; ++i) {
                        for (size_t j = i + 1; j < n; ++j) {
                            double diff = neighborhood[i] - neighborhood[j];
                            contrast += diff * diff;
                        }
                    }
                    int pairs = static_cast<int>(n * (n - 1) / 2);
                    out[idx] = (pairs > 0) ? contrast / pairs : 0.0;
                    break;
                }
                case 4: { // Angular Second Moment (homogeneity)
                    std::map<int, int> freq;
                    for (double v : neighborhood) {
                        freq[static_cast<int>(std::round(v))]++;
                    }
                    double asm_val = 0.0;
                    for (const auto& kv : freq) {
                        double p = static_cast<double>(kv.second) / n;
                        asm_val += p * p;
                    }
                    out[idx] = asm_val;
                    break;
                }
                }
            }
            if (r % 100 == 0) reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(TextureModule)

} // namespace aplaceholder
