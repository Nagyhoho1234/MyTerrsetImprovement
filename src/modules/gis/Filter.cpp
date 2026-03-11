#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>

namespace aplaceholder {

class FilterModule : public Module {
public:
    QString name() const override { return "FILTER"; }
    QString description() const override {
        return "Digital spatial filters for raster images. "
               "Applies neighborhood-based convolution filters including "
               "smoothing, statistical, and edge detection operations.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::output("output", "Output filtered image"),
            ParameterDef::combo("filter_type", "Filter type",
                {"Mean", "Median", "Mode", "Min", "Max", "Laplacian"}, 0,
                "Type of spatial filter to apply"),
            ParameterDef::integer("kernel_size", "Kernel size", 3, 3, 25,
                "Size of the filter kernel (must be odd)"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        int filterType = parameter("filter_type").toInt();
        int kernelSize = parameter("kernel_size").toInt();

        // Enforce odd kernel size
        if (kernelSize % 2 == 0) {
            kernelSize++;
        }
        int half = kernelSize / 2;

        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        // Build Laplacian kernel if needed
        // Standard Laplacian: center = -4 (or scaled for larger kernels),
        // neighbors = 1. For arbitrary kernel sizes, use a discrete Laplacian
        // approximation: center = -(kernelSize*kernelSize - 1), all others = 1
        // But standard 3x3 Laplacian is: [0 1 0; 1 -4 1; 0 1 0]
        // For larger kernels, we use: all neighbors = 1, center = -(count of neighbors)

        reportProgress(0.0, "Applying filter...");
        std::vector<double> neighborhood;
        neighborhood.reserve(static_cast<size_t>(kernelSize) * kernelSize);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;

                if (hasND && data[idx] == noData) {
                    out[idx] = noData;
                    continue;
                }

                if (filterType == 5) {
                    // Laplacian edge detection
                    // Use discrete Laplacian: sum of neighbors minus center * count
                    double sum = 0.0;
                    int count = 0;
                    double center = data[idx];

                    for (int kr = -half; kr <= half; ++kr) {
                        for (int kc = -half; kc <= half; ++kc) {
                            if (kr == 0 && kc == 0) continue;
                            int nr = r + kr;
                            int nc = c + kc;
                            if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
                                continue;
                            size_t nIdx = static_cast<size_t>(nr) * cols + nc;
                            if (hasND && data[nIdx] == noData)
                                continue;
                            sum += data[nIdx];
                            count++;
                        }
                    }

                    // Laplacian = sum_of_neighbors - count * center
                    out[idx] = (count > 0) ? (sum - count * center) : 0.0;
                } else {
                    // Collect valid neighborhood values
                    neighborhood.clear();
                    for (int kr = -half; kr <= half; ++kr) {
                        for (int kc = -half; kc <= half; ++kc) {
                            int nr = r + kr;
                            int nc = c + kc;
                            if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
                                continue;
                            size_t nIdx = static_cast<size_t>(nr) * cols + nc;
                            if (hasND && data[nIdx] == noData)
                                continue;
                            neighborhood.push_back(data[nIdx]);
                        }
                    }

                    if (neighborhood.empty()) {
                        out[idx] = noData;
                        continue;
                    }

                    switch (filterType) {
                    case 0: { // Mean
                        double sum = 0.0;
                        for (double v : neighborhood) sum += v;
                        out[idx] = sum / neighborhood.size();
                        break;
                    }
                    case 1: { // Median
                        std::sort(neighborhood.begin(), neighborhood.end());
                        size_t n = neighborhood.size();
                        if (n % 2 == 1)
                            out[idx] = neighborhood[n / 2];
                        else
                            out[idx] = (neighborhood[n / 2 - 1] + neighborhood[n / 2]) / 2.0;
                        break;
                    }
                    case 2: { // Mode
                        // Find the most frequent value (integer-based)
                        std::map<int64_t, int> freq;
                        for (double v : neighborhood) {
                            int64_t key = static_cast<int64_t>(std::round(v));
                            freq[key]++;
                        }
                        int maxFreq = 0;
                        int64_t modeVal = 0;
                        for (const auto& kv : freq) {
                            if (kv.second > maxFreq) {
                                maxFreq = kv.second;
                                modeVal = kv.first;
                            }
                        }
                        out[idx] = static_cast<double>(modeVal);
                        break;
                    }
                    case 3: { // Min
                        out[idx] = *std::min_element(neighborhood.begin(), neighborhood.end());
                        break;
                    }
                    case 4: { // Max
                        out[idx] = *std::max_element(neighborhood.begin(), neighborhood.end());
                        break;
                    }
                    }
                }
            }
            if (r % 100 == 0)
                reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(FilterModule)

} // namespace aplaceholder
