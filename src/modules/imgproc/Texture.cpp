#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class TextureModule : public Module {
public:
    QString name() const override { return "TEXTURE"; }
    QString description() const override {
        return "GLCM texture measures. Computes Grey-Level Co-occurrence Matrix and derives "
               "texture statistics: contrast, dissimilarity, homogeneity, energy, entropy, correlation.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output texture image"),
            ParameterDef::combo("measure", "Texture measure",
                {"Contrast", "Dissimilarity", "Homogeneity", "Energy", "Entropy", "Correlation"},
                0, "GLCM-derived texture measure to compute"),
            ParameterDef::integer("window_size", "Window size (pixels)", 7, 3, 101,
                "Size of the moving window (must be odd)"),
            ParameterDef::combo("direction", "Direction (degrees)",
                {"0", "45", "90", "135"}, 0,
                "Direction of pixel pair offset for GLCM computation"),
            ParameterDef::integer("num_levels", "Number of grey levels", 16, 2, 256,
                "Quantization levels for GLCM computation"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input image");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        int measureIdx = parameter("measure").toInt();
        int winSize = parameter("window_size").toInt();
        int dirIdx = parameter("direction").toInt();
        int numLevels = parameter("num_levels").toInt();

        // Ensure window size is odd
        if (winSize % 2 == 0) winSize++;
        int halfWin = winSize / 2;

        // Direction offsets (dx, dy) for 0, 45, 90, 135 degrees
        int dx = 0, dy = 0;
        switch (dirIdx) {
        case 0: dx = 1; dy = 0; break;   // 0 degrees (horizontal right)
        case 1: dx = 1; dy = -1; break;  // 45 degrees (upper-right diagonal)
        case 2: dx = 0; dy = -1; break;  // 90 degrees (vertical up)
        case 3: dx = -1; dy = -1; break; // 135 degrees (upper-left diagonal)
        }

        reportProgress(0.05, "Quantizing input...");

        // ------------------------------------------------------------------
        // 1. Quantize input to [0, numLevels-1]
        // ------------------------------------------------------------------
        const auto& src = raster->data(0);
        auto stats = raster->computeStats(0);
        double srcMin = stats.min;
        double srcMax = stats.max;
        double srcRange = srcMax - srcMin;

        std::vector<int> quantized(total);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && src[i] == noData) {
                quantized[i] = -1; // mark as NoData
            } else {
                double norm = (srcRange > 0) ? (src[i] - srcMin) / srcRange : 0.0;
                int q = static_cast<int>(norm * (numLevels - 1) + 0.5);
                quantized[i] = std::clamp(q, 0, numLevels - 1);
            }
        }

        reportProgress(0.1, "Computing GLCM texture...");

        // ------------------------------------------------------------------
        // 2. Create output
        // ------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(noData);

        auto& out = output.data(0);

        // ------------------------------------------------------------------
        // 3. For each pixel, compute GLCM within the window and derive measure
        // ------------------------------------------------------------------
        std::vector<double> glcm(numLevels * numLevels);

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                int64_t idx = static_cast<int64_t>(row) * cols + col;

                if (hasND && quantized[idx] < 0) {
                    out[idx] = noData;
                    continue;
                }

                // Build GLCM for this window
                std::fill(glcm.begin(), glcm.end(), 0.0);
                int pairCount = 0;

                for (int wy = -halfWin; wy <= halfWin; ++wy) {
                    for (int wx = -halfWin; wx <= halfWin; ++wx) {
                        int r1 = row + wy;
                        int c1 = col + wx;
                        int r2 = r1 + dy;
                        int c2 = c1 + dx;

                        if (r1 < 0 || r1 >= rows || c1 < 0 || c1 >= cols) continue;
                        if (r2 < 0 || r2 >= rows || c2 < 0 || c2 >= cols) continue;

                        int q1 = quantized[static_cast<int64_t>(r1) * cols + c1];
                        int q2 = quantized[static_cast<int64_t>(r2) * cols + c2];

                        if (q1 < 0 || q2 < 0) continue;

                        // Symmetric GLCM
                        glcm[q1 * numLevels + q2] += 1.0;
                        glcm[q2 * numLevels + q1] += 1.0;
                        pairCount += 2;
                    }
                }

                if (pairCount == 0) {
                    out[idx] = 0.0;
                    continue;
                }

                // Normalize GLCM
                double invCount = 1.0 / pairCount;
                for (int k = 0; k < numLevels * numLevels; ++k)
                    glcm[k] *= invCount;

                // Compute the requested measure
                double result = 0.0;

                if (measureIdx == 5) {
                    // Correlation needs mean and stddev of marginals
                    double muI = 0.0, muJ = 0.0;
                    double sigI = 0.0, sigJ = 0.0;

                    // Marginal means
                    for (int i = 0; i < numLevels; ++i) {
                        double pI = 0.0, pJ = 0.0;
                        for (int j = 0; j < numLevels; ++j) {
                            pI += glcm[i * numLevels + j];
                            pJ += glcm[j * numLevels + i];
                        }
                        muI += i * pI;
                        muJ += i * pJ;
                    }

                    // Marginal variances
                    for (int i = 0; i < numLevels; ++i) {
                        double pI = 0.0, pJ = 0.0;
                        for (int j = 0; j < numLevels; ++j) {
                            pI += glcm[i * numLevels + j];
                            pJ += glcm[j * numLevels + i];
                        }
                        sigI += (i - muI) * (i - muI) * pI;
                        sigJ += (i - muJ) * (i - muJ) * pJ;
                    }

                    sigI = std::sqrt(sigI);
                    sigJ = std::sqrt(sigJ);

                    if (sigI > 1e-12 && sigJ > 1e-12) {
                        for (int i = 0; i < numLevels; ++i)
                            for (int j = 0; j < numLevels; ++j)
                                result += (i - muI) * (j - muJ) * glcm[i * numLevels + j];
                        result /= (sigI * sigJ);
                    }
                } else {
                    for (int i = 0; i < numLevels; ++i) {
                        for (int j = 0; j < numLevels; ++j) {
                            double p = glcm[i * numLevels + j];
                            if (p <= 0.0) continue;

                            switch (measureIdx) {
                            case 0: // Contrast
                                result += (i - j) * (i - j) * p;
                                break;
                            case 1: // Dissimilarity
                                result += std::abs(i - j) * p;
                                break;
                            case 2: // Homogeneity (Inverse Difference Moment)
                                result += p / (1.0 + (i - j) * (i - j));
                                break;
                            case 3: // Energy (Angular Second Moment)
                                result += p * p;
                                break;
                            case 4: // Entropy
                                result -= p * std::log(p);
                                break;
                            }
                        }
                    }
                }

                out[idx] = result;
            }

            if (row % 50 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(row) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(TextureModule)

} // namespace aplaceholder
