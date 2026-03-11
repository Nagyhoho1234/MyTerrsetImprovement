#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class DenoiseModule : public Module {
public:
    QString name() const override { return "DENOISE"; }
    QString description() const override {
        return "Image denoising using the Lee filter for speckle reduction in SAR imagery. "
               "Applies an adaptive filter: pixel = mean + K*(pixel - mean) where "
               "K = variance_local / (variance_local + variance_noise).";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output denoised image"),
            ParameterDef::integer("window_size", "Filter window size (odd)", 5, 3, 99,
                "Size of the moving window (must be odd)"),
            ParameterDef::real("noise_variance", "Noise variance (0 = auto-estimate)", 0.0,
                0.0, 999999.0,
                "Noise variance. If 0, it is auto-estimated from the image."),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString outputPath = parameter("output").toString();
        int winSize = parameter("window_size").toInt();
        double noiseVar = parameter("noise_variance").toDouble();

        if (winSize % 2 == 0) {
            setError("Window size must be odd.");
            return false;
        }

        auto input = GdalIO::read(inputPath);
        if (!input) {
            setError("Failed to read input image: " + inputPath);
            return false;
        }

        int cols = input->cols();
        int rows = input->rows();
        int numBands = input->bands();
        int half = winSize / 2;

        bool hasND = input->hasNoData();
        double inputND = input->noDataValue();
        double outNoData = -9999.0;

        Raster output(cols, rows, numBands, DataType::Float32);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(outNoData);

        for (int b = 0; b < numBands; ++b) {
            reportProgress(static_cast<double>(b) / numBands,
                           QString("Processing band %1...").arg(b + 1));

            const auto& inData = input->data(b);
            auto& outData = output.data(b);

            // Auto-estimate noise variance if needed
            double sigmaN = noiseVar;
            if (sigmaN <= 0.0) {
                // Estimate from the whole image using the mean of local variances
                // in small windows at the image edges
                double sumVar = 0.0;
                int count = 0;
                for (int r = half; r < rows - half; r += winSize) {
                    for (int c = half; c < cols - half; c += winSize) {
                        double sum = 0.0, sumSq = 0.0;
                        int n = 0;
                        for (int wr = -half; wr <= half; ++wr) {
                            for (int wc = -half; wc <= half; ++wc) {
                                double v = inData[static_cast<int64_t>(r + wr) * cols + (c + wc)];
                                if (hasND && v == inputND) continue;
                                sum += v;
                                sumSq += v * v;
                                ++n;
                            }
                        }
                        if (n > 1) {
                            double mean = sum / n;
                            double var = (sumSq - n * mean * mean) / (n - 1);
                            sumVar += var;
                            ++count;
                        }
                    }
                }
                sigmaN = (count > 0) ? sumVar / count : 1.0;
            }

            // Apply Lee filter
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    int64_t idx = static_cast<int64_t>(r) * cols + c;
                    double center = inData[idx];

                    if (hasND && center == inputND) {
                        outData[idx] = outNoData;
                        continue;
                    }

                    // Compute local statistics
                    double sum = 0.0, sumSq = 0.0;
                    int n = 0;

                    int rMin = std::max(0, r - half);
                    int rMax = std::min(rows - 1, r + half);
                    int cMin = std::max(0, c - half);
                    int cMax = std::min(cols - 1, c + half);

                    for (int wr = rMin; wr <= rMax; ++wr) {
                        for (int wc = cMin; wc <= cMax; ++wc) {
                            double v = inData[static_cast<int64_t>(wr) * cols + wc];
                            if (hasND && v == inputND) continue;
                            sum += v;
                            sumSq += v * v;
                            ++n;
                        }
                    }

                    if (n < 2) {
                        outData[idx] = center;
                        continue;
                    }

                    double localMean = sum / n;
                    double localVar = (sumSq - n * localMean * localMean) / (n - 1);
                    if (localVar < 0.0) localVar = 0.0;

                    // Lee filter coefficient
                    double K = localVar / (localVar + sigmaN);
                    outData[idx] = localMean + K * (center - localMean);
                }

                if (r % 500 == 0)
                    reportProgress((static_cast<double>(b) + static_cast<double>(r) / rows) / numBands);
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
    }
};

REGISTER_MODULE(DenoiseModule)

} // namespace aplaceholder
