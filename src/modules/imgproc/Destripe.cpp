#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class DestripeModule : public Module {
public:
    QString name() const override { return "DESTRIPE"; }
    QString description() const override {
        return "Remove striping artifacts from satellite imagery. Detects and corrects "
               "line-by-line radiometric variations by normalizing each scan line to "
               "match the overall scene mean and standard deviation.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image (with striping)"),
            ParameterDef::output("output", "Output destriped image"),
            ParameterDef::combo("direction", "Stripe direction",
                                {"horizontal", "vertical"}, 0,
                                "Horizontal for scanner-type sensors (e.g., Landsat MSS/TM), "
                                "vertical for pushbroom sensors (e.g., SPOT)."),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString outputPath = parameter("output").toString();
        int dirIdx = parameter("direction").toInt();
        bool horizontal = (dirIdx == 0);

        auto input = GdalIO::read(inputPath);
        if (!input) {
            setError("Failed to read input image: " + inputPath);
            return false;
        }

        int cols = input->cols();
        int rows = input->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        bool hasND = input->hasNoData();
        double inputND = input->noDataValue();
        const auto& inData = input->data(0);

        // Step 1: Compute overall scene mean and standard deviation
        double globalSum = 0.0;
        double globalSumSq = 0.0;
        int64_t globalCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            double val = inData[i];
            if (hasND && val == inputND)
                continue;
            globalSum += val;
            globalSumSq += val * val;
            ++globalCount;
        }

        if (globalCount == 0) {
            setError("Image contains no valid pixels");
            return false;
        }

        double globalMean = globalSum / globalCount;
        double globalVariance = (globalSumSq / globalCount) - (globalMean * globalMean);
        double globalStddev = std::sqrt(std::max(globalVariance, 0.0));

        // Step 2: Compute per-line mean and standard deviation
        // For horizontal stripes: each row is a line
        // For vertical stripes: each column is a line
        int numLines = horizontal ? rows : cols;
        int lineLength = horizontal ? cols : rows;

        std::vector<double> lineMean(numLines, 0.0);
        std::vector<double> lineStddev(numLines, 1.0);
        std::vector<int64_t> lineCount(numLines, 0);

        for (int line = 0; line < numLines; ++line) {
            double sum = 0.0;
            double sumSq = 0.0;
            int64_t count = 0;

            for (int pos = 0; pos < lineLength; ++pos) {
                int col = horizontal ? pos : line;
                int row = horizontal ? line : pos;
                int64_t idx = static_cast<int64_t>(row) * cols + col;
                double val = inData[idx];
                if (hasND && val == inputND)
                    continue;
                sum += val;
                sumSq += val * val;
                ++count;
            }

            if (count > 0) {
                lineMean[line] = sum / count;
                double var = (sumSq / count) - (lineMean[line] * lineMean[line]);
                lineStddev[line] = std::sqrt(std::max(var, 0.0));
                lineCount[line] = count;
            }
        }

        reportProgress(0.4, "Computing correction...");

        // Step 3: Rescale each line to match global mean/stddev
        // corrected = globalMean + globalStddev * (original - lineMean) / lineStddev
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        double noData = -9999.0;
        output.setNoDataValue(noData);
        auto& outData = output.data(0);

        for (int line = 0; line < numLines; ++line) {
            for (int pos = 0; pos < lineLength; ++pos) {
                int col = horizontal ? pos : line;
                int row = horizontal ? line : pos;
                int64_t idx = static_cast<int64_t>(row) * cols + col;
                double val = inData[idx];

                if (hasND && val == inputND) {
                    outData[idx] = noData;
                    continue;
                }

                if (lineCount[line] == 0 || lineStddev[line] < 1e-10) {
                    // If the line has no variance, just shift to global mean
                    outData[idx] = globalMean;
                } else {
                    outData[idx] = globalMean +
                                   globalStddev * (val - lineMean[line]) / lineStddev[line];
                }
            }

            if (line % 100 == 0)
                reportProgress(0.4 + 0.6 * static_cast<double>(line) / numLines);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
    }
};

REGISTER_MODULE(DestripeModule)

} // namespace aplaceholder
