#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>

namespace aplaceholder {

class PurifyModule : public Module {
public:
    QString name() const override { return "PURIFY"; }
    QString description() const override {
        return "Purify training signatures. Removes outlier pixels from training site "
               "signatures by iteratively discarding pixels beyond a specified number "
               "of standard deviations from the class mean in spectral space.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input multi-band image"),
            ParameterDef::file("training", "Training site image (class labels)",
                "Integer raster where each value identifies a training class (0 = background)"),
            ParameterDef::output("output", "Output purified training site image"),
            ParameterDef::real("std_threshold", "Standard deviation threshold", 2.0, 0.5, 10.0,
                "Pixels beyond this many std devs from class mean are removed"),
            ParameterDef::integer("max_iterations", "Maximum iterations", 5, 1, 50,
                "Maximum number of purification passes"),
            ParameterDef::real("min_retain", "Minimum retention fraction", 0.5, 0.1, 1.0,
                "Stop if fewer than this fraction of original pixels remain"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        auto training = GdalIO::read(parameter("training").toString());
        if (!raster || !training) {
            setError("Failed to read input image or training sites");
            return false;
        }

        if (raster->cols() != training->cols() || raster->rows() != training->rows()) {
            setError("Input image and training image must have the same dimensions");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int nBands = raster->bands();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double stdThresh = parameter("std_threshold").toDouble();
        int maxIter = parameter("max_iterations").toInt();
        double minRetain = parameter("min_retain").toDouble();

        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        std::vector<const std::vector<double>*> bands(nBands);
        for (int b = 0; b < nBands; ++b)
            bands[b] = &raster->data(b);
        const auto& trainData = training->data(0);

        reportProgress(0.05, "Identifying training classes...");

        // Identify unique classes (non-zero values)
        std::vector<int> classIds;
        {
            std::vector<bool> seen(10001, false);
            for (int64_t i = 0; i < total; ++i) {
                int cls = static_cast<int>(trainData[i]);
                if (cls > 0 && cls <= 10000 && !seen[cls]) {
                    seen[cls] = true;
                    classIds.push_back(cls);
                }
            }
            std::sort(classIds.begin(), classIds.end());
        }

        if (classIds.empty()) {
            setError("No training classes found (all values are 0 or nodata)");
            return false;
        }

        // Mutable copy of training labels
        std::vector<double> labels(trainData.begin(), trainData.end());

        // Track original counts per class
        std::vector<int64_t> origCounts(classIds.size(), 0);
        for (size_t c = 0; c < classIds.size(); ++c)
            for (int64_t i = 0; i < total; ++i)
                if (static_cast<int>(labels[i]) == classIds[c])
                    origCounts[c]++;

        reportProgress(0.10, "Purifying training signatures...");

        for (int iter = 0; iter < maxIter; ++iter) {
            bool anyRemoved = false;

            for (size_t ci = 0; ci < classIds.size(); ++ci) {
                int cls = classIds[ci];

                // Compute class mean and stddev per band
                std::vector<double> mean(nBands, 0.0);
                std::vector<double> sd(nBands, 0.0);
                int64_t count = 0;

                for (int64_t i = 0; i < total; ++i) {
                    if (static_cast<int>(labels[i]) != cls) continue;
                    if (hasND && (*bands[0])[i] == noData) continue;
                    count++;
                    for (int b = 0; b < nBands; ++b)
                        mean[b] += (*bands[b])[i];
                }

                if (count < 3) continue;

                for (int b = 0; b < nBands; ++b)
                    mean[b] /= count;

                for (int64_t i = 0; i < total; ++i) {
                    if (static_cast<int>(labels[i]) != cls) continue;
                    if (hasND && (*bands[0])[i] == noData) continue;
                    for (int b = 0; b < nBands; ++b) {
                        double d = (*bands[b])[i] - mean[b];
                        sd[b] += d * d;
                    }
                }
                for (int b = 0; b < nBands; ++b)
                    sd[b] = std::sqrt(sd[b] / (count - 1));

                // Remove outliers
                int64_t removedThisClass = 0;
                for (int64_t i = 0; i < total; ++i) {
                    if (static_cast<int>(labels[i]) != cls) continue;
                    if (hasND && (*bands[0])[i] == noData) continue;

                    bool outlier = false;
                    for (int b = 0; b < nBands; ++b) {
                        if (sd[b] < 1e-10) continue;
                        if (std::abs((*bands[b])[i] - mean[b]) > stdThresh * sd[b]) {
                            outlier = true;
                            break;
                        }
                    }

                    if (outlier) {
                        int64_t currentCount = count - removedThisClass;
                        if (static_cast<double>(currentCount - 1) / origCounts[ci] >= minRetain) {
                            labels[i] = 0.0;
                            removedThisClass++;
                            anyRemoved = true;
                        }
                    }
                }
            }

            if (!anyRemoved) break;

            reportProgress(0.10 + 0.80 * (iter + 1.0) / maxIter,
                QString("Iteration %1 complete").arg(iter + 1));
        }

        reportProgress(0.92, "Writing purified training output...");

        Raster output(cols, rows, 1, DataType::Int16);
        output.setGeoTransform(training->geoTransform());
        output.setProjection(training->projection());
        output.setNoDataValue(0);

        auto& dst = output.data(0);
        for (int64_t i = 0; i < total; ++i)
            dst[i] = labels[i];

        reportProgress(1.0, "Done.");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PurifyModule)

} // namespace aplaceholder
