#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <memory>
#include <numeric>
#include <algorithm>

namespace aplaceholder {

class SimWeightModule : public Module {
public:
    QString name() const override { return "SIM_WEIGHT"; }
    QString description() const override {
        return "Similarity-weighted combination of multiple input rasters. Each input "
               "is weighted by its spatial similarity (Pearson correlation) to a reference "
               "raster. Higher-similarity inputs receive more weight in the combination. "
               "Output = sum(w_i * input_i) / sum(w_i) where w_i = max(0, corr(input_i, ref)).";
    }
    QString category() const override { return "Habitat & Biodiversity"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("inputs", "Input rasters (comma-separated)",
                "Comma-separated paths to continuous rasters to combine"),
            ParameterDef::file("reference", "Reference raster",
                "Reference raster used to compute similarity weights"),
            ParameterDef::output("output", "Output weighted combination raster",
                "Continuous raster of similarity-weighted combination"),
        };
    }

    bool execute() override {
        // 1. Load reference raster
        auto reference = GdalIO::read(parameter("reference").toString());
        if (!reference) {
            setError("Failed to read reference raster.");
            return false;
        }

        int cols = reference->cols();
        int rows = reference->rows();
        int64_t total = reference->cellCount();

        reportProgress(0.03, "Reference raster loaded.");

        // 2. Parse and load input rasters
        QString inputsStr = parameter("inputs").toString();
        QStringList inputPaths = inputsStr.split(',', Qt::SkipEmptyParts);
        for (int i = 0; i < inputPaths.size(); ++i)
            inputPaths[i] = inputPaths[i].trimmed();

        int numInputs = inputPaths.size();
        if (numInputs < 1) {
            setError("At least one input raster is required.");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> inputs;
        for (int i = 0; i < numInputs; ++i) {
            auto r = GdalIO::read(inputPaths[i]);
            if (!r) {
                setError("Failed to read input raster: " + inputPaths[i]);
                return false;
            }
            if (r->cols() != cols || r->rows() != rows) {
                setError(QString("Dimension mismatch: reference is %1x%2 but input %3 is %4x%5")
                             .arg(cols).arg(rows).arg(i)
                             .arg(r->cols()).arg(r->rows()));
                return false;
            }
            inputs.push_back(std::move(r));
        }

        reportProgress(0.10, QString("%1 input rasters loaded.").arg(numInputs));

        // 3. Compute Pearson correlation between each input and the reference
        //    Use only pixels that are valid across ALL rasters
        double noData = reference->noDataValue();
        bool hasND_ref = reference->hasNoData();
        double ndRef = reference->noDataValue();

        std::vector<const std::vector<double>*> iData;
        std::vector<bool> hasNDVec;
        std::vector<double> ndVec;
        for (int i = 0; i < numInputs; ++i) {
            iData.push_back(&inputs[i]->data(0));
            hasNDVec.push_back(inputs[i]->hasNoData());
            ndVec.push_back(inputs[i]->noDataValue());
        }

        const auto& refData = reference->data(0);

        // Build valid pixel mask
        std::vector<bool> valid(total, true);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND_ref && refData[i] == ndRef) {
                valid[i] = false;
                continue;
            }
            for (int r = 0; r < numInputs; ++r) {
                if (hasNDVec[r] && (*iData[r])[i] == ndVec[r]) {
                    valid[i] = false;
                    break;
                }
            }
        }

        // Compute mean of reference over valid pixels
        double sumRef = 0.0;
        int64_t validCount = 0;
        for (int64_t i = 0; i < total; ++i) {
            if (!valid[i]) continue;
            sumRef += refData[i];
            ++validCount;
        }

        if (validCount == 0) {
            setError("No valid pixels found across all inputs and reference.");
            return false;
        }

        double meanRef = sumRef / validCount;

        // Compute correlation weights
        std::vector<double> weights(numInputs, 0.0);

        for (int r = 0; r < numInputs; ++r) {
            reportProgress(0.10 + 0.40 * static_cast<double>(r) / numInputs,
                           QString("Computing similarity for input %1/%2...").arg(r + 1).arg(numInputs));

            const auto& d = *iData[r];

            // Compute mean of input r
            double sumI = 0.0;
            for (int64_t i = 0; i < total; ++i) {
                if (!valid[i]) continue;
                sumI += d[i];
            }
            double meanI = sumI / validCount;

            // Pearson correlation
            double sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
            for (int64_t i = 0; i < total; ++i) {
                if (!valid[i]) continue;
                double dx = d[i] - meanI;
                double dy = refData[i] - meanRef;
                sumXY += dx * dy;
                sumX2 += dx * dx;
                sumY2 += dy * dy;
            }

            double denom = std::sqrt(sumX2 * sumY2);
            double corr = (denom > 0.0) ? sumXY / denom : 0.0;

            // Use only positive correlations as weights
            weights[r] = std::max(0.0, corr);
        }

        reportProgress(0.55, "Similarity weights computed.");

        // Log weights
        for (int r = 0; r < numInputs; ++r) {
            reportProgress(0.55, QString("  Input %1 weight: %2").arg(r).arg(weights[r], 0, 'f', 4));
        }

        // Check if all weights are zero
        double totalWeight = 0.0;
        for (int r = 0; r < numInputs; ++r)
            totalWeight += weights[r];

        if (totalWeight <= 0.0) {
            // Fall back to equal weights
            for (int r = 0; r < numInputs; ++r)
                weights[r] = 1.0;
            totalWeight = static_cast<double>(numInputs);
            reportProgress(0.56, "All correlations non-positive; using equal weights.");
        }

        // 4. Compute weighted combination
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(reference->geoTransform());
        output.setProjection(reference->projection());
        output.setNoDataValue(noData);

        auto& outData = output.data(0);

        for (int64_t i = 0; i < total; ++i) {
            if (!valid[i]) {
                outData[i] = noData;
                continue;
            }

            double weightedSum = 0.0;
            for (int r = 0; r < numInputs; ++r)
                weightedSum += weights[r] * (*iData[r])[i];

            outData[i] = weightedSum / totalWeight;

            if (i % 500000 == 0)
                reportProgress(0.55 + 0.35 * static_cast<double>(i) / total);
        }

        reportProgress(0.92, "Writing output...");

        // 5. Write output
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        auto stats = output.computeStats(0);
        reportProgress(1.0,
            QString("Done. Output — min: %1, max: %2, mean: %3, valid pixels: %4")
                .arg(stats.min, 0, 'f', 4)
                .arg(stats.max, 0, 'f', 4)
                .arg(stats.mean, 0, 'f', 4)
                .arg(validCount));

        return true;
    }
};

REGISTER_MODULE(SimWeightModule)

} // namespace aplaceholder
