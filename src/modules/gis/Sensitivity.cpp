#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>

namespace aplaceholder {

class SensitivityModule : public Module {
public:
    QString name() const override { return "SENSITIVITY"; }
    QString description() const override {
        return "One-at-a-time sensitivity analysis. Varies each input weight "
               "by a specified percentage and measures the effect on "
               "a weighted sum model output.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("inputs", "Input raster images (comma-separated paths)"),
            ParameterDef::file("weights", "Weights for each input (comma-separated values)"),
            ParameterDef::output("output_report", "Output sensitivity report text file"),
            ParameterDef::real("variation_pct", "Variation percentage", 10.0, 0.1, 100.0,
                "Percentage by which to vary each weight (e.g., 10 means ±10%)"),
        };
    }

    bool execute() override {
        QStringList inputPaths = parameter("inputs").toString().split(",", Qt::SkipEmptyParts);
        QStringList weightStrs = parameter("weights").toString().split(",", Qt::SkipEmptyParts);
        double variationPct = parameter("variation_pct").toDouble();
        if (variationPct <= 0.0) variationPct = 10.0;

        if (inputPaths.size() != weightStrs.size()) {
            setError("Number of inputs must match number of weights");
            return false;
        }

        int numInputs = inputPaths.size();
        if (numInputs < 1) {
            setError("At least one input is required");
            return false;
        }

        // Parse weights
        std::vector<double> weights(numInputs);
        for (int i = 0; i < numInputs; ++i) {
            bool ok;
            weights[i] = weightStrs[i].trimmed().toDouble(&ok);
            if (!ok) {
                setError("Invalid weight value: " + weightStrs[i].trimmed());
                return false;
            }
        }

        reportProgress(0.0, "Loading input rasters...");

        // Load rasters
        std::vector<std::unique_ptr<Raster>> rasters(numInputs);
        for (int i = 0; i < numInputs; ++i) {
            rasters[i] = GdalIO::read(inputPaths[i].trimmed());
            if (!rasters[i]) {
                setError("Failed to read raster: " + inputPaths[i].trimmed());
                return false;
            }
        }

        // Validate dimensions
        int cols = rasters[0]->cols(), rows = rasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        for (int i = 1; i < numInputs; ++i) {
            if (rasters[i]->cols() != cols || rasters[i]->rows() != rows) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
        }

        reportProgress(0.1, "Computing baseline weighted sum...");

        // Compute baseline mean of weighted sum
        auto computeWeightedMean = [&](const std::vector<double>& w) -> double {
            double totalSum = 0.0;
            int64_t validCount = 0;
            for (int64_t p = 0; p < total; ++p) {
                bool valid = true;
                double val = 0.0;
                for (int i = 0; i < numInputs; ++i) {
                    const auto& data = rasters[i]->data(0);
                    if (rasters[i]->hasNoData() && data[p] == rasters[i]->noDataValue()) {
                        valid = false;
                        break;
                    }
                    val += w[i] * data[p];
                }
                if (valid) {
                    totalSum += val;
                    ++validCount;
                }
            }
            return (validCount > 0) ? totalSum / validCount : 0.0;
        };

        double baselineMean = computeWeightedMean(weights);

        reportProgress(0.3, "Running sensitivity analysis...");

        double variation = variationPct / 100.0;

        // For each input, vary its weight ±variation and measure change
        struct SensitivityResult {
            QString inputName;
            double weight;
            double plusMean;
            double minusMean;
            double plusChange;
            double minusChange;
            double sensitivity;
        };

        std::vector<SensitivityResult> results;
        for (int i = 0; i < numInputs; ++i) {
            std::vector<double> wPlus = weights;
            std::vector<double> wMinus = weights;
            wPlus[i] = weights[i] * (1.0 + variation);
            wMinus[i] = weights[i] * (1.0 - variation);

            double plusMean = computeWeightedMean(wPlus);
            double minusMean = computeWeightedMean(wMinus);

            double plusChange = (baselineMean != 0.0) ?
                (plusMean - baselineMean) / std::abs(baselineMean) * 100.0 : 0.0;
            double minusChange = (baselineMean != 0.0) ?
                (minusMean - baselineMean) / std::abs(baselineMean) * 100.0 : 0.0;
            double sensitivity = (plusChange - minusChange) / 2.0;

            results.push_back({inputPaths[i].trimmed(), weights[i],
                               plusMean, minusMean, plusChange, minusChange, sensitivity});

            reportProgress(0.3 + 0.6 * static_cast<double>(i + 1) / numInputs);
        }

        reportProgress(0.9, "Writing report...");

        // Write report
        QString reportPath = parameter("output_report").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        outFile << "One-at-a-Time Sensitivity Analysis\n";
        outFile << "====================================\n\n";
        outFile << "Variation: ±" << variationPct << "%\n";
        outFile << "Number of inputs: " << numInputs << "\n";
        outFile << "Baseline weighted mean: " << baselineMean << "\n\n";

        outFile << "Results:\n";
        outFile << "--------\n\n";

        for (int i = 0; i < numInputs; ++i) {
            const auto& res = results[i];
            outFile << "Input " << (i + 1) << ": " << res.inputName.toStdString() << "\n";
            outFile << "  Weight: " << res.weight << "\n";
            outFile << "  +" << variationPct << "% weight -> mean = " << res.plusMean
                    << " (change: " << res.plusChange << "%)\n";
            outFile << "  -" << variationPct << "% weight -> mean = " << res.minusMean
                    << " (change: " << res.minusChange << "%)\n";
            outFile << "  Sensitivity index: " << res.sensitivity << "\n\n";
        }

        // Rank by absolute sensitivity
        outFile << "Sensitivity Ranking (most to least sensitive):\n";
        std::vector<int> order(numInputs);
        for (int i = 0; i < numInputs; ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return std::abs(results[a].sensitivity) > std::abs(results[b].sensitivity);
        });
        for (int rank = 0; rank < numInputs; ++rank) {
            int i = order[rank];
            outFile << "  " << (rank + 1) << ". Input " << (i + 1) << " (sensitivity: "
                    << results[i].sensitivity << ")\n";
        }

        outFile.close();

        reportProgress(1.0, QString("Sensitivity analysis complete for %1 inputs").arg(numInputs));

        return true;
    }
};

REGISTER_MODULE(SensitivityModule)

} // namespace aplaceholder
