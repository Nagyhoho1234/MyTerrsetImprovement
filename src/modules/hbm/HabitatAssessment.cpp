#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

#include <cmath>
#include <vector>
#include <memory>

namespace aplaceholder {

class HabitatAssessmentModule : public Module {
public:
    QString name() const override { return "HABITAT_ASSESSMENT"; }
    QString description() const override {
        return "Habitat suitability assessment using weighted linear combination (WLC) "
               "of multiple factor rasters. Each factor is assumed to be pre-standardized "
               "to a 0-1 suitability scale. The module computes the weighted sum and "
               "optionally applies a binary threshold for suitable/unsuitable classification.";
    }
    QString category() const override { return "Habitat & Biodiversity"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_rasters", "Input factor rasters (comma-separated paths)",
                "Suitability factor rasters (e.g. land cover suitability, distance to water, "
                "elevation suitability, slope suitability). Each raster should contain values "
                "in the range 0.0 to 1.0."),
            ParameterDef::file("weights", "Weights (comma-separated values)",
                "Comma-separated numeric weights for each input raster. "
                "Weights are normalized so they sum to 1.0."),
            ParameterDef::output("output", "Output suitability raster",
                "Output raster with weighted suitability scores (0-1) or binary "
                "suitable/unsuitable if threshold is applied"),
            ParameterDef::real("threshold", "Suitability threshold",
                0.0, 0.0, 1.0,
                "If > 0, output is binary: 1 = suitable (score >= threshold), "
                "0 = unsuitable. Set to 0 for continuous suitability output."),
        };
    }

    bool execute() override {
        // --- Parse parameters --------------------------------------------------
        const QString inputPathsStr = parameter("input_rasters").toString();
        const QString weightsStr    = parameter("weights").toString();
        const QString outputPath    = parameter("output").toString();
        const double threshold      = parameter("threshold").toDouble();

        if (inputPathsStr.isEmpty() || weightsStr.isEmpty() || outputPath.isEmpty()) {
            setError("Input rasters, weights, and output path are all required.");
            return false;
        }

        // Split comma-separated paths and weights
        QStringList inputPaths = inputPathsStr.split(",", Qt::SkipEmptyParts);
        QStringList weightTokens = weightsStr.split(",", Qt::SkipEmptyParts);

        for (auto& p : inputPaths) p = p.trimmed();
        for (auto& w : weightTokens) w = w.trimmed();

        if (inputPaths.size() != weightTokens.size()) {
            setError(QString("Number of input rasters (%1) does not match number of "
                             "weights (%2).")
                         .arg(inputPaths.size())
                         .arg(weightTokens.size()));
            return false;
        }

        const int numFactors = inputPaths.size();
        if (numFactors == 0) {
            setError("At least one input factor raster is required.");
            return false;
        }

        // Parse and normalize weights
        std::vector<double> weights(numFactors);
        double weightSum = 0.0;
        for (int i = 0; i < numFactors; ++i) {
            bool ok = false;
            weights[i] = weightTokens[i].toDouble(&ok);
            if (!ok || weights[i] < 0.0) {
                setError(QString("Invalid weight value: '%1'. Weights must be "
                                 "non-negative numbers.")
                             .arg(weightTokens[i]));
                return false;
            }
            weightSum += weights[i];
        }
        if (weightSum <= 0.0) {
            setError("Sum of weights must be greater than zero.");
            return false;
        }
        for (int i = 0; i < numFactors; ++i) {
            weights[i] /= weightSum;
        }

        // --- Load all input rasters -------------------------------------------
        reportProgress(0.0, "Loading input rasters...");

        std::vector<std::unique_ptr<Raster>> rasters(numFactors);
        rasters[0] = GdalIO::read(inputPaths[0]);
        if (!rasters[0]) {
            setError("Failed to read input raster: " + inputPaths[0]);
            return false;
        }

        const int cols = rasters[0]->cols();
        const int rows = rasters[0]->rows();
        const bool hasNoData = rasters[0]->hasNoData();
        const double noData  = rasters[0]->noDataValue();

        for (int i = 1; i < numFactors; ++i) {
            rasters[i] = GdalIO::read(inputPaths[i]);
            if (!rasters[i]) {
                setError("Failed to read input raster: " + inputPaths[i]);
                return false;
            }
            if (rasters[i]->cols() != cols || rasters[i]->rows() != rows) {
                setError(QString("Raster dimension mismatch: '%1' is %2x%3 but "
                                 "expected %4x%5.")
                             .arg(inputPaths[i])
                             .arg(rasters[i]->cols())
                             .arg(rasters[i]->rows())
                             .arg(cols)
                             .arg(rows));
                return false;
            }
        }

        reportProgress(0.1, "Input rasters loaded. Computing suitability...");

        // --- Create output raster ---------------------------------------------
        auto output = std::make_unique<Raster>(cols, rows, 1, DataType::Float32);
        output->setGeoTransform(rasters[0]->geoTransform());
        output->setProjection(rasters[0]->projection());
        output->setNoDataValue(noData);
        output->allocate();

        auto& outData = output->data(0);

        // Collect references to input data vectors for fast access
        std::vector<const std::vector<double>*> inData(numFactors);
        for (int i = 0; i < numFactors; ++i) {
            inData[i] = &rasters[i]->data(0);
        }

        const int64_t totalCells = static_cast<int64_t>(rows) * cols;
        const bool applyThreshold = (threshold > 0.0);

        // --- Weighted linear combination --------------------------------------
        for (int r = 0; r < rows; ++r) {
            if (r % 100 == 0) {
                reportProgress(0.1 + 0.85 * static_cast<double>(r) / rows,
                               QString("Processing row %1 / %2").arg(r + 1).arg(rows));
            }

            for (int c = 0; c < cols; ++c) {
                const int64_t idx = static_cast<int64_t>(r) * cols + c;

                // Check for nodata in any layer
                bool isNoData = false;
                for (int f = 0; f < numFactors; ++f) {
                    double val = (*inData[f])[idx];
                    if (hasNoData && val == noData) {
                        isNoData = true;
                        break;
                    }
                    // Also check per-raster nodata if different
                    if (rasters[f]->hasNoData() && val == rasters[f]->noDataValue()) {
                        isNoData = true;
                        break;
                    }
                }

                if (isNoData) {
                    outData[idx] = noData;
                    continue;
                }

                // Weighted linear combination: S = sum(wi * fi)
                double suitability = 0.0;
                for (int f = 0; f < numFactors; ++f) {
                    double val = (*inData[f])[idx];
                    // Clamp factor values to [0, 1]
                    if (val < 0.0) val = 0.0;
                    if (val > 1.0) val = 1.0;
                    suitability += weights[f] * val;
                }

                if (applyThreshold) {
                    outData[idx] = (suitability >= threshold) ? 1.0 : 0.0;
                } else {
                    outData[idx] = suitability;
                }
            }
        }

        // --- Write output ------------------------------------------------------
        reportProgress(0.95, "Writing output raster...");
        if (!GdalIO::write(*output, outputPath)) {
            setError("Failed to write output raster: " + outputPath);
            return false;
        }

        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(HabitatAssessmentModule)

} // namespace aplaceholder
