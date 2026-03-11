#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class MdChoiceModule : public Module {
public:
    QString name() const override { return "MDCHOICE"; }
    QString description() const override {
        return "Multi-dimensional choice model for decision support. "
               "Takes multiple suitability maps and corresponding weights, "
               "computes weighted suitability scores, and produces a ranked "
               "allocation raster assigning each pixel to the highest-scoring alternative.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("suitability_maps", "Suitability raster images (semicolon-separated list)"),
            ParameterDef::file("weights", "Weight values per map (semicolon-separated, must sum to 1.0)"),
            ParameterDef::output("output_allocation", "Output ranked allocation raster"),
            ParameterDef::output("output_score", "Output best weighted suitability score raster"),
        };
    }

    bool execute() override {
        // Parse semicolon-separated suitability map paths
        QString mapsStr = parameter("suitability_maps").toString();
        QStringList mapPaths = mapsStr.split(";", Qt::SkipEmptyParts);

        // Parse semicolon-separated weights
        QString weightsStr = parameter("weights").toString();
        QStringList weightTokens = weightsStr.split(";", Qt::SkipEmptyParts);

        int numAlternatives = mapPaths.size();
        if (numAlternatives < 2) {
            setError("At least 2 suitability maps are required");
            return false;
        }

        if (weightTokens.size() != numAlternatives) {
            setError("Number of weights must match number of suitability maps");
            return false;
        }

        // Parse weights
        std::vector<double> weights(numAlternatives);
        double weightSum = 0.0;
        for (int k = 0; k < numAlternatives; ++k) {
            bool ok = false;
            weights[k] = weightTokens[k].trimmed().toDouble(&ok);
            if (!ok || weights[k] < 0.0) {
                setError("Invalid weight value: " + weightTokens[k]);
                return false;
            }
            weightSum += weights[k];
        }

        if (std::abs(weightSum - 1.0) > 0.01) {
            setError(QString("Weights must sum to 1.0 (current sum = %1)").arg(weightSum));
            return false;
        }

        reportProgress(0.0, "Reading suitability maps...");

        // Read all suitability rasters
        std::vector<std::shared_ptr<Raster>> rasters(numAlternatives);
        for (int k = 0; k < numAlternatives; ++k) {
            rasters[k] = GdalIO::read(mapPaths[k].trimmed());
            if (!rasters[k]) {
                setError("Failed to read suitability map: " + mapPaths[k].trimmed());
                return false;
            }
        }

        // Verify all rasters have the same dimensions
        int cols = rasters[0]->cols(), rows = rasters[0]->rows();
        for (int k = 1; k < numAlternatives; ++k) {
            if (rasters[k]->cols() != cols || rasters[k]->rows() != rows) {
                setError("All suitability maps must have the same dimensions");
                return false;
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;

        reportProgress(0.3, "Computing weighted scores and allocation...");

        double outNoData = -9999.0;

        Raster allocation(cols, rows, 1, DataType::Float64);
        allocation.setGeoTransform(rasters[0]->geoTransform());
        allocation.setProjection(rasters[0]->projection());
        allocation.setNoDataValue(outNoData);
        auto& allocData = allocation.data(0);

        Raster score(cols, rows, 1, DataType::Float64);
        score.setGeoTransform(rasters[0]->geoTransform());
        score.setProjection(rasters[0]->projection());
        score.setNoDataValue(outNoData);
        auto& scoreData = score.data(0);

        for (int64_t i = 0; i < total; ++i) {
            bool isNoData = false;

            // Check if any raster has NoData at this pixel
            for (int k = 0; k < numAlternatives; ++k) {
                if (rasters[k]->hasNoData() && rasters[k]->data(0)[i] == rasters[k]->noDataValue()) {
                    isNoData = true;
                    break;
                }
            }

            if (isNoData) {
                allocData[i] = outNoData;
                scoreData[i] = outNoData;
                continue;
            }

            // Find the alternative with the highest weighted suitability
            double bestScore = -1e300;
            int bestAlt = 0;

            for (int k = 0; k < numAlternatives; ++k) {
                double weighted = weights[k] * rasters[k]->data(0)[i];
                if (weighted > bestScore) {
                    bestScore = weighted;
                    bestAlt = k + 1; // 1-based allocation class
                }
            }

            allocData[i] = static_cast<double>(bestAlt);
            scoreData[i] = bestScore;
        }

        reportProgress(0.7, "Writing allocation raster...");

        if (!GdalIO::write(allocation, parameter("output_allocation").toString())) {
            setError("Failed to write allocation raster");
            return false;
        }

        reportProgress(0.9, "Writing score raster...");

        if (!GdalIO::write(score, parameter("output_score").toString())) {
            setError("Failed to write score raster");
            return false;
        }

        reportProgress(1.0,
            QString("Allocation complete: %1 alternatives, %2 pixels processed")
                .arg(numAlternatives)
                .arg(total));

        return true;
    }
};

REGISTER_MODULE(MdChoiceModule)

} // namespace aplaceholder
