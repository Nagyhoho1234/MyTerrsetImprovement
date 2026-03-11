#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class KnnModule : public Module {
public:
    QString name() const override { return "KNN"; }
    QString description() const override {
        return "K-Nearest Neighbor classifier. Classifies each pixel based on the majority "
               "class among its K nearest training examples in band space (Euclidean distance).";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster where pixel values indicate class IDs (0 = unclassified)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("k", "K (number of neighbors)", 5, 1, 1000,
                "Number of nearest neighbors used for majority voting"),
        };
    }

    bool execute() override {
        // --------------------------------------------------------------------
        // 1. Read input bands
        // --------------------------------------------------------------------
        QString bandsParam = parameter("input_bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths) p = p.trimmed();

        if (bandPaths.isEmpty()) {
            setError("No input bands specified");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> bandRasters;
        for (const auto& path : bandPaths) {
            auto r = GdalIO::read(path);
            if (!r) {
                setError("Failed to read band image: " + path);
                return false;
            }
            bandRasters.push_back(std::move(r));
        }

        int numBands = static_cast<int>(bandRasters.size());
        int cols = bandRasters[0]->cols();
        int rows = bandRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = bandRasters[0]->hasNoData();
        double noData = bandRasters[0]->noDataValue();

        for (int b = 1; b < numBands; ++b) {
            if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                setError("All input bands must have the same dimensions");
                return false;
            }
        }

        // --------------------------------------------------------------------
        // 2. Read training raster
        // --------------------------------------------------------------------
        auto trainRaster = GdalIO::read(parameter("training_raster").toString());
        if (!trainRaster) {
            setError("Failed to read training raster");
            return false;
        }
        if (trainRaster->cols() != cols || trainRaster->rows() != rows) {
            setError("Training raster dimensions must match input bands");
            return false;
        }

        int K = parameter("k").toInt();

        // --------------------------------------------------------------------
        // 3. Collect training samples
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        const auto& trainData = trainRaster->data(0);

        struct TrainSample {
            std::vector<double> features;
            int classId;
        };

        std::vector<TrainSample> trainSamples;

        for (int64_t i = 0; i < total; ++i) {
            int cls = static_cast<int>(trainData[i]);
            if (cls <= 0) continue;
            if (hasND && (*bands[0])[i] == noData) continue;

            TrainSample s;
            s.features.resize(numBands);
            for (int b = 0; b < numBands; ++b)
                s.features[b] = (*bands[b])[i];
            s.classId = cls;
            trainSamples.push_back(std::move(s));
        }

        int nTrain = static_cast<int>(trainSamples.size());
        if (nTrain == 0) {
            setError("No valid training samples found");
            return false;
        }

        if (K > nTrain) K = nTrain;

        reportProgress(0.1, "Classifying pixels using KNN...");

        // --------------------------------------------------------------------
        // 4. Classify each pixel by K-nearest neighbor majority vote
        // --------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);
        std::fill(out.begin(), out.end(), 0.0);

        // Use a partial sort approach for K nearest neighbors
        struct DistClass {
            double dist;
            int classId;
        };

        std::vector<DistClass> distances(nTrain);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                out[i] = 0;
                continue;
            }

            // Compute squared Euclidean distance to all training samples
            for (int t = 0; t < nTrain; ++t) {
                double distSq = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = (*bands[b])[i] - trainSamples[t].features[b];
                    distSq += diff * diff;
                }
                distances[t].dist = distSq;
                distances[t].classId = trainSamples[t].classId;
            }

            // Partial sort to find K nearest
            std::partial_sort(distances.begin(), distances.begin() + K, distances.end(),
                [](const DistClass& a, const DistClass& b) { return a.dist < b.dist; });

            // Majority vote among K nearest
            std::vector<std::pair<int, int>> votes; // (classId, count)
            for (int k = 0; k < K; ++k) {
                int cls = distances[k].classId;
                bool found = false;
                for (auto& v : votes) {
                    if (v.first == cls) {
                        v.second++;
                        found = true;
                        break;
                    }
                }
                if (!found)
                    votes.push_back({cls, 1});
            }

            int bestClass = votes[0].first;
            int bestCount = votes[0].second;
            for (size_t v = 1; v < votes.size(); ++v) {
                if (votes[v].second > bestCount) {
                    bestCount = votes[v].second;
                    bestClass = votes[v].first;
                }
            }

            out[i] = bestClass;

            if (i % 100000 == 0)
                reportProgress(0.1 + 0.85 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(KnnModule)

} // namespace aplaceholder
