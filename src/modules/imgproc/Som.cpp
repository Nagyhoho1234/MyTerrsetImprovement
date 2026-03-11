#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>

namespace aplaceholder {

class SomModule : public Module {
public:
    QString name() const override { return "SOM"; }
    QString description() const override {
        return "Self-Organizing Map (Kohonen) classifier. 2D grid of neurons with competitive "
               "learning; after training, each neuron is assigned to its dominant class from "
               "training data.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster where pixel values indicate class IDs (0 = unclassified)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("grid_width", "SOM grid width", 10, 2, 100,
                "Width of the 2D SOM neuron grid"),
            ParameterDef::integer("grid_height", "SOM grid height", 10, 2, 100,
                "Height of the 2D SOM neuron grid"),
            ParameterDef::integer("max_epochs", "Maximum epochs", 100, 1, 10000,
                "Maximum number of training epochs"),
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

        int gridW = parameter("grid_width").toInt();
        int gridH = parameter("grid_height").toInt();
        int maxEpochs = parameter("max_epochs").toInt();
        int numNeurons = gridW * gridH;

        // --------------------------------------------------------------------
        // 3. Collect training samples
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        const auto& trainData = trainRaster->data(0);

        struct Sample {
            std::vector<double> features;
            int classId;
        };

        std::vector<Sample> samples;
        std::vector<int> classIds;

        for (int64_t i = 0; i < total; ++i) {
            int cls = static_cast<int>(trainData[i]);
            if (cls <= 0) continue;
            if (hasND && (*bands[0])[i] == noData) continue;

            Sample s;
            s.features.resize(numBands);
            for (int b = 0; b < numBands; ++b)
                s.features[b] = (*bands[b])[i];
            s.classId = cls;
            samples.push_back(std::move(s));

            if (std::find(classIds.begin(), classIds.end(), cls) == classIds.end())
                classIds.push_back(cls);
        }

        std::sort(classIds.begin(), classIds.end());
        int numClasses = static_cast<int>(classIds.size());
        int nSamples = static_cast<int>(samples.size());

        if (numClasses < 2) {
            setError("Need at least 2 classes in training data");
            return false;
        }

        reportProgress(0.05, "Initializing SOM grid...");

        // --------------------------------------------------------------------
        // 4. Initialize SOM neuron weight vectors
        // --------------------------------------------------------------------
        std::mt19937 rng(42);

        // Compute data range for initialization
        std::vector<double> bandMin(numBands, std::numeric_limits<double>::max());
        std::vector<double> bandMax(numBands, -std::numeric_limits<double>::max());
        for (const auto& s : samples) {
            for (int b = 0; b < numBands; ++b) {
                bandMin[b] = std::min(bandMin[b], s.features[b]);
                bandMax[b] = std::max(bandMax[b], s.features[b]);
            }
        }

        // Neuron weights: neurons[neuronIdx][band]
        std::vector<std::vector<double>> neurons(numNeurons, std::vector<double>(numBands));
        for (int n = 0; n < numNeurons; ++n) {
            for (int b = 0; b < numBands; ++b) {
                std::uniform_real_distribution<double> dist(bandMin[b], bandMax[b]);
                neurons[n][b] = dist(rng);
            }
        }

        // Grid positions for neighborhood calculation
        auto neuronRow = [gridW](int idx) { return idx / gridW; };
        auto neuronCol = [gridW](int idx) { return idx % gridW; };

        // --------------------------------------------------------------------
        // 5. Train SOM with competitive learning
        // --------------------------------------------------------------------
        double initLR = 0.5;
        double initRadius = std::max(gridW, gridH) / 2.0;
        double timeConst = maxEpochs / std::log(initRadius + 1.0);

        std::vector<int> sampleOrder(nSamples);
        std::iota(sampleOrder.begin(), sampleOrder.end(), 0);

        for (int epoch = 0; epoch < maxEpochs; ++epoch) {
            double lr = initLR * std::exp(-static_cast<double>(epoch) / maxEpochs);
            double radius = initRadius * std::exp(-static_cast<double>(epoch) / timeConst);
            double radiusSq = radius * radius;
            if (radiusSq < 1e-10) radiusSq = 1e-10;

            std::shuffle(sampleOrder.begin(), sampleOrder.end(), rng);

            for (int si = 0; si < nSamples; ++si) {
                const auto& sample = samples[sampleOrder[si]];

                // Find Best Matching Unit (BMU)
                double bestDist = std::numeric_limits<double>::max();
                int bmu = 0;
                for (int n = 0; n < numNeurons; ++n) {
                    double distSq = 0.0;
                    for (int b = 0; b < numBands; ++b) {
                        double diff = sample.features[b] - neurons[n][b];
                        distSq += diff * diff;
                    }
                    if (distSq < bestDist) {
                        bestDist = distSq;
                        bmu = n;
                    }
                }

                // Update BMU and its neighbors
                int bmuR = neuronRow(bmu);
                int bmuC = neuronCol(bmu);

                for (int n = 0; n < numNeurons; ++n) {
                    int nr = neuronRow(n);
                    int nc = neuronCol(n);
                    double gridDistSq = static_cast<double>((nr - bmuR) * (nr - bmuR) +
                                                            (nc - bmuC) * (nc - bmuC));

                    if (gridDistSq > radiusSq * 4.0) continue; // skip distant neurons

                    double influence = std::exp(-gridDistSq / (2.0 * radiusSq));
                    double updateRate = lr * influence;

                    for (int b = 0; b < numBands; ++b)
                        neurons[n][b] += updateRate * (sample.features[b] - neurons[n][b]);
                }
            }

            if (epoch % 10 == 0)
                reportProgress(0.05 + 0.6 * static_cast<double>(epoch) / maxEpochs,
                               QString("SOM epoch %1 of %2").arg(epoch + 1).arg(maxEpochs));
        }

        reportProgress(0.65, "Labeling SOM neurons from training data...");

        // --------------------------------------------------------------------
        // 6. Label each neuron by dominant class from training samples
        // --------------------------------------------------------------------
        // For each neuron, count how many training samples of each class map to it
        std::vector<std::vector<int>> neuronClassCounts(numNeurons, std::vector<int>(numClasses, 0));

        for (int i = 0; i < nSamples; ++i) {
            double bestDist = std::numeric_limits<double>::max();
            int bmu = 0;
            for (int n = 0; n < numNeurons; ++n) {
                double distSq = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = samples[i].features[b] - neurons[n][b];
                    distSq += diff * diff;
                }
                if (distSq < bestDist) {
                    bestDist = distSq;
                    bmu = n;
                }
            }

            int cidx = static_cast<int>(std::find(classIds.begin(), classIds.end(),
                                                    samples[i].classId) - classIds.begin());
            neuronClassCounts[bmu][cidx]++;
        }

        // Assign each neuron its dominant class (or 0 if no samples mapped)
        std::vector<int> neuronLabel(numNeurons, 0);
        for (int n = 0; n < numNeurons; ++n) {
            int bestCount = 0;
            int bestClass = 0;
            for (int ci = 0; ci < numClasses; ++ci) {
                if (neuronClassCounts[n][ci] > bestCount) {
                    bestCount = neuronClassCounts[n][ci];
                    bestClass = classIds[ci];
                }
            }
            neuronLabel[n] = bestClass;
        }

        // For neurons with no training samples, assign nearest labeled neuron's class
        for (int n = 0; n < numNeurons; ++n) {
            if (neuronLabel[n] != 0) continue;
            double bestDist = std::numeric_limits<double>::max();
            for (int n2 = 0; n2 < numNeurons; ++n2) {
                if (neuronLabel[n2] == 0) continue;
                double distSq = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = neurons[n][b] - neurons[n2][b];
                    distSq += diff * diff;
                }
                if (distSq < bestDist) {
                    bestDist = distSq;
                    neuronLabel[n] = neuronLabel[n2];
                }
            }
        }

        // --------------------------------------------------------------------
        // 7. Classify all pixels
        // --------------------------------------------------------------------
        reportProgress(0.7, "Classifying pixels...");

        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);
        std::fill(out.begin(), out.end(), 0.0);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                out[i] = 0;
                continue;
            }

            double bestDist = std::numeric_limits<double>::max();
            int bmu = 0;
            for (int n = 0; n < numNeurons; ++n) {
                double distSq = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = (*bands[b])[i] - neurons[n][b];
                    distSq += diff * diff;
                }
                if (distSq < bestDist) {
                    bestDist = distSq;
                    bmu = n;
                }
            }

            out[i] = neuronLabel[bmu];

            if (i % 1000000 == 0)
                reportProgress(0.7 + 0.25 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SomModule)

} // namespace aplaceholder
