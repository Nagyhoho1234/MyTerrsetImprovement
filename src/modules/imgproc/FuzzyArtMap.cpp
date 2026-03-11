#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <QStringList>

namespace aplaceholder {

class FuzzyArtMapModule : public Module {
public:
    QString name() const override { return "FUZZYARTMAP"; }
    QString description() const override {
        return "Fuzzy ARTMAP neural network classifier. Uses vigilance-controlled "
               "adaptive resonance for supervised classification of multispectral imagery.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Band images (comma-separated)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::file("training_raster", "Training site raster",
                "Raster with class labels for training pixels"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::real("vigilance", "Vigilance parameter", 0.75, 0.0, 1.0,
                "Controls granularity of category creation. Higher = more categories, "
                "more specific matching."),
            ParameterDef::real("learning_rate", "Learning rate", 1.0, 0.0, 1.0,
                "Rate at which neuron weights are updated (1.0 = fast learning)."),
        };
    }

    bool execute() override {
        // Parse band paths
        QStringList bandPaths = parameter("bands").toString().split(",", Qt::SkipEmptyParts);
        int numBands = bandPaths.size();
        if (numBands == 0) {
            setError("No band images specified");
            return false;
        }

        // Read training raster
        auto trainRaster = GdalIO::read(parameter("training_raster").toString());
        if (!trainRaster) {
            setError("Failed to read training raster");
            return false;
        }

        int cols = trainRaster->cols();
        int rows = trainRaster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double trainND = trainRaster->noDataValue();
        bool trainHasND = trainRaster->hasNoData();
        const auto& trainData = trainRaster->data(0);

        // Read band rasters
        std::vector<std::unique_ptr<Raster>> bandRasters;
        for (int b = 0; b < numBands; ++b) {
            auto br = GdalIO::read(bandPaths[b].trimmed());
            if (!br) {
                setError("Failed to read band raster: " + bandPaths[b].trimmed());
                return false;
            }
            if (br->cols() != cols || br->rows() != rows) {
                setError("Band dimensions do not match training raster: " +
                         bandPaths[b].trimmed());
                return false;
            }
            bandRasters.push_back(std::move(br));
        }

        double vigilance = parameter("vigilance").toDouble();
        double learningRate = parameter("learning_rate").toDouble();

        reportProgress(0.05, "Normalizing input data...");

        // Compute global min/max per band for normalization to [0, 1]
        std::vector<double> bandMin(numBands, 1e308);
        std::vector<double> bandMax(numBands, -1e308);
        for (int b = 0; b < numBands; ++b) {
            const auto& bd = bandRasters[b]->data(0);
            bool hasND = bandRasters[b]->hasNoData();
            double nd = bandRasters[b]->noDataValue();
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && bd[i] == nd) continue;
                if (bd[i] < bandMin[b]) bandMin[b] = bd[i];
                if (bd[i] > bandMax[b]) bandMax[b] = bd[i];
            }
        }

        // Complement coding: input vector becomes [x, 1-x] (dimension = 2*numBands)
        int M = 2 * numBands;

        // Helper to get normalized complement-coded pixel vector
        auto getPixelVector = [&](int64_t idx) -> std::vector<double> {
            std::vector<double> vec(M);
            for (int b = 0; b < numBands; ++b) {
                double val = bandRasters[b]->data(0)[idx];
                double range = bandMax[b] - bandMin[b];
                double norm = (range > 0.0) ? (val - bandMin[b]) / range : 0.0;
                norm = std::max(0.0, std::min(1.0, norm));
                vec[b] = norm;
                vec[b + numBands] = 1.0 - norm;
            }
            return vec;
        };

        // Fuzzy ARTMAP neuron structure
        struct Neuron {
            std::vector<double> weights; // weight vector (size M), initialized to 1
            int classLabel;
        };

        std::vector<Neuron> neurons;

        // Fuzzy AND (component-wise minimum)
        auto fuzzyAnd = [&](const std::vector<double>& a,
                            const std::vector<double>& b) -> std::vector<double> {
            std::vector<double> result(M);
            for (int j = 0; j < M; ++j)
                result[j] = std::min(a[j], b[j]);
            return result;
        };

        // L1 norm
        auto norm1 = [&](const std::vector<double>& v) -> double {
            double s = 0.0;
            for (int j = 0; j < M; ++j)
                s += v[j];
            return s;
        };

        double alpha = 0.001; // choice parameter (small positive constant)

        reportProgress(0.1, "Training Fuzzy ARTMAP network...");

        // Phase 1: Training - learn from training pixels
        // Collect training pixel indices
        std::vector<int64_t> trainPixels;
        for (int64_t i = 0; i < total; ++i) {
            if (trainHasND && trainData[i] == trainND) continue;
            // Check if any band is nodata
            bool skip = false;
            for (int b = 0; b < numBands; ++b) {
                if (bandRasters[b]->hasNoData() &&
                    bandRasters[b]->data(0)[i] == bandRasters[b]->noDataValue()) {
                    skip = true;
                    break;
                }
            }
            if (!skip) trainPixels.push_back(i);
        }

        int64_t numTrain = static_cast<int64_t>(trainPixels.size());
        if (numTrain == 0) {
            setError("No valid training pixels found");
            return false;
        }

        for (int64_t t = 0; t < numTrain; ++t) {
            int64_t idx = trainPixels[t];
            int targetClass = static_cast<int>(trainData[idx]);
            std::vector<double> input = getPixelVector(idx);
            double inputNorm = norm1(input);

            // Try to find a matching neuron
            // Compute choice function for each neuron, try in order of highest choice
            struct CandidateScore {
                int neuronIdx;
                double choiceVal;
            };
            std::vector<CandidateScore> candidates;
            for (int n = 0; n < static_cast<int>(neurons.size()); ++n) {
                std::vector<double> match = fuzzyAnd(input, neurons[n].weights);
                double matchNorm = norm1(match);
                double choiceVal = matchNorm / (alpha + norm1(neurons[n].weights));
                candidates.push_back({n, choiceVal});
            }

            // Sort by choice value descending
            std::sort(candidates.begin(), candidates.end(),
                      [](const CandidateScore& a, const CandidateScore& b) {
                          return a.choiceVal > b.choiceVal;
                      });

            bool resonance = false;
            for (auto& cand : candidates) {
                int n = cand.neuronIdx;
                std::vector<double> match = fuzzyAnd(input, neurons[n].weights);
                double matchNorm = norm1(match);

                // Vigilance test
                double matchRatio = (inputNorm > 0.0) ? matchNorm / inputNorm : 0.0;
                if (matchRatio >= vigilance) {
                    // Check class match (map field)
                    if (neurons[n].classLabel == targetClass) {
                        // Resonance: update weights
                        for (int j = 0; j < M; ++j) {
                            neurons[n].weights[j] = learningRate * match[j] +
                                                     (1.0 - learningRate) * neurons[n].weights[j];
                        }
                        resonance = true;
                        break;
                    }
                    // Match reset: class mismatch, raise vigilance and try next neuron
                }
            }

            if (!resonance) {
                // Create new neuron
                Neuron newNeuron;
                newNeuron.weights = input;
                newNeuron.classLabel = targetClass;
                neurons.push_back(std::move(newNeuron));
            }

            if (t % (numTrain / 20 + 1) == 0)
                reportProgress(0.1 + 0.5 * static_cast<double>(t) / numTrain,
                               "Training... (" + QString::number(neurons.size()) + " neurons)");
        }

        reportProgress(0.6, "Classifying with " + QString::number(neurons.size()) + " neurons...");

        // Phase 2: Classification - classify all pixels
        double outND = -9999.0;
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(trainRaster->geoTransform());
        output.setProjection(trainRaster->projection());
        output.setNoDataValue(outND);
        auto& outData = output.data(0);

        for (int64_t i = 0; i < total; ++i) {
            // Check nodata
            bool skip = false;
            for (int b = 0; b < numBands; ++b) {
                if (bandRasters[b]->hasNoData() &&
                    bandRasters[b]->data(0)[i] == bandRasters[b]->noDataValue()) {
                    skip = true;
                    break;
                }
            }
            if (skip) {
                outData[i] = outND;
                continue;
            }

            std::vector<double> input = getPixelVector(i);

            // Find best matching neuron (highest choice function)
            int bestNeuron = -1;
            double bestChoice = -1.0;
            for (int n = 0; n < static_cast<int>(neurons.size()); ++n) {
                std::vector<double> match = fuzzyAnd(input, neurons[n].weights);
                double matchNorm = norm1(match);
                double choiceVal = matchNorm / (alpha + norm1(neurons[n].weights));
                if (choiceVal > bestChoice) {
                    bestChoice = choiceVal;
                    bestNeuron = n;
                }
            }

            if (bestNeuron >= 0) {
                outData[i] = static_cast<double>(neurons[bestNeuron].classLabel);
            } else {
                outData[i] = outND;
            }

            if (i % (total / 20 + 1) == 0)
                reportProgress(0.6 + 0.35 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Classification complete.");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(FuzzyArtMapModule)

} // namespace aplaceholder
