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

class MlpModule : public Module {
public:
    QString name() const override { return "MLP"; }
    QString description() const override {
        return "Multi-Layer Perceptron neural network classifier. Single hidden layer with "
               "sigmoid activation and backpropagation training.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster where pixel values indicate class IDs (0 = unclassified)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("hidden_nodes", "Number of hidden layer nodes", 20, 1, 1000,
                "Number of neurons in the hidden layer"),
            ParameterDef::real("learning_rate", "Learning rate", 0.01, 0.0001, 10.0,
                "Step size for weight updates during backpropagation"),
            ParameterDef::integer("max_epochs", "Maximum epochs", 500, 1, 100000,
                "Maximum number of training iterations over all samples"),
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

        int nHidden = parameter("hidden_nodes").toInt();
        double lr = parameter("learning_rate").toDouble();
        int maxEpochs = parameter("max_epochs").toInt();

        // --------------------------------------------------------------------
        // 3. Collect training samples and normalize
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

        if (numClasses < 2) {
            setError("Need at least 2 classes in training data");
            return false;
        }

        int nSamples = static_cast<int>(samples.size());

        // Compute normalization: min-max to [0, 1]
        std::vector<double> bandMin(numBands, std::numeric_limits<double>::max());
        std::vector<double> bandMax(numBands, -std::numeric_limits<double>::max());

        for (const auto& s : samples) {
            for (int b = 0; b < numBands; ++b) {
                bandMin[b] = std::min(bandMin[b], s.features[b]);
                bandMax[b] = std::max(bandMax[b], s.features[b]);
            }
        }

        std::vector<double> bandRange(numBands);
        for (int b = 0; b < numBands; ++b) {
            bandRange[b] = bandMax[b] - bandMin[b];
            if (bandRange[b] < 1e-10) bandRange[b] = 1.0;
        }

        // Normalize training samples
        for (auto& s : samples) {
            for (int b = 0; b < numBands; ++b)
                s.features[b] = (s.features[b] - bandMin[b]) / bandRange[b];
        }

        // Build class index map
        std::vector<int> classToIndex(classIds.back() + 1, -1);
        for (int ci = 0; ci < numClasses; ++ci)
            classToIndex[classIds[ci]] = ci;

        reportProgress(0.05, "Training MLP neural network...");

        // --------------------------------------------------------------------
        // 4. Initialize network weights
        // --------------------------------------------------------------------
        // Architecture: numBands -> nHidden -> numClasses
        // Weights: W1[nHidden][numBands], b1[nHidden]
        //          W2[numClasses][nHidden], b2[numClasses]

        std::mt19937 rng(42);
        std::uniform_real_distribution<double> initDist(-0.5, 0.5);

        auto randInit = [&]() { return initDist(rng); };

        // Hidden layer weights
        std::vector<std::vector<double>> W1(nHidden, std::vector<double>(numBands));
        std::vector<double> b1(nHidden, 0.0);
        for (int h = 0; h < nHidden; ++h) {
            for (int b = 0; b < numBands; ++b)
                W1[h][b] = randInit();
            b1[h] = randInit();
        }

        // Output layer weights
        std::vector<std::vector<double>> W2(numClasses, std::vector<double>(nHidden));
        std::vector<double> b2(numClasses, 0.0);
        for (int c = 0; c < numClasses; ++c) {
            for (int h = 0; h < nHidden; ++h)
                W2[c][h] = randInit();
            b2[c] = randInit();
        }

        auto sigmoid = [](double x) -> double {
            return 1.0 / (1.0 + std::exp(-std::clamp(x, -500.0, 500.0)));
        };

        // --------------------------------------------------------------------
        // 5. Train with backpropagation (stochastic gradient descent)
        // --------------------------------------------------------------------
        std::vector<int> indices(nSamples);
        std::iota(indices.begin(), indices.end(), 0);

        std::vector<double> hiddenOut(nHidden);
        std::vector<double> outputOut(numClasses);
        std::vector<double> outputDelta(numClasses);
        std::vector<double> hiddenDelta(nHidden);
        std::vector<double> target(numClasses, 0.0);

        for (int epoch = 0; epoch < maxEpochs; ++epoch) {
            std::shuffle(indices.begin(), indices.end(), rng);

            for (int si = 0; si < nSamples; ++si) {
                const auto& sample = samples[indices[si]];

                // Forward pass: hidden layer
                for (int h = 0; h < nHidden; ++h) {
                    double sum = b1[h];
                    for (int b = 0; b < numBands; ++b)
                        sum += W1[h][b] * sample.features[b];
                    hiddenOut[h] = sigmoid(sum);
                }

                // Forward pass: output layer
                for (int c = 0; c < numClasses; ++c) {
                    double sum = b2[c];
                    for (int h = 0; h < nHidden; ++h)
                        sum += W2[c][h] * hiddenOut[h];
                    outputOut[c] = sigmoid(sum);
                }

                // Target vector: 1 for correct class, 0 for others
                std::fill(target.begin(), target.end(), 0.0);
                target[classToIndex[sample.classId]] = 1.0;

                // Output layer deltas
                for (int c = 0; c < numClasses; ++c) {
                    double err = target[c] - outputOut[c];
                    outputDelta[c] = err * outputOut[c] * (1.0 - outputOut[c]);
                }

                // Hidden layer deltas
                for (int h = 0; h < nHidden; ++h) {
                    double errSum = 0.0;
                    for (int c = 0; c < numClasses; ++c)
                        errSum += outputDelta[c] * W2[c][h];
                    hiddenDelta[h] = errSum * hiddenOut[h] * (1.0 - hiddenOut[h]);
                }

                // Update output layer weights
                for (int c = 0; c < numClasses; ++c) {
                    for (int h = 0; h < nHidden; ++h)
                        W2[c][h] += lr * outputDelta[c] * hiddenOut[h];
                    b2[c] += lr * outputDelta[c];
                }

                // Update hidden layer weights
                for (int h = 0; h < nHidden; ++h) {
                    for (int b = 0; b < numBands; ++b)
                        W1[h][b] += lr * hiddenDelta[h] * sample.features[b];
                    b1[h] += lr * hiddenDelta[h];
                }
            }

            if (epoch % 50 == 0)
                reportProgress(0.05 + 0.65 * static_cast<double>(epoch) / maxEpochs,
                               QString("Epoch %1 of %2").arg(epoch + 1).arg(maxEpochs));
        }

        // --------------------------------------------------------------------
        // 6. Classify all pixels
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

            // Normalize pixel
            std::vector<double> pixel(numBands);
            for (int b = 0; b < numBands; ++b)
                pixel[b] = ((*bands[b])[i] - bandMin[b]) / bandRange[b];

            // Forward pass: hidden layer
            for (int h = 0; h < nHidden; ++h) {
                double sum = b1[h];
                for (int b = 0; b < numBands; ++b)
                    sum += W1[h][b] * pixel[b];
                hiddenOut[h] = sigmoid(sum);
            }

            // Forward pass: output layer
            double bestVal = -1.0;
            int bestClass = 0;
            for (int c = 0; c < numClasses; ++c) {
                double sum = b2[c];
                for (int h = 0; h < nHidden; ++h)
                    sum += W2[c][h] * hiddenOut[h];
                double val = sigmoid(sum);
                if (val > bestVal) {
                    bestVal = val;
                    bestClass = classIds[c];
                }
            }

            out[i] = bestClass;

            if (i % 1000000 == 0)
                reportProgress(0.7 + 0.25 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(MlpModule)

} // namespace aplaceholder
