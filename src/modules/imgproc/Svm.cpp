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

class SvmModule : public Module {
public:
    QString name() const override { return "SVM"; }
    QString description() const override {
        return "Support Vector Machine classifier. Uses a simplified linear SVM with "
               "one-vs-rest approach to classify pixels based on training site data.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster where pixel values indicate class IDs (0 = unclassified)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::real("C_param", "C (regularization parameter)", 1.0, 0.001, 999999.0,
                "Controls tolerance for misclassification; higher C = narrower margins"),
            ParameterDef::integer("max_iterations", "Maximum iterations (SMO)", 100, 1, 100000,
                "Maximum iterations for the simplified SMO training loop"),
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

        double C = parameter("C_param").toDouble();
        int maxIter = parameter("max_iterations").toInt();

        // --------------------------------------------------------------------
        // 3. Collect band data and training samples
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        const auto& trainData = trainRaster->data(0);

        // Collect training samples: feature vectors and class labels
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

        reportProgress(0.1, "Training SVM classifiers (one-vs-rest)...");

        // --------------------------------------------------------------------
        // 4. Train one-vs-rest linear SVMs using simplified SMO
        // --------------------------------------------------------------------
        // For each class c, train a binary SVM: class c = +1, all others = -1
        // Linear SVM: decision function = w . x + bias
        // Simplified SMO: iterate over all alphas, check KKT, update pairs

        int nSamples = static_cast<int>(samples.size());

        struct SvmModel {
            int classId;
            std::vector<double> weights; // numBands
            double bias;
        };

        std::vector<SvmModel> models;
        std::mt19937 rng(42);

        for (int ci = 0; ci < numClasses; ++ci) {
            int targetClass = classIds[ci];

            // Binary labels: +1 for target class, -1 for others
            std::vector<double> y(nSamples);
            for (int i = 0; i < nSamples; ++i)
                y[i] = (samples[i].classId == targetClass) ? 1.0 : -1.0;

            // Alphas (Lagrange multipliers)
            std::vector<double> alpha(nSamples, 0.0);
            double bias = 0.0;

            // Precompute linear kernel cache would be too large; compute on-the-fly
            // Simplified SMO: iterate, pick pairs, update
            double tol = 1e-3;

            for (int iter = 0; iter < maxIter; ++iter) {
                int numChanged = 0;

                for (int i = 0; i < nSamples; ++i) {
                    // Compute f(x_i) = sum_j alpha_j * y_j * K(x_j, x_i) + bias
                    double fi = bias;
                    for (int j = 0; j < nSamples; ++j) {
                        if (alpha[j] == 0.0) continue;
                        double dot = 0.0;
                        for (int b = 0; b < numBands; ++b)
                            dot += samples[j].features[b] * samples[i].features[b];
                        fi += alpha[j] * y[j] * dot;
                    }

                    double Ei = fi - y[i];

                    // Check KKT conditions
                    if ((y[i] * Ei < -tol && alpha[i] < C) ||
                        (y[i] * Ei > tol && alpha[i] > 0)) {

                        // Pick random j != i
                        int j = i;
                        while (j == i)
                            j = std::uniform_int_distribution<int>(0, nSamples - 1)(rng);

                        // Compute f(x_j)
                        double fj = bias;
                        for (int k = 0; k < nSamples; ++k) {
                            if (alpha[k] == 0.0) continue;
                            double dot = 0.0;
                            for (int b = 0; b < numBands; ++b)
                                dot += samples[k].features[b] * samples[j].features[b];
                            fj += alpha[k] * y[k] * dot;
                        }
                        double Ej = fj - y[j];

                        double alphaIOld = alpha[i];
                        double alphaJOld = alpha[j];

                        // Compute bounds L, H
                        double L, H;
                        if (y[i] != y[j]) {
                            L = std::max(0.0, alpha[j] - alpha[i]);
                            H = std::min(C, C + alpha[j] - alpha[i]);
                        } else {
                            L = std::max(0.0, alpha[i] + alpha[j] - C);
                            H = std::min(C, alpha[i] + alpha[j]);
                        }

                        if (std::abs(L - H) < 1e-10) continue;

                        // Compute eta = 2 * K(i,j) - K(i,i) - K(j,j)
                        double kii = 0.0, kjj = 0.0, kij = 0.0;
                        for (int b = 0; b < numBands; ++b) {
                            kii += samples[i].features[b] * samples[i].features[b];
                            kjj += samples[j].features[b] * samples[j].features[b];
                            kij += samples[i].features[b] * samples[j].features[b];
                        }
                        double eta = 2.0 * kij - kii - kjj;
                        if (eta >= 0) continue;

                        // Update alpha[j]
                        alpha[j] = alphaJOld - y[j] * (Ei - Ej) / eta;
                        alpha[j] = std::min(H, std::max(L, alpha[j]));

                        if (std::abs(alpha[j] - alphaJOld) < 1e-5) continue;

                        // Update alpha[i]
                        alpha[i] = alphaIOld + y[i] * y[j] * (alphaJOld - alpha[j]);

                        // Update bias
                        double b1 = bias - Ei
                                    - y[i] * (alpha[i] - alphaIOld) * kii
                                    - y[j] * (alpha[j] - alphaJOld) * kij;
                        double b2 = bias - Ej
                                    - y[i] * (alpha[i] - alphaIOld) * kij
                                    - y[j] * (alpha[j] - alphaJOld) * kjj;

                        if (alpha[i] > 0 && alpha[i] < C)
                            bias = b1;
                        else if (alpha[j] > 0 && alpha[j] < C)
                            bias = b2;
                        else
                            bias = (b1 + b2) / 2.0;

                        numChanged++;
                    }
                }

                if (numChanged == 0) break;
            }

            // Extract weight vector: w = sum_i alpha_i * y_i * x_i
            SvmModel model;
            model.classId = targetClass;
            model.weights.resize(numBands, 0.0);
            model.bias = bias;

            for (int i = 0; i < nSamples; ++i) {
                if (alpha[i] == 0.0) continue;
                for (int b = 0; b < numBands; ++b)
                    model.weights[b] += alpha[i] * y[i] * samples[i].features[b];
            }

            models.push_back(std::move(model));

            reportProgress(0.1 + 0.6 * static_cast<double>(ci + 1) / numClasses,
                           QString("Trained SVM for class %1 of %2").arg(ci + 1).arg(numClasses));
        }

        // --------------------------------------------------------------------
        // 5. Classify all pixels
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

            double bestScore = -std::numeric_limits<double>::max();
            int bestClass = classIds[0];

            for (const auto& model : models) {
                double score = model.bias;
                for (int b = 0; b < numBands; ++b)
                    score += model.weights[b] * (*bands[b])[i];

                if (score > bestScore) {
                    bestScore = score;
                    bestClass = model.classId;
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

REGISTER_MODULE(SvmModule)

} // namespace aplaceholder
