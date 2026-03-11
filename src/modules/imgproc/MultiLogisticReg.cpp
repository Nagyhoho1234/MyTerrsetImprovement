#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class MultiLogisticRegModule : public Module {
public:
    QString name() const override { return "MULTILOGISTICREG"; }
    QString description() const override {
        return "Multinomial logistic regression classifier. One-vs-rest logistic regression "
               "using Iteratively Reweighted Least Squares (IRLS) for parameter estimation.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster where pixel values indicate class IDs (0 = unclassified)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("max_iterations", "Maximum IRLS iterations", 50, 1, 1000,
                "Maximum number of IRLS iterations per class"),
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

        int maxIter = parameter("max_iterations").toInt();

        // --------------------------------------------------------------------
        // 3. Collect training samples
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        const auto& trainData = trainRaster->data(0);

        // Features include intercept: [1, x1, x2, ..., xB]
        int nFeatures = numBands + 1; // +1 for intercept

        struct Sample {
            std::vector<double> features; // length = nFeatures (with intercept)
            int classId;
        };

        std::vector<Sample> samples;
        std::vector<int> classIds;

        for (int64_t i = 0; i < total; ++i) {
            int cls = static_cast<int>(trainData[i]);
            if (cls <= 0) continue;
            if (hasND && (*bands[0])[i] == noData) continue;

            Sample s;
            s.features.resize(nFeatures);
            s.features[0] = 1.0; // intercept
            for (int b = 0; b < numBands; ++b)
                s.features[b + 1] = (*bands[b])[i];
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

        reportProgress(0.05, "Training logistic regression models...");

        // --------------------------------------------------------------------
        // 4. Train one-vs-rest logistic regression using IRLS
        // --------------------------------------------------------------------
        // For each class: binary logistic regression with IRLS
        // p(y=1|x) = sigmoid(beta^T * x)
        // IRLS update: beta_new = (X^T W X)^{-1} X^T W z
        // where W = diag(p * (1 - p)), z = X*beta + W^{-1}(y - p)

        struct LogRegModel {
            int classId;
            std::vector<double> beta; // nFeatures
        };

        std::vector<LogRegModel> models;

        auto sigmoid = [](double x) -> double {
            return 1.0 / (1.0 + std::exp(-std::clamp(x, -500.0, 500.0)));
        };

        for (int ci = 0; ci < numClasses; ++ci) {
            int targetClass = classIds[ci];

            // Binary labels
            std::vector<double> y(nSamples);
            for (int i = 0; i < nSamples; ++i)
                y[i] = (samples[i].classId == targetClass) ? 1.0 : 0.0;

            // Initialize beta to zeros
            std::vector<double> beta(nFeatures, 0.0);

            std::vector<double> p(nSamples);
            std::vector<double> w(nSamples);

            for (int iter = 0; iter < maxIter; ++iter) {
                // Compute predictions p = sigmoid(X * beta)
                for (int i = 0; i < nSamples; ++i) {
                    double z = 0.0;
                    for (int f = 0; f < nFeatures; ++f)
                        z += samples[i].features[f] * beta[f];
                    p[i] = sigmoid(z);
                    w[i] = p[i] * (1.0 - p[i]);
                    if (w[i] < 1e-10) w[i] = 1e-10; // prevent singularity
                }

                // Compute X^T W X [nFeatures x nFeatures]
                std::vector<std::vector<double>> XtWX(nFeatures, std::vector<double>(nFeatures, 0.0));
                for (int i = 0; i < nSamples; ++i) {
                    for (int f1 = 0; f1 < nFeatures; ++f1)
                        for (int f2 = f1; f2 < nFeatures; ++f2) {
                            double val = samples[i].features[f1] * w[i] * samples[i].features[f2];
                            XtWX[f1][f2] += val;
                            if (f1 != f2) XtWX[f2][f1] += val;
                        }
                }

                // Regularize
                for (int f = 0; f < nFeatures; ++f)
                    XtWX[f][f] += 1e-6;

                // Compute X^T (W * X * beta + (y - p)) = X^T W z
                // where z_i = X_i * beta + (y_i - p_i) / w_i
                std::vector<double> XtWz(nFeatures, 0.0);
                for (int i = 0; i < nSamples; ++i) {
                    double zi = 0.0;
                    for (int f = 0; f < nFeatures; ++f)
                        zi += samples[i].features[f] * beta[f];
                    zi += (y[i] - p[i]) / w[i];

                    for (int f = 0; f < nFeatures; ++f)
                        XtWz[f] += samples[i].features[f] * w[i] * zi;
                }

                // Solve XtWX * beta_new = XtWz using Cholesky
                std::vector<std::vector<double>> L(nFeatures, std::vector<double>(nFeatures, 0.0));
                for (int i = 0; i < nFeatures; ++i) {
                    for (int j = 0; j <= i; ++j) {
                        double sum = XtWX[i][j];
                        for (int k = 0; k < j; ++k)
                            sum -= L[i][k] * L[j][k];
                        if (i == j) {
                            if (sum <= 0) sum = 1e-10;
                            L[i][j] = std::sqrt(sum);
                        } else {
                            L[i][j] = sum / L[j][j];
                        }
                    }
                }

                // Forward substitution: L * yy = XtWz
                std::vector<double> yy(nFeatures);
                for (int i = 0; i < nFeatures; ++i) {
                    double sum = XtWz[i];
                    for (int k = 0; k < i; ++k)
                        sum -= L[i][k] * yy[k];
                    yy[i] = sum / L[i][i];
                }

                // Back substitution: L^T * beta_new = yy
                std::vector<double> betaNew(nFeatures);
                for (int i = nFeatures - 1; i >= 0; --i) {
                    double sum = yy[i];
                    for (int k = i + 1; k < nFeatures; ++k)
                        sum -= L[k][i] * betaNew[k];
                    betaNew[i] = sum / L[i][i];
                }

                // Check convergence
                double maxDiff = 0.0;
                for (int f = 0; f < nFeatures; ++f)
                    maxDiff = std::max(maxDiff, std::abs(betaNew[f] - beta[f]));

                beta = std::move(betaNew);

                if (maxDiff < 1e-6) break;
            }

            LogRegModel model;
            model.classId = targetClass;
            model.beta = std::move(beta);
            models.push_back(std::move(model));

            reportProgress(0.05 + 0.65 * static_cast<double>(ci + 1) / numClasses,
                           QString("Trained model for class %1 of %2").arg(ci + 1).arg(numClasses));
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

            // Build feature vector with intercept
            std::vector<double> features(nFeatures);
            features[0] = 1.0;
            for (int b = 0; b < numBands; ++b)
                features[b + 1] = (*bands[b])[i];

            double bestProb = -1.0;
            int bestClass = classIds[0];

            for (const auto& model : models) {
                double z = 0.0;
                for (int f = 0; f < nFeatures; ++f)
                    z += model.beta[f] * features[f];
                double prob = sigmoid(z);

                if (prob > bestProb) {
                    bestProb = prob;
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

REGISTER_MODULE(MultiLogisticRegModule)

} // namespace aplaceholder
