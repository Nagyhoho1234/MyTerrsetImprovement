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

class RbfnnModule : public Module {
public:
    QString name() const override { return "RBFNN"; }
    QString description() const override {
        return "Radial Basis Function Neural Network. Hidden layer centers determined by K-Means "
               "clustering, Gaussian activation, and least-squares output weights.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster where pixel values indicate class IDs (0 = unclassified)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("num_centers", "Number of RBF centers", 50, 2, 10000,
                "Number of hidden layer neurons (K-Means cluster centers)"),
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

        int numCenters = parameter("num_centers").toInt();

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

        if (numCenters > nSamples) numCenters = nSamples;

        // Build class index map
        std::vector<int> classToIndex(classIds.back() + 1, -1);
        for (int ci = 0; ci < numClasses; ++ci)
            classToIndex[classIds[ci]] = ci;

        reportProgress(0.05, "Computing RBF centers via K-Means...");

        // --------------------------------------------------------------------
        // 4. K-Means on training data to find RBF centers
        // --------------------------------------------------------------------
        std::mt19937 rng(42);
        std::vector<std::vector<double>> centers(numCenters, std::vector<double>(numBands));

        // Initialize centers from random training samples
        {
            std::vector<int> perm(nSamples);
            std::iota(perm.begin(), perm.end(), 0);
            std::shuffle(perm.begin(), perm.end(), rng);
            for (int c = 0; c < numCenters; ++c)
                centers[c] = samples[perm[c % nSamples]].features;
        }

        std::vector<int> assignment(nSamples, 0);

        for (int iter = 0; iter < 50; ++iter) {
            // Assign samples to nearest center
            bool changed = false;
            for (int i = 0; i < nSamples; ++i) {
                double bestDist = std::numeric_limits<double>::max();
                int bestC = 0;
                for (int c = 0; c < numCenters; ++c) {
                    double distSq = 0.0;
                    for (int b = 0; b < numBands; ++b) {
                        double diff = samples[i].features[b] - centers[c][b];
                        distSq += diff * diff;
                    }
                    if (distSq < bestDist) {
                        bestDist = distSq;
                        bestC = c;
                    }
                }
                if (assignment[i] != bestC) {
                    assignment[i] = bestC;
                    changed = true;
                }
            }
            if (!changed) break;

            // Update centers
            std::vector<std::vector<double>> sums(numCenters, std::vector<double>(numBands, 0.0));
            std::vector<int> counts(numCenters, 0);
            for (int i = 0; i < nSamples; ++i) {
                int c = assignment[i];
                for (int b = 0; b < numBands; ++b)
                    sums[c][b] += samples[i].features[b];
                counts[c]++;
            }
            for (int c = 0; c < numCenters; ++c) {
                if (counts[c] > 0) {
                    for (int b = 0; b < numBands; ++b)
                        centers[c][b] = sums[c][b] / counts[c];
                }
            }
        }

        reportProgress(0.3, "Computing RBF widths and activation matrix...");

        // --------------------------------------------------------------------
        // 5. Compute RBF widths (sigma per center)
        // --------------------------------------------------------------------
        // Use average distance to the nearest other center as sigma
        std::vector<double> sigma(numCenters);
        for (int c = 0; c < numCenters; ++c) {
            double minDist = std::numeric_limits<double>::max();
            for (int c2 = 0; c2 < numCenters; ++c2) {
                if (c2 == c) continue;
                double distSq = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = centers[c][b] - centers[c2][b];
                    distSq += diff * diff;
                }
                if (distSq < minDist) minDist = distSq;
            }
            sigma[c] = std::sqrt(minDist);
            if (sigma[c] < 1e-10) sigma[c] = 1.0;
        }

        // --------------------------------------------------------------------
        // 6. Build activation matrix H[nSamples][numCenters] and solve
        //    least-squares for output weights W[numCenters][numClasses]
        // --------------------------------------------------------------------
        // H[i][c] = exp(-||x_i - center_c||^2 / (2 * sigma_c^2))
        // Solve: H * W = T (target matrix) using normal equations:
        //   W = (H^T H)^{-1} H^T T

        // Compute H^T H [numCenters x numCenters] and H^T T [numCenters x numClasses]
        std::vector<std::vector<double>> HtH(numCenters, std::vector<double>(numCenters, 0.0));
        std::vector<std::vector<double>> HtT(numCenters, std::vector<double>(numClasses, 0.0));

        // Process sample-by-sample to avoid storing full H matrix
        std::vector<double> hi(numCenters);
        for (int i = 0; i < nSamples; ++i) {
            // Compute activations for this sample
            for (int c = 0; c < numCenters; ++c) {
                double distSq = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = samples[i].features[b] - centers[c][b];
                    distSq += diff * diff;
                }
                hi[c] = std::exp(-distSq / (2.0 * sigma[c] * sigma[c]));
            }

            // Accumulate H^T H
            for (int c1 = 0; c1 < numCenters; ++c1)
                for (int c2 = c1; c2 < numCenters; ++c2) {
                    double val = hi[c1] * hi[c2];
                    HtH[c1][c2] += val;
                    if (c1 != c2) HtH[c2][c1] += val;
                }

            // Target: 1.0 for correct class column, 0.0 otherwise
            int targetIdx = classToIndex[samples[i].classId];
            for (int c = 0; c < numCenters; ++c)
                HtT[c][targetIdx] += hi[c];
        }

        // Regularize HtH
        for (int c = 0; c < numCenters; ++c)
            HtH[c][c] += 1e-6;

        reportProgress(0.5, "Solving for output weights (Cholesky)...");

        // Cholesky decomposition of HtH = L * L^T
        std::vector<std::vector<double>> L(numCenters, std::vector<double>(numCenters, 0.0));
        for (int i = 0; i < numCenters; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = HtH[i][j];
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

        // Solve L * y = HtT for each class column, then L^T * W = y
        std::vector<std::vector<double>> W(numCenters, std::vector<double>(numClasses, 0.0));

        for (int cl = 0; cl < numClasses; ++cl) {
            // Forward substitution: L * y = HtT[:,cl]
            std::vector<double> y(numCenters, 0.0);
            for (int i = 0; i < numCenters; ++i) {
                double sum = HtT[i][cl];
                for (int k = 0; k < i; ++k)
                    sum -= L[i][k] * y[k];
                y[i] = sum / L[i][i];
            }

            // Back substitution: L^T * w = y
            std::vector<double> w(numCenters, 0.0);
            for (int i = numCenters - 1; i >= 0; --i) {
                double sum = y[i];
                for (int k = i + 1; k < numCenters; ++k)
                    sum -= L[k][i] * w[k];
                w[i] = sum / L[i][i];
            }

            for (int c = 0; c < numCenters; ++c)
                W[c][cl] = w[c];
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

            // Compute RBF activations
            for (int c = 0; c < numCenters; ++c) {
                double distSq = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = (*bands[b])[i] - centers[c][b];
                    distSq += diff * diff;
                }
                hi[c] = std::exp(-distSq / (2.0 * sigma[c] * sigma[c]));
            }

            // Compute output for each class
            double bestVal = -std::numeric_limits<double>::max();
            int bestClass = classIds[0];
            for (int cl = 0; cl < numClasses; ++cl) {
                double val = 0.0;
                for (int c = 0; c < numCenters; ++c)
                    val += hi[c] * W[c][cl];
                if (val > bestVal) {
                    bestVal = val;
                    bestClass = classIds[cl];
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

REGISTER_MODULE(RbfnnModule)

} // namespace aplaceholder
