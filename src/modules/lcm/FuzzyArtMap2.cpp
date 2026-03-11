#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>

namespace aplaceholder {

class FuzzyArtMap2Module : public Module {
public:
    QString name() const override { return "FUZZY_ARTMAP2"; }
    QString description() const override {
        return "Fuzzy ARTMAP neural network for land change transition potential modeling. "
               "Uses complement coding and vigilance-controlled category learning.";
    }
    QString category() const override { return "Land Change Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::raster("change_raster", "Change raster (binary 0/1)",
                "Binary raster where 1 = change, 0 = no change"),
            ParameterDef::file("drivers", "Driver variables (comma-separated)",
                "Raster layers of explanatory/driver variables"),
            ParameterDef::output("output", "Output transition potential map"),
            ParameterDef::real("vigilance", "Vigilance parameter", 0.75, 0.0, 1.0,
                "Controls category granularity (higher = more specific categories)"),
        };
    }

    bool execute() override {
        QString changePath = parameter("change_raster").toString();
        QStringList driverFiles = parameter("drivers").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : driverFiles) f = f.trimmed();
        int nDrivers = driverFiles.size();
        double vigilance = parameter("vigilance").toDouble();
        QString outPath = parameter("output").toString();

        if (nDrivers < 1) {
            setError("At least one driver variable is required");
            return false;
        }

        reportProgress(0.0, "Reading change raster...");

        auto changeRaster = GdalIO::read(changePath);
        if (!changeRaster) { setError("Failed to read change raster"); return false; }

        int cols = changeRaster->cols(), rows = changeRaster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = changeRaster->noDataValue();
        bool hasND = changeRaster->hasNoData();

        reportProgress(0.05, "Reading driver variables...");

        std::vector<std::unique_ptr<Raster>> drivers;
        drivers.reserve(nDrivers);
        for (int d = 0; d < nDrivers; ++d) {
            auto r = GdalIO::read(driverFiles[d]);
            if (!r) { setError(QString("Failed to read: %1").arg(driverFiles[d])); return false; }
            if (r->cols() != cols || r->rows() != rows) {
                setError("All rasters must have the same dimensions");
                return false;
            }
            drivers.push_back(std::move(r));
        }

        reportProgress(0.1, "Normalizing and extracting training data...");

        // Compute min/max for normalization
        std::vector<double> dMin(nDrivers, 1e30), dMax(nDrivers, -1e30);
        for (int64_t i = 0; i < total; ++i) {
            double changeVal = changeRaster->data(0)[i];
            if (hasND && changeVal == noData) continue;
            bool valid = true;
            for (int d = 0; d < nDrivers; ++d) {
                double val = drivers[d]->data(0)[i];
                if (hasND && val == noData) { valid = false; break; }
            }
            if (!valid) continue;
            for (int d = 0; d < nDrivers; ++d) {
                double val = drivers[d]->data(0)[i];
                dMin[d] = std::min(dMin[d], val);
                dMax[d] = std::max(dMax[d], val);
            }
        }

        for (int d = 0; d < nDrivers; ++d)
            if (dMax[d] - dMin[d] < 1e-15) dMax[d] = dMin[d] + 1.0;

        // Collect training samples with complement coding
        // Input dimension M = 2 * nDrivers (complement coded)
        int M = 2 * nDrivers;
        struct Sample {
            std::vector<double> input; // complement coded
            int label; // 0 or 1
        };
        std::vector<Sample> trainData;

        for (int64_t i = 0; i < total; ++i) {
            double changeVal = changeRaster->data(0)[i];
            if (hasND && changeVal == noData) continue;
            bool valid = true;
            std::vector<double> normalized(nDrivers);
            for (int d = 0; d < nDrivers; ++d) {
                double val = drivers[d]->data(0)[i];
                if (hasND && val == noData) { valid = false; break; }
                normalized[d] = (val - dMin[d]) / (dMax[d] - dMin[d]);
            }
            if (!valid) continue;

            Sample s;
            s.input.resize(M);
            for (int d = 0; d < nDrivers; ++d) {
                s.input[d] = normalized[d];
                s.input[nDrivers + d] = 1.0 - normalized[d]; // complement
            }
            s.label = (changeVal != 0.0) ? 1 : 0;
            trainData.push_back(std::move(s));
        }

        int64_t nSamples = trainData.size();
        if (nSamples < 10) {
            setError("Insufficient valid training samples");
            return false;
        }

        reportProgress(0.2, "Training Fuzzy ARTMAP network...");

        // Fuzzy ARTMAP parameters
        double alpha = 0.001; // choice parameter (small)
        double beta = 1.0;    // learning rate (fast learning)
        double eps = -0.001;  // match tracking epsilon

        // Category templates: weights for ARTa
        struct Category {
            std::vector<double> w; // weight vector (initialized to 1)
            int mapLabel;          // associated output class
        };
        std::vector<Category> categories;

        // Shuffle training order
        std::mt19937 rng(42);
        std::vector<int64_t> order(nSamples);
        for (int64_t i = 0; i < nSamples; ++i) order[i] = i;

        // Training epochs
        int numEpochs = 5;
        for (int epoch = 0; epoch < numEpochs; ++epoch) {
            std::shuffle(order.begin(), order.end(), rng);

            for (int64_t si = 0; si < nSamples; ++si) {
                const auto& sample = trainData[order[si]];
                const auto& I = sample.input;
                int target = sample.label;

                // Compute input norm |I|
                double normI = 0.0;
                for (int d = 0; d < M; ++d) normI += I[d];

                // Find winning category (choice function)
                double currentVigilance = vigilance;
                bool resonance = false;

                while (!resonance) {
                    int bestJ = -1;
                    double bestT = -1.0;

                    for (int j = 0; j < (int)categories.size(); ++j) {
                        // Fuzzy AND: min(I, w)
                        double fuzzyAndNorm = 0.0;
                        for (int d = 0; d < M; ++d)
                            fuzzyAndNorm += std::min(I[d], categories[j].w[d]);

                        // Category norm
                        double wNorm = 0.0;
                        for (int d = 0; d < M; ++d) wNorm += categories[j].w[d];

                        // Choice function T = |I ^ w| / (alpha + |w|)
                        double T = fuzzyAndNorm / (alpha + wNorm);

                        // Vigilance test: |I ^ w| / |I| >= rho
                        double matchVal = fuzzyAndNorm / (normI + 1e-15);

                        if (matchVal >= currentVigilance && T > bestT) {
                            // Check if map field matches
                            if (categories[j].mapLabel == target) {
                                bestT = T;
                                bestJ = j;
                            } else {
                                // Match tracking: raise vigilance
                                currentVigilance = matchVal + eps;
                                if (currentVigilance > 1.0) break;
                            }
                        }
                    }

                    if (bestJ >= 0) {
                        // Update weights: w_new = beta * (I ^ w) + (1-beta) * w
                        for (int d = 0; d < M; ++d) {
                            double fuzzyAnd = std::min(I[d], categories[bestJ].w[d]);
                            categories[bestJ].w[d] = beta * fuzzyAnd + (1.0 - beta) * categories[bestJ].w[d];
                        }
                        resonance = true;
                    } else {
                        // Create new category
                        Category newCat;
                        newCat.w = I; // fast learning: set to input
                        newCat.mapLabel = target;
                        categories.push_back(std::move(newCat));
                        resonance = true;
                    }
                }
            }

            reportProgress(0.2 + 0.4 * (epoch + 1.0) / numEpochs);
        }

        reportProgress(0.6, "Computing transition potential surface...");

        // Apply trained network to all pixels
        Raster potentialOut(cols, rows, 1, DataType::Float64);
        potentialOut.setGeoTransform(changeRaster->geoTransform());
        potentialOut.setProjection(changeRaster->projection());
        potentialOut.setNoDataValue(noData);
        auto& potData = potentialOut.data(0);

        for (int64_t i = 0; i < total; ++i) {
            bool valid = true;
            std::vector<double> I(M);
            for (int d = 0; d < nDrivers; ++d) {
                double val = drivers[d]->data(0)[i];
                if (hasND && val == noData) { valid = false; break; }
                double norm = (val - dMin[d]) / (dMax[d] - dMin[d]);
                norm = std::clamp(norm, 0.0, 1.0);
                I[d] = norm;
                I[nDrivers + d] = 1.0 - norm;
            }

            if (!valid) {
                potData[i] = noData;
                continue;
            }

            // Find best matching category and compute soft output
            // Use weighted vote across categories matching class 1
            double totalActivation = 0.0;
            double changeActivation = 0.0;
            double normI = 0.0;
            for (int d = 0; d < M; ++d) normI += I[d];

            for (int j = 0; j < (int)categories.size(); ++j) {
                double fuzzyAndNorm = 0.0;
                for (int d = 0; d < M; ++d)
                    fuzzyAndNorm += std::min(I[d], categories[j].w[d]);

                double wNorm = 0.0;
                for (int d = 0; d < M; ++d) wNorm += categories[j].w[d];

                double T = fuzzyAndNorm / (alpha + wNorm);
                double matchVal = fuzzyAndNorm / (normI + 1e-15);

                if (matchVal >= vigilance * 0.5) { // relaxed vigilance for prediction
                    totalActivation += T;
                    if (categories[j].mapLabel == 1)
                        changeActivation += T;
                }
            }

            potData[i] = (totalActivation > 1e-15) ? changeActivation / totalActivation : 0.0;

            if (i % 500000 == 0)
                reportProgress(0.6 + 0.35 * static_cast<double>(i) / total);
        }

        reportProgress(0.95, "Writing output...");
        if (!GdalIO::write(potentialOut, outPath)) {
            setError("Failed to write transition potential output");
            return false;
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(FuzzyArtMap2Module)

} // namespace aplaceholder
