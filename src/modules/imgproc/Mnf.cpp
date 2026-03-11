#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class MnfModule : public Module {
public:
    QString name() const override { return "MNF"; }
    QString description() const override {
        return "Minimum Noise Fraction (MNF) transform. Similar to PCA but operates on "
               "noise-whitened data for dimensionality reduction of noisy hyperspectral "
               "imagery. Estimates the noise covariance from spatial differencing, whitens "
               "the noise, then applies PCA to the whitened data.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::output("output_prefix", "Output filename prefix",
                "Output rasters will be named prefix_mnf1, prefix_mnf2, etc."),
            ParameterDef::integer("num_components", "Number of MNF components", 3, 1, 999,
                "Number of output MNF components to retain"),
        };
    }

    bool execute() override {
        QString bandsParam = parameter("bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths)
            p = p.trimmed();

        QString prefix = parameter("output_prefix").toString();
        int numComp = parameter("num_components").toInt();
        int numBands = bandPaths.size();

        if (numBands < 2) {
            setError("At least 2 input bands are required.");
            return false;
        }
        if (numComp > numBands)
            numComp = numBands;

        // Read all bands
        reportProgress(0.0, "Reading input bands...");
        std::vector<std::unique_ptr<Raster>> rasters(numBands);
        int cols = 0, rows = 0;

        for (int b = 0; b < numBands; ++b) {
            rasters[b] = GdalIO::read(bandPaths[b]);
            if (!rasters[b]) {
                setError("Failed to read band: " + bandPaths[b]);
                return false;
            }
            if (b == 0) {
                cols = rasters[b]->cols();
                rows = rasters[b]->rows();
            } else if (rasters[b]->cols() != cols || rasters[b]->rows() != rows) {
                setError("All bands must have the same dimensions.");
                return false;
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = rasters[0]->hasNoData();
        double noData = rasters[0]->noDataValue();
        double outNoData = -9999.0;

        // Collect data pointers
        std::vector<const std::vector<double>*> data(numBands);
        for (int b = 0; b < numBands; ++b)
            data[b] = &rasters[b]->data(0);

        // Step 1: Estimate noise covariance using spatial differencing
        reportProgress(0.1, "Estimating noise covariance...");
        std::vector<double> noiseCov(numBands * numBands, 0.0);
        int64_t noiseSamples = 0;

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols - 1; ++c) {
                int64_t i1 = static_cast<int64_t>(r) * cols + c;
                int64_t i2 = i1 + 1;

                bool valid = true;
                if (hasND) {
                    for (int b = 0; b < numBands && valid; ++b) {
                        if ((*data[b])[i1] == noData || (*data[b])[i2] == noData)
                            valid = false;
                    }
                }
                if (!valid) continue;

                for (int bi = 0; bi < numBands; ++bi) {
                    double di = (*data[bi])[i1] - (*data[bi])[i2];
                    for (int bj = bi; bj < numBands; ++bj) {
                        double dj = (*data[bj])[i1] - (*data[bj])[i2];
                        noiseCov[bi * numBands + bj] += di * dj;
                    }
                }
                ++noiseSamples;
            }
        }

        if (noiseSamples < numBands) {
            setError("Not enough valid pixel pairs to estimate noise covariance.");
            return false;
        }

        for (int bi = 0; bi < numBands; ++bi) {
            for (int bj = bi; bj < numBands; ++bj) {
                noiseCov[bi * numBands + bj] /= (2.0 * noiseSamples);
                noiseCov[bj * numBands + bi] = noiseCov[bi * numBands + bj];
            }
        }

        // Step 2: Compute data covariance
        reportProgress(0.3, "Computing data covariance...");
        std::vector<double> means(numBands, 0.0);
        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            bool valid = true;
            if (hasND) {
                for (int b = 0; b < numBands && valid; ++b)
                    if ((*data[b])[i] == noData) valid = false;
            }
            if (!valid) continue;
            for (int b = 0; b < numBands; ++b)
                means[b] += (*data[b])[i];
            ++validCount;
        }

        if (validCount < numBands) {
            setError("Not enough valid pixels for MNF computation.");
            return false;
        }

        for (int b = 0; b < numBands; ++b)
            means[b] /= validCount;

        std::vector<double> dataCov(numBands * numBands, 0.0);
        for (int64_t i = 0; i < total; ++i) {
            bool valid = true;
            if (hasND) {
                for (int b = 0; b < numBands && valid; ++b)
                    if ((*data[b])[i] == noData) valid = false;
            }
            if (!valid) continue;
            for (int bi = 0; bi < numBands; ++bi) {
                double di = (*data[bi])[i] - means[bi];
                for (int bj = bi; bj < numBands; ++bj) {
                    dataCov[bi * numBands + bj] += di * ((*data[bj])[i] - means[bj]);
                }
            }
        }

        for (int bi = 0; bi < numBands; ++bi) {
            for (int bj = bi; bj < numBands; ++bj) {
                dataCov[bi * numBands + bj] /= (validCount - 1);
                dataCov[bj * numBands + bi] = dataCov[bi * numBands + bj];
            }
        }

        // Step 3: Power iteration to find top MNF components
        // MNF components are eigenvectors of Sigma_noise^{-1} * Sigma_data
        // We use simplified iterative approach: solve via noise-whitened PCA
        reportProgress(0.5, "Computing MNF transform...");

        // Add small regularization to noise covariance diagonal
        for (int b = 0; b < numBands; ++b)
            noiseCov[b * numBands + b] += 1e-10;

        // Cholesky-like noise whitening: approximate with diagonal sqrt
        // For simplicity, use diagonal elements of noise covariance for whitening
        std::vector<double> noiseStd(numBands);
        for (int b = 0; b < numBands; ++b)
            noiseStd[b] = std::sqrt(noiseCov[b * numBands + b]);

        // Whiten the data covariance
        std::vector<double> whiteCov(numBands * numBands);
        for (int bi = 0; bi < numBands; ++bi) {
            for (int bj = 0; bj < numBands; ++bj) {
                whiteCov[bi * numBands + bj] = dataCov[bi * numBands + bj] /
                    (noiseStd[bi] * noiseStd[bj]);
            }
        }

        // Power iteration for eigenvectors of whitened covariance
        std::vector<std::vector<double>> eigVecs(numComp, std::vector<double>(numBands));
        std::vector<double> eigVals(numComp, 0.0);

        std::vector<double> workCov = whiteCov;

        for (int comp = 0; comp < numComp; ++comp) {
            // Initialize with random-ish vector
            std::vector<double> vec(numBands);
            for (int b = 0; b < numBands; ++b)
                vec[b] = 1.0 / std::sqrt(numBands) + 0.01 * b;

            // Power iteration
            for (int iter = 0; iter < 200; ++iter) {
                std::vector<double> newVec(numBands, 0.0);
                for (int i = 0; i < numBands; ++i)
                    for (int j = 0; j < numBands; ++j)
                        newVec[i] += workCov[i * numBands + j] * vec[j];

                double norm = 0.0;
                for (int b = 0; b < numBands; ++b)
                    norm += newVec[b] * newVec[b];
                norm = std::sqrt(norm);

                if (norm < 1e-15) break;

                for (int b = 0; b < numBands; ++b)
                    vec[b] = newVec[b] / norm;

                eigVals[comp] = norm;
            }

            // Apply noise whitening to the eigenvector
            for (int b = 0; b < numBands; ++b)
                eigVecs[comp][b] = vec[b] / noiseStd[b];

            // Deflate
            for (int i = 0; i < numBands; ++i)
                for (int j = 0; j < numBands; ++j)
                    workCov[i * numBands + j] -= eigVals[comp] * vec[i] * vec[j];
        }

        // Step 4: Project data onto MNF components and write outputs
        for (int comp = 0; comp < numComp; ++comp) {
            reportProgress(0.6 + 0.35 * static_cast<double>(comp) / numComp,
                           QString("Computing MNF component %1...").arg(comp + 1));

            Raster output(cols, rows, 1, DataType::Float32);
            output.setGeoTransform(rasters[0]->geoTransform());
            output.setProjection(rasters[0]->projection());
            output.setNoDataValue(outNoData);
            auto& out = output.data(0);

            for (int64_t i = 0; i < total; ++i) {
                bool valid = true;
                if (hasND) {
                    for (int b = 0; b < numBands && valid; ++b)
                        if ((*data[b])[i] == noData) valid = false;
                }
                if (!valid) {
                    out[i] = outNoData;
                    continue;
                }

                double result = 0.0;
                for (int b = 0; b < numBands; ++b)
                    result += eigVecs[comp][b] * ((*data[b])[i] - means[b]);
                out[i] = result;
            }

            QString outPath = prefix + "_mnf" + QString::number(comp + 1);
            if (!GdalIO::write(output, outPath)) {
                setError("Failed to write MNF component: " + outPath);
                return false;
            }
        }

        reportProgress(1.0, "MNF transform complete.");
        return true;
    }
};

REGISTER_MODULE(MnfModule)

} // namespace aplaceholder
