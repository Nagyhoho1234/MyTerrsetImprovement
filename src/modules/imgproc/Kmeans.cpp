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

class KmeansModule : public Module {
public:
    QString name() const override { return "KMEANS"; }
    QString description() const override {
        return "K-Means clustering. Random initialization with iterative assignment and centroid "
               "update until convergence or maximum iterations reached.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("num_classes", "Number of clusters (K)", 10, 2, 256,
                "Number of clusters to partition the data into"),
            ParameterDef::integer("max_iterations", "Maximum iterations", 50, 1, 1000,
                "Maximum number of K-Means iterations"),
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

        int K = parameter("num_classes").toInt();
        int maxIter = parameter("max_iterations").toInt();

        // --------------------------------------------------------------------
        // 2. Collect band data and valid pixel list
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        std::vector<int64_t> validPixels;
        validPixels.reserve(total);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) continue;
            validPixels.push_back(i);
        }

        if (validPixels.empty()) {
            setError("No valid pixels in input bands");
            return false;
        }

        int64_t nValid = static_cast<int64_t>(validPixels.size());

        reportProgress(0.05, "Initializing K-Means centroids...");

        // --------------------------------------------------------------------
        // 3. Initialize centroids randomly from valid pixels
        // --------------------------------------------------------------------
        std::mt19937 rng(42);
        std::vector<std::vector<double>> centroids(K, std::vector<double>(numBands));

        {
            std::vector<int64_t> perm(nValid);
            std::iota(perm.begin(), perm.end(), 0);
            std::shuffle(perm.begin(), perm.end(), rng);
            for (int c = 0; c < K; ++c) {
                int64_t pi = validPixels[perm[c % nValid]];
                for (int b = 0; b < numBands; ++b)
                    centroids[c][b] = (*bands[b])[pi];
            }
        }

        std::vector<int> assignment(nValid, -1);

        // --------------------------------------------------------------------
        // 4. K-Means iterations
        // --------------------------------------------------------------------
        for (int iter = 0; iter < maxIter; ++iter) {
            reportProgress(0.05 + 0.85 * static_cast<double>(iter) / maxIter,
                           QString("K-Means iteration %1 of %2...").arg(iter + 1).arg(maxIter));

            // Assign pixels to nearest centroid
            bool changed = false;
            std::vector<std::vector<double>> sums(K, std::vector<double>(numBands, 0.0));
            std::vector<int64_t> counts(K, 0);

            for (int64_t vi = 0; vi < nValid; ++vi) {
                int64_t pi = validPixels[vi];
                double bestDist = std::numeric_limits<double>::max();
                int bestC = 0;

                for (int c = 0; c < K; ++c) {
                    double distSq = 0.0;
                    for (int b = 0; b < numBands; ++b) {
                        double diff = (*bands[b])[pi] - centroids[c][b];
                        distSq += diff * diff;
                    }
                    if (distSq < bestDist) {
                        bestDist = distSq;
                        bestC = c;
                    }
                }

                if (assignment[vi] != bestC) {
                    assignment[vi] = bestC;
                    changed = true;
                }

                for (int b = 0; b < numBands; ++b)
                    sums[bestC][b] += (*bands[b])[pi];
                counts[bestC]++;
            }

            // Update centroids
            for (int c = 0; c < K; ++c) {
                if (counts[c] > 0) {
                    for (int b = 0; b < numBands; ++b)
                        centroids[c][b] = sums[c][b] / counts[c];
                }
            }

            if (!changed) break;
        }

        // --------------------------------------------------------------------
        // 5. Create output
        // --------------------------------------------------------------------
        reportProgress(0.9, "Writing output...");

        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);
        std::fill(out.begin(), out.end(), 0.0);

        for (int64_t vi = 0; vi < nValid; ++vi) {
            int64_t pi = validPixels[vi];
            out[pi] = assignment[vi] + 1; // 1-based class ID
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(KmeansModule)

} // namespace aplaceholder
