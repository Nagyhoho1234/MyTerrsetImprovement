#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <random>
#include <limits>
#include <numeric>

namespace aplaceholder {

class ClusterModule : public Module {
public:
    QString name() const override { return "CLUSTER"; }
    QString description() const override {
        return "Unsupervised clustering of multi-band imagery. Groups pixels into "
               "clusters using iterative partitioning algorithms.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input multi-band image"),
            ParameterDef::output("output", "Output cluster image"),
            ParameterDef::integer("num_clusters", "Number of clusters", 10, 2, 256,
                "Target number of clusters to generate"),
            ParameterDef::integer("max_iterations", "Maximum iterations", 50, 1, 1000,
                "Maximum number of clustering iterations"),
            ParameterDef::combo("method", "Clustering method",
                {"K-Means", "ISODATA"}, 0,
                "Algorithm to use for unsupervised clustering"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input image");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int numBands = raster->bands();
        int K = parameter("num_clusters").toInt();
        int maxIter = parameter("max_iterations").toInt();
        int method = parameter("method").toInt();

        if (method == 1) {
            // ISODATA: not yet implemented, fall back to K-Means with a note
            reportProgress(0.0, "ISODATA not yet implemented; using K-Means...");
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        // Collect band data pointers
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &raster->data(b);

        // Find min/max per band for initialization (excluding nodata)
        std::vector<double> bandMin(numBands, std::numeric_limits<double>::max());
        std::vector<double> bandMax(numBands, std::numeric_limits<double>::lowest());
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) continue;
            for (int b = 0; b < numBands; ++b) {
                double v = (*bands[b])[i];
                if (v < bandMin[b]) bandMin[b] = v;
                if (v > bandMax[b]) bandMax[b] = v;
            }
        }

        // Initialize K cluster centers uniformly spread across data range
        // Use deterministic seeding for reproducibility
        std::mt19937 rng(42);
        std::vector<std::vector<double>> centers(K, std::vector<double>(numBands));
        for (int k = 0; k < K; ++k) {
            for (int b = 0; b < numBands; ++b) {
                // Spread centers evenly plus small jitter
                double frac = (K > 1) ? static_cast<double>(k) / (K - 1) : 0.5;
                double range = bandMax[b] - bandMin[b];
                std::uniform_real_distribution<double> jitter(-range * 0.02, range * 0.02);
                centers[k][b] = bandMin[b] + frac * range + jitter(rng);
            }
        }

        // Pixel-to-cluster assignment
        std::vector<int> assignment(total, -1);

        // K-Means iterations
        for (int iter = 0; iter < maxIter; ++iter) {
            reportProgress(static_cast<double>(iter) / maxIter,
                           QString("Iteration %1 of %2").arg(iter + 1).arg(maxIter));

            bool changed = false;

            // Step 1: Assign each pixel to nearest center
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && (*bands[0])[i] == noData) {
                    assignment[i] = -1;
                    continue;
                }

                double bestDist = std::numeric_limits<double>::max();
                int bestK = 0;
                for (int k = 0; k < K; ++k) {
                    double dist = 0.0;
                    for (int b = 0; b < numBands; ++b) {
                        double d = (*bands[b])[i] - centers[k][b];
                        dist += d * d;
                    }
                    if (dist < bestDist) {
                        bestDist = dist;
                        bestK = k;
                    }
                }

                if (assignment[i] != bestK) {
                    assignment[i] = bestK;
                    changed = true;
                }
            }

            if (!changed) break; // Converged

            // Step 2: Recompute centers as mean of assigned pixels
            std::vector<std::vector<double>> newCenters(K, std::vector<double>(numBands, 0.0));
            std::vector<int64_t> counts(K, 0);

            for (int64_t i = 0; i < total; ++i) {
                int k = assignment[i];
                if (k < 0) continue;
                counts[k]++;
                for (int b = 0; b < numBands; ++b)
                    newCenters[k][b] += (*bands[b])[i];
            }

            for (int k = 0; k < K; ++k) {
                if (counts[k] > 0) {
                    for (int b = 0; b < numBands; ++b)
                        centers[k][b] = newCenters[k][b] / counts[k];
                }
                // If a cluster is empty, keep the old center (or re-seed)
            }
        }

        // Write output: cluster labels as integer raster (1-based class IDs)
        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(-9999);

        auto& out = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            if (assignment[i] < 0)
                out[i] = -9999;
            else
                out[i] = assignment[i] + 1; // 1-based cluster IDs
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ClusterModule)

} // namespace aplaceholder
