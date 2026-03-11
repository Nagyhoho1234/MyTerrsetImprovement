#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>

namespace aplaceholder {

class IsoClustModule : public Module {
public:
    QString name() const override { return "ISOCLUST"; }
    QString description() const override {
        return "ISODATA unsupervised clustering. Iteratively assigns pixels to clusters, "
               "recomputes means, and splits/merges clusters based on variance and distance criteria.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("num_classes", "Number of clusters", 10, 2, 256,
                "Initial/target number of clusters"),
            ParameterDef::integer("max_iterations", "Maximum iterations", 20, 1, 200,
                "Maximum number of ISODATA iterations"),
            ParameterDef::real("merge_distance", "Merge distance threshold", 5.0, 0.0, 999999.0,
                "Clusters closer than this Euclidean distance will be merged"),
            ParameterDef::real("split_stddev", "Split standard deviation threshold", 10.0, 0.0, 999999.0,
                "Clusters with stddev above this in any band will be split"),
            ParameterDef::integer("min_members", "Minimum cluster members", 10, 1, 999999,
                "Clusters with fewer members are dissolved"),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Read input bands
        // ------------------------------------------------------------------
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
        double mergeDist = parameter("merge_distance").toDouble();
        double splitStddev = parameter("split_stddev").toDouble();
        int minMembers = parameter("min_members").toInt();

        // ------------------------------------------------------------------
        // 2. Collect band data and build valid pixel list
        // ------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // Identify valid pixels
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

        reportProgress(0.05, "Initializing clusters...");

        // ------------------------------------------------------------------
        // 3. Initialize cluster centroids by sampling pixels evenly
        // ------------------------------------------------------------------
        struct Cluster {
            std::vector<double> mean;   // per band
            std::vector<double> sumVal; // accumulator per band
            std::vector<double> sumSq;  // sum of squares per band
            int64_t count;
        };

        std::vector<Cluster> clusters(K);
        for (int c = 0; c < K; ++c) {
            int64_t idx = validPixels[(int64_t)c * nValid / K];
            clusters[c].mean.resize(numBands);
            clusters[c].sumVal.resize(numBands, 0.0);
            clusters[c].sumSq.resize(numBands, 0.0);
            clusters[c].count = 0;
            for (int b = 0; b < numBands; ++b)
                clusters[c].mean[b] = (*bands[b])[idx];
        }

        // Assignment array: which cluster each valid pixel belongs to
        std::vector<int> assignment(nValid, -1);

        // ------------------------------------------------------------------
        // 4. ISODATA iterations
        // ------------------------------------------------------------------
        for (int iter = 0; iter < maxIter; ++iter) {
            reportProgress(0.05 + 0.85 * static_cast<double>(iter) / maxIter,
                           QString("Iteration %1 of %2...").arg(iter + 1).arg(maxIter));

            int numClusters = static_cast<int>(clusters.size());
            if (numClusters == 0) {
                setError("All clusters dissolved; cannot continue");
                return false;
            }

            // Reset accumulators
            for (auto& cl : clusters) {
                std::fill(cl.sumVal.begin(), cl.sumVal.end(), 0.0);
                std::fill(cl.sumSq.begin(), cl.sumSq.end(), 0.0);
                cl.count = 0;
            }

            // 4a. Assign pixels to nearest cluster
            bool changed = false;
            for (int64_t vi = 0; vi < nValid; ++vi) {
                int64_t pi = validPixels[vi];
                double bestDist = std::numeric_limits<double>::max();
                int bestCluster = 0;

                for (int c = 0; c < numClusters; ++c) {
                    double distSq = 0.0;
                    for (int b = 0; b < numBands; ++b) {
                        double diff = (*bands[b])[pi] - clusters[c].mean[b];
                        distSq += diff * diff;
                    }
                    if (distSq < bestDist) {
                        bestDist = distSq;
                        bestCluster = c;
                    }
                }

                if (assignment[vi] != bestCluster) {
                    assignment[vi] = bestCluster;
                    changed = true;
                }

                // Accumulate statistics
                for (int b = 0; b < numBands; ++b) {
                    double val = (*bands[b])[pi];
                    clusters[bestCluster].sumVal[b] += val;
                    clusters[bestCluster].sumSq[b] += val * val;
                }
                clusters[bestCluster].count++;
            }

            // 4b. Recompute means
            for (auto& cl : clusters) {
                if (cl.count > 0) {
                    for (int b = 0; b < numBands; ++b)
                        cl.mean[b] = cl.sumVal[b] / cl.count;
                }
            }

            // 4c. Remove clusters with too few members
            {
                std::vector<Cluster> kept;
                std::vector<int> oldToNew(numClusters, -1);
                for (int c = 0; c < numClusters; ++c) {
                    if (clusters[c].count >= minMembers) {
                        oldToNew[c] = static_cast<int>(kept.size());
                        kept.push_back(std::move(clusters[c]));
                    }
                }
                if (!kept.empty()) {
                    clusters = std::move(kept);
                    // Remap assignments (dissolved pixels will be reassigned next iteration)
                    for (int64_t vi = 0; vi < nValid; ++vi) {
                        int old = assignment[vi];
                        if (old >= 0 && old < numClusters && oldToNew[old] >= 0)
                            assignment[vi] = oldToNew[old];
                        else
                            assignment[vi] = 0; // fallback, will be reassigned
                    }
                }
                numClusters = static_cast<int>(clusters.size());
            }

            // 4d. Merge close clusters
            {
                bool merged = true;
                while (merged) {
                    merged = false;
                    int nc = static_cast<int>(clusters.size());
                    for (int a = 0; a < nc && !merged; ++a) {
                        for (int b2 = a + 1; b2 < nc && !merged; ++b2) {
                            double distSq = 0.0;
                            for (int b = 0; b < numBands; ++b) {
                                double diff = clusters[a].mean[b] - clusters[b2].mean[b];
                                distSq += diff * diff;
                            }
                            if (std::sqrt(distSq) < mergeDist) {
                                // Merge b2 into a (weighted average)
                                int64_t totalCount = clusters[a].count + clusters[b2].count;
                                if (totalCount > 0) {
                                    for (int b = 0; b < numBands; ++b) {
                                        clusters[a].mean[b] =
                                            (clusters[a].mean[b] * clusters[a].count +
                                             clusters[b2].mean[b] * clusters[b2].count) / totalCount;
                                        clusters[a].sumVal[b] += clusters[b2].sumVal[b];
                                        clusters[a].sumSq[b] += clusters[b2].sumSq[b];
                                    }
                                    clusters[a].count = totalCount;
                                }
                                // Remap assignments from b2 to a
                                for (int64_t vi = 0; vi < nValid; ++vi) {
                                    if (assignment[vi] == b2)
                                        assignment[vi] = a;
                                    else if (assignment[vi] > b2)
                                        assignment[vi]--;
                                }
                                clusters.erase(clusters.begin() + b2);
                                merged = true;
                            }
                        }
                    }
                }
            }

            // 4e. Split clusters with high variance
            {
                int nc = static_cast<int>(clusters.size());
                std::vector<Cluster> toAdd;
                for (int c = 0; c < nc; ++c) {
                    if (clusters[c].count < 2 * minMembers) continue;

                    for (int b = 0; b < numBands; ++b) {
                        double variance = (clusters[c].sumSq[b] / clusters[c].count)
                                          - (clusters[c].mean[b] * clusters[c].mean[b]);
                        if (variance < 0) variance = 0;
                        double stddev = std::sqrt(variance);

                        if (stddev > splitStddev) {
                            // Split this cluster along band b
                            Cluster newCluster;
                            newCluster.mean.resize(numBands);
                            newCluster.sumVal.resize(numBands, 0.0);
                            newCluster.sumSq.resize(numBands, 0.0);
                            newCluster.count = 0;

                            for (int bb = 0; bb < numBands; ++bb)
                                newCluster.mean[bb] = clusters[c].mean[bb];

                            // Offset the two clusters slightly along the split band
                            clusters[c].mean[b] -= stddev * 0.5;
                            newCluster.mean[b] += stddev * 0.5;

                            toAdd.push_back(std::move(newCluster));
                            break; // only split once per cluster per iteration
                        }
                    }
                }
                for (auto& cl : toAdd)
                    clusters.push_back(std::move(cl));
            }

            // Early termination if no assignments changed
            if (!changed && iter > 0) break;
        }

        reportProgress(0.9, "Final assignment...");

        // ------------------------------------------------------------------
        // 5. Final assignment pass
        // ------------------------------------------------------------------
        int numClusters = static_cast<int>(clusters.size());
        for (int64_t vi = 0; vi < nValid; ++vi) {
            int64_t pi = validPixels[vi];
            double bestDist = std::numeric_limits<double>::max();
            int bestCluster = 0;

            for (int c = 0; c < numClusters; ++c) {
                double distSq = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = (*bands[b])[pi] - clusters[c].mean[b];
                    distSq += diff * diff;
                }
                if (distSq < bestDist) {
                    bestDist = distSq;
                    bestCluster = c;
                }
            }
            assignment[vi] = bestCluster;
        }

        // ------------------------------------------------------------------
        // 6. Create output
        // ------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);
        // Initialize to unclassified
        std::fill(out.begin(), out.end(), 0.0);

        // Cluster IDs are 1-based in output
        for (int64_t vi = 0; vi < nValid; ++vi) {
            int64_t pi = validPixels[vi];
            out[pi] = assignment[vi] + 1; // 1-based class ID
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(IsoClustModule)

} // namespace aplaceholder
