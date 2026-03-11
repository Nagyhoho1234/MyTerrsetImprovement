#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>

namespace aplaceholder {

class ChainClusterModule : public Module {
public:
    QString name() const override { return "CHAINCLUSTER"; }
    QString description() const override {
        return "Chain clustering classifier. Groups pixels into clusters using a sequential "
               "chain-based agglomerative approach: each pixel is assigned to the nearest "
               "existing cluster or starts a new cluster if beyond the distance threshold.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input multi-band image"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::real("distance_threshold", "Distance threshold", 10.0, 0.001, 100000.0,
                "Maximum spectral distance to join an existing cluster"),
            ParameterDef::integer("max_clusters", "Maximum number of clusters", 50, 2, 10000,
                "Upper limit on number of clusters to create"),
            ParameterDef::integer("min_cluster_size", "Minimum cluster size (pixels)", 10, 1, 1000000,
                "Clusters smaller than this are merged into the nearest neighbor"),
            ParameterDef::boolean("merge_small", "Merge small clusters", true,
                "Merge clusters with fewer pixels than minimum into nearest cluster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input image");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int nBands = raster->bands();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double distThresh = parameter("distance_threshold").toDouble();
        int maxClusters = parameter("max_clusters").toInt();
        int minSize = parameter("min_cluster_size").toInt();
        bool mergeSmall = parameter("merge_small").toBool();

        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        std::vector<const std::vector<double>*> bands(nBands);
        for (int b = 0; b < nBands; ++b)
            bands[b] = &raster->data(b);

        reportProgress(0.05, "Running chain clustering...");

        // Cluster centroids: centroids[k][b] = mean of band b for cluster k
        std::vector<std::vector<double>> centroids;
        std::vector<int64_t> clusterCounts;
        std::vector<int> labels(total, -1);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                labels[i] = -1;
                continue;
            }

            std::vector<double> pixel(nBands);
            for (int b = 0; b < nBands; ++b)
                pixel[b] = (*bands[b])[i];

            double minDist = std::numeric_limits<double>::max();
            int nearestCluster = -1;
            int nClusters = static_cast<int>(centroids.size());

            for (int k = 0; k < nClusters; ++k) {
                double dist = 0.0;
                for (int b = 0; b < nBands; ++b) {
                    double d = pixel[b] - centroids[k][b];
                    dist += d * d;
                }
                dist = std::sqrt(dist);
                if (dist < minDist) {
                    minDist = dist;
                    nearestCluster = k;
                }
            }

            if (nearestCluster >= 0 && minDist <= distThresh) {
                labels[i] = nearestCluster;
                int64_t count = clusterCounts[nearestCluster] + 1;
                for (int b = 0; b < nBands; ++b) {
                    centroids[nearestCluster][b] =
                        centroids[nearestCluster][b] * (count - 1.0) / count +
                        pixel[b] / count;
                }
                clusterCounts[nearestCluster] = count;
            } else if (nClusters < maxClusters) {
                centroids.push_back(pixel);
                clusterCounts.push_back(1);
                labels[i] = nClusters;
            } else {
                labels[i] = nearestCluster;
                if (nearestCluster >= 0) {
                    int64_t count = clusterCounts[nearestCluster] + 1;
                    for (int b = 0; b < nBands; ++b) {
                        centroids[nearestCluster][b] =
                            centroids[nearestCluster][b] * (count - 1.0) / count +
                            pixel[b] / count;
                    }
                    clusterCounts[nearestCluster] = count;
                }
            }

            if (i % 500000 == 0)
                reportProgress(0.05 + 0.70 * static_cast<double>(i) / total);
        }

        reportProgress(0.80, "Merging small clusters...");

        if (mergeSmall) {
            int nClusters = static_cast<int>(centroids.size());
            std::vector<int> remap(nClusters);
            std::iota(remap.begin(), remap.end(), 0);

            for (int k = 0; k < nClusters; ++k) {
                if (clusterCounts[k] < minSize && clusterCounts[k] > 0) {
                    double minDist = std::numeric_limits<double>::max();
                    int nearest = -1;
                    for (int j = 0; j < nClusters; ++j) {
                        if (j == k || clusterCounts[j] < minSize) continue;
                        double dist = 0.0;
                        for (int b = 0; b < nBands; ++b) {
                            double d = centroids[k][b] - centroids[j][b];
                            dist += d * d;
                        }
                        if (dist < minDist) {
                            minDist = dist;
                            nearest = j;
                        }
                    }
                    if (nearest >= 0)
                        remap[k] = nearest;
                }
            }

            for (int64_t i = 0; i < total; ++i) {
                if (labels[i] >= 0)
                    labels[i] = remap[labels[i]];
            }
        }

        reportProgress(0.90, "Writing output...");

        Raster output(cols, rows, 1, DataType::Int16);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(-1);

        auto& dst = output.data(0);
        for (int64_t i = 0; i < total; ++i)
            dst[i] = static_cast<double>(labels[i]);

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ChainClusterModule)

} // namespace aplaceholder
