#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <queue>
#include <unordered_map>
#include <numeric>

namespace aplaceholder {

class SegmentModule : public Module {
public:
    QString name() const override { return "SEGMENT"; }
    QString description() const override {
        return "Image segmentation. Groups spatially contiguous pixels with similar "
               "spectral values into homogeneous regions.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output segmented image"),
            ParameterDef::real("similarity_threshold", "Similarity threshold", 10.0, 0.0, 1000.0,
                "Maximum spectral distance for merging adjacent pixels into segments"),
            ParameterDef::integer("min_segment_size", "Minimum segment size", 25, 1, 100000,
                "Minimum number of pixels per segment"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int numBands = raster->bands();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double threshold = parameter("similarity_threshold").toDouble();
        int minSize = parameter("min_segment_size").toInt();

        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        // Label array: -1 = unassigned
        std::vector<int> labels(total, -1);
        int nextLabel = 0;

        // 4-connected neighbor offsets
        const int dx[] = {-1, 1, 0, 0};
        const int dy[] = {0, 0, -1, 1};

        reportProgress(0.0, "Region growing...");

        // Phase 1: Region growing
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                if (labels[idx] >= 0) continue;

                // Check for nodata in any band
                bool isNoData = false;
                if (hasND) {
                    for (int b = 0; b < numBands; ++b) {
                        if (raster->data(b)[idx] == noData) {
                            isNoData = true;
                            break;
                        }
                    }
                }
                if (isNoData) {
                    labels[idx] = -2; // mark as nodata
                    continue;
                }

                // Start new region from this seed pixel
                int currentLabel = nextLabel++;
                labels[idx] = currentLabel;

                std::queue<int64_t> frontier;
                frontier.push(idx);

                while (!frontier.empty()) {
                    int64_t cur = frontier.front();
                    frontier.pop();
                    int cr = static_cast<int>(cur / cols);
                    int cc = static_cast<int>(cur % cols);

                    for (int d = 0; d < 4; ++d) {
                        int nr = cr + dy[d];
                        int nc = cc + dx[d];
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;

                        int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                        if (labels[nIdx] >= 0 || labels[nIdx] == -2) continue;

                        // Check nodata for neighbor
                        bool nIsNoData = false;
                        if (hasND) {
                            for (int b = 0; b < numBands; ++b) {
                                if (raster->data(b)[nIdx] == noData) {
                                    nIsNoData = true;
                                    break;
                                }
                            }
                        }
                        if (nIsNoData) {
                            labels[nIdx] = -2;
                            continue;
                        }

                        // Compute spectral distance between seed pixel and neighbor
                        double distSq = 0.0;
                        for (int b = 0; b < numBands; ++b) {
                            double diff = raster->data(b)[idx] - raster->data(b)[nIdx];
                            distSq += diff * diff;
                        }
                        double dist = std::sqrt(distSq);

                        if (dist <= threshold) {
                            labels[nIdx] = currentLabel;
                            frontier.push(nIdx);
                        }
                    }
                }
            }
            if (r % 100 == 0)
                reportProgress(0.5 * r / rows);
        }

        reportProgress(0.5, "Merging small segments...");

        // Phase 2: Merge regions smaller than min_size into nearest neighbor
        // Compute segment sizes
        std::unordered_map<int, int64_t> segSizes;
        for (int64_t i = 0; i < total; ++i) {
            if (labels[i] >= 0)
                segSizes[labels[i]]++;
        }

        // Compute segment mean spectral values
        std::unordered_map<int, std::vector<double>> segSums;
        for (auto& [label, size] : segSizes) {
            segSums[label].resize(numBands, 0.0);
        }
        for (int64_t i = 0; i < total; ++i) {
            if (labels[i] >= 0) {
                for (int b = 0; b < numBands; ++b) {
                    segSums[labels[i]][b] += raster->data(b)[i];
                }
            }
        }

        // Iteratively merge small segments into nearest neighboring segment
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto& [label, size] : segSizes) {
                if (size >= minSize || size == 0) continue;

                // Find neighboring segment labels
                std::unordered_map<int, bool> neighborLabels;
                for (int64_t i = 0; i < total; ++i) {
                    if (labels[i] != label) continue;
                    int r = static_cast<int>(i / cols);
                    int c = static_cast<int>(i % cols);
                    for (int d = 0; d < 4; ++d) {
                        int nr = r + dy[d];
                        int nc = c + dx[d];
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                        int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                        if (labels[nIdx] >= 0 && labels[nIdx] != label) {
                            neighborLabels[labels[nIdx]] = true;
                        }
                    }
                }

                if (neighborLabels.empty()) continue;

                // Find nearest neighbor by mean spectral distance
                int bestNeighbor = -1;
                double bestDist = 1e308;
                for (auto& [nLabel, _] : neighborLabels) {
                    double distSq = 0.0;
                    for (int b = 0; b < numBands; ++b) {
                        double meanA = segSums[label][b] / size;
                        double meanB = segSums[nLabel][b] / segSizes[nLabel];
                        double diff = meanA - meanB;
                        distSq += diff * diff;
                    }
                    if (distSq < bestDist) {
                        bestDist = distSq;
                        bestNeighbor = nLabel;
                    }
                }

                if (bestNeighbor < 0) continue;

                // Merge: relabel all pixels of this segment to bestNeighbor
                for (int64_t i = 0; i < total; ++i) {
                    if (labels[i] == label)
                        labels[i] = bestNeighbor;
                }
                for (int b = 0; b < numBands; ++b) {
                    segSums[bestNeighbor][b] += segSums[label][b];
                }
                segSizes[bestNeighbor] += size;
                segSums[label].assign(numBands, 0.0);
                size = 0;
                changed = true;
            }
        }

        reportProgress(0.9, "Writing output...");

        // Create output raster with segment labels
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(noData);

        auto& dst = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            dst[i] = (labels[i] >= 0) ? static_cast<double>(labels[i]) : noData;
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SegmentModule)

} // namespace aplaceholder
