#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <QStringList>

namespace aplaceholder {

class SegmentationModule : public Module {
public:
    QString name() const override { return "SEGMENTATION"; }
    QString description() const override {
        return "Advanced mean-shift segmentation. Performs iterative density-based "
               "clustering in joint spatial-spectral feature space, then labels "
               "connected components as segments with optional small-region merging.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Band images (comma-separated)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::output("output", "Output segment image"),
            ParameterDef::real("spatial_bandwidth", "Spatial bandwidth", 7.0, 1.0, 200.0,
                "Spatial search radius in pixels for mean-shift kernel"),
            ParameterDef::real("spectral_bandwidth", "Spectral bandwidth", 15.0, 1.0, 500.0,
                "Spectral distance threshold for mean-shift kernel"),
            ParameterDef::integer("min_region_size", "Minimum region size", 20, 1, 100000,
                "Minimum number of pixels per segment; smaller regions are merged"),
        };
    }

    bool execute() override {
        // Parse band paths
        QStringList bandPaths = parameter("bands").toString().split(",", Qt::SkipEmptyParts);
        int numBands = bandPaths.size();
        if (numBands == 0) {
            setError("No band images specified");
            return false;
        }

        // Read first band for dimensions
        auto firstBand = GdalIO::read(bandPaths[0].trimmed());
        if (!firstBand) {
            setError("Failed to read first band raster");
            return false;
        }

        int cols = firstBand->cols();
        int rows = firstBand->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Read all band data into a flat structure
        std::vector<std::vector<double>> bandData(numBands);
        bandData[0] = firstBand->data(0);
        auto geoTransform = firstBand->geoTransform();
        auto projection = firstBand->projection();
        bool hasND = firstBand->hasNoData();
        double noData = firstBand->noDataValue();
        firstBand.reset();

        for (int b = 1; b < numBands; ++b) {
            auto br = GdalIO::read(bandPaths[b].trimmed());
            if (!br) {
                setError("Failed to read band raster: " + bandPaths[b].trimmed());
                return false;
            }
            if (br->cols() != cols || br->rows() != rows) {
                setError("Band dimensions do not match: " + bandPaths[b].trimmed());
                return false;
            }
            bandData[b] = br->data(0);
        }

        double hs = parameter("spatial_bandwidth").toDouble();
        double hr = parameter("spectral_bandwidth").toDouble();
        int minRegionSize = parameter("min_region_size").toInt();

        double hs2 = hs * hs;
        double hr2 = hr * hr;
        int spatialRadius = static_cast<int>(std::ceil(hs));

        reportProgress(0.05, "Running mean-shift...");

        // Mean-shift filtering: for each pixel, iteratively shift to local mode
        // Store converged spectral values per pixel
        std::vector<std::vector<double>> converged(total, std::vector<double>(numBands));

        for (int64_t i = 0; i < total; ++i) {
            int py = static_cast<int>(i / cols);
            int px = static_cast<int>(i % cols);

            // Check nodata
            bool isND = false;
            if (hasND) {
                for (int b = 0; b < numBands; ++b) {
                    if (bandData[b][i] == noData) {
                        isND = true;
                        break;
                    }
                }
            }
            if (isND) {
                for (int b = 0; b < numBands; ++b)
                    converged[i][b] = noData;
                continue;
            }

            // Initialize mode at current pixel location + spectral values
            double sx = static_cast<double>(px);
            double sy = static_cast<double>(py);
            std::vector<double> sv(numBands);
            for (int b = 0; b < numBands; ++b)
                sv[b] = bandData[b][i];

            // Iterate mean-shift
            const int maxIter = 20;
            const double convergenceThreshold = 0.5;

            for (int iter = 0; iter < maxIter; ++iter) {
                double sumWeight = 0.0;
                double newSx = 0.0, newSy = 0.0;
                std::vector<double> newSv(numBands, 0.0);

                int rMin = std::max(0, static_cast<int>(sy) - spatialRadius);
                int rMax = std::min(rows - 1, static_cast<int>(sy) + spatialRadius);
                int cMin = std::max(0, static_cast<int>(sx) - spatialRadius);
                int cMax = std::min(cols - 1, static_cast<int>(sx) + spatialRadius);

                for (int ry = rMin; ry <= rMax; ++ry) {
                    for (int rx = cMin; rx <= cMax; ++rx) {
                        int64_t nIdx = static_cast<int64_t>(ry) * cols + rx;

                        // Check nodata
                        bool nIsND = false;
                        if (hasND) {
                            for (int b = 0; b < numBands; ++b) {
                                if (bandData[b][nIdx] == noData) {
                                    nIsND = true;
                                    break;
                                }
                            }
                        }
                        if (nIsND) continue;

                        // Spatial distance
                        double dsx = rx - sx;
                        double dsy = ry - sy;
                        double spatDist2 = dsx * dsx + dsy * dsy;
                        if (spatDist2 > hs2) continue;

                        // Spectral distance
                        double specDist2 = 0.0;
                        for (int b = 0; b < numBands; ++b) {
                            double diff = bandData[b][nIdx] - sv[b];
                            specDist2 += diff * diff;
                        }
                        if (specDist2 > hr2) continue;

                        // Epanechnikov-like kernel (uniform within bandwidth)
                        double weight = 1.0;
                        sumWeight += weight;
                        newSx += weight * rx;
                        newSy += weight * ry;
                        for (int b = 0; b < numBands; ++b)
                            newSv[b] += weight * bandData[b][nIdx];
                    }
                }

                if (sumWeight == 0.0) break;

                newSx /= sumWeight;
                newSy /= sumWeight;
                for (int b = 0; b < numBands; ++b)
                    newSv[b] /= sumWeight;

                // Check convergence
                double shift = (newSx - sx) * (newSx - sx) + (newSy - sy) * (newSy - sy);
                for (int b = 0; b < numBands; ++b) {
                    double d = newSv[b] - sv[b];
                    shift += d * d;
                }

                sx = newSx;
                sy = newSy;
                sv = newSv;

                if (shift < convergenceThreshold) break;
            }

            converged[i] = sv;

            if (i % (total / 40 + 1) == 0)
                reportProgress(0.05 + 0.6 * static_cast<double>(i) / total,
                               "Mean-shift filtering...");
        }

        reportProgress(0.65, "Labeling connected components...");

        // Phase 2: Label connected components based on converged spectral similarity
        std::vector<int> labels(total, -1);
        int nextLabel = 0;
        const int dx[] = {-1, 1, 0, 0};
        const int dy[] = {0, 0, -1, 1};

        for (int64_t i = 0; i < total; ++i) {
            if (labels[i] >= 0) continue;

            // Check nodata
            bool isND = false;
            if (hasND) {
                for (int b = 0; b < numBands; ++b) {
                    if (converged[i][b] == noData) {
                        isND = true;
                        break;
                    }
                }
            }
            if (isND) {
                labels[i] = -2;
                continue;
            }

            int currentLabel = nextLabel++;
            labels[i] = currentLabel;

            std::vector<int64_t> stack;
            stack.push_back(i);

            while (!stack.empty()) {
                int64_t cur = stack.back();
                stack.pop_back();
                int cr = static_cast<int>(cur / cols);
                int cc = static_cast<int>(cur % cols);

                for (int d = 0; d < 4; ++d) {
                    int nr = cr + dy[d];
                    int nc = cc + dx[d];
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;

                    int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                    if (labels[nIdx] >= 0 || labels[nIdx] == -2) continue;

                    // Check nodata
                    bool nIsND = false;
                    if (hasND) {
                        for (int b = 0; b < numBands; ++b) {
                            if (converged[nIdx][b] == noData) {
                                nIsND = true;
                                break;
                            }
                        }
                    }
                    if (nIsND) {
                        labels[nIdx] = -2;
                        continue;
                    }

                    // Check if converged spectral values are close enough
                    double specDist2 = 0.0;
                    for (int b = 0; b < numBands; ++b) {
                        double diff = converged[cur][b] - converged[nIdx][b];
                        specDist2 += diff * diff;
                    }

                    if (specDist2 <= hr2 * 0.25) { // half-bandwidth for component labeling
                        labels[nIdx] = currentLabel;
                        stack.push_back(nIdx);
                    }
                }
            }
        }

        reportProgress(0.8, "Merging small regions...");

        // Phase 3: Merge small regions
        std::unordered_map<int, int64_t> segSizes;
        for (int64_t i = 0; i < total; ++i) {
            if (labels[i] >= 0)
                segSizes[labels[i]]++;
        }

        // Compute mean converged values per segment
        std::unordered_map<int, std::vector<double>> segMeans;
        for (auto& [label, size] : segSizes) {
            segMeans[label].resize(numBands, 0.0);
        }
        for (int64_t i = 0; i < total; ++i) {
            if (labels[i] >= 0) {
                for (int b = 0; b < numBands; ++b)
                    segMeans[labels[i]][b] += converged[i][b];
            }
        }
        for (auto& [label, means] : segMeans) {
            int64_t sz = segSizes[label];
            if (sz > 0) {
                for (int b = 0; b < numBands; ++b)
                    means[b] /= sz;
            }
        }

        // Iteratively merge small segments
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto& [label, size] : segSizes) {
                if (size >= minRegionSize || size == 0) continue;

                // Find neighboring segments
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
                        if (labels[nIdx] >= 0 && labels[nIdx] != label)
                            neighborLabels[labels[nIdx]] = true;
                    }
                }

                if (neighborLabels.empty()) continue;

                int bestNeighbor = -1;
                double bestDist = 1e308;
                for (auto& [nLabel, _] : neighborLabels) {
                    double dist2 = 0.0;
                    for (int b = 0; b < numBands; ++b) {
                        double diff = segMeans[label][b] - segMeans[nLabel][b];
                        dist2 += diff * diff;
                    }
                    if (dist2 < bestDist) {
                        bestDist = dist2;
                        bestNeighbor = nLabel;
                    }
                }

                if (bestNeighbor < 0) continue;

                // Merge into best neighbor
                int64_t totalSize = size + segSizes[bestNeighbor];
                for (int b = 0; b < numBands; ++b) {
                    segMeans[bestNeighbor][b] =
                        (segMeans[bestNeighbor][b] * segSizes[bestNeighbor] +
                         segMeans[label][b] * size) / totalSize;
                }
                for (int64_t i = 0; i < total; ++i) {
                    if (labels[i] == label)
                        labels[i] = bestNeighbor;
                }
                segSizes[bestNeighbor] = totalSize;
                size = 0;
                changed = true;
            }
        }

        reportProgress(0.95, "Writing output...");

        // Create output raster
        double outND = -9999.0;
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(geoTransform);
        output.setProjection(projection);
        output.setNoDataValue(outND);

        auto& outData = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            outData[i] = (labels[i] >= 0) ? static_cast<double>(labels[i]) : outND;
        }

        reportProgress(1.0, "Segmentation complete.");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SegmentationModule)

} // namespace aplaceholder
