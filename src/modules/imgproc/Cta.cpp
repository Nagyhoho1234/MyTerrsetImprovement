#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>
#include <memory>

namespace aplaceholder {

class CtaModule : public Module {
public:
    QString name() const override { return "CTA"; }
    QString description() const override {
        return "Classification Tree Analysis. Builds a binary decision tree by recursive "
               "splitting on band thresholds using Gini impurity criterion.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster where pixel values indicate class IDs (0 = unclassified)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("max_depth", "Maximum tree depth", 10, 1, 50,
                "Maximum depth of the decision tree"),
            ParameterDef::integer("min_samples", "Minimum samples per leaf", 5, 1, 10000,
                "Minimum number of samples required to create a leaf node"),
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

        int maxDepth = parameter("max_depth").toInt();
        int minSamples = parameter("min_samples").toInt();

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

        if (numClasses < 2) {
            setError("Need at least 2 classes in training data");
            return false;
        }

        reportProgress(0.1, "Building decision tree...");

        // --------------------------------------------------------------------
        // 4. Build decision tree
        // --------------------------------------------------------------------
        struct TreeNode {
            bool isLeaf = false;
            int classId = 0;          // for leaf nodes
            int splitBand = -1;        // for internal nodes
            double splitThreshold = 0; // for internal nodes
            int left = -1;             // child index: <= threshold
            int right = -1;            // child index: > threshold
        };

        std::vector<TreeNode> tree;
        tree.reserve(1024);

        // Gini impurity for a set of samples
        auto gini = [&classIds, numClasses](const std::vector<int>& indices,
                                             const std::vector<Sample>& allSamples) -> double {
            if (indices.empty()) return 0.0;
            std::vector<int> counts(numClasses, 0);
            for (int idx : indices) {
                int cls = allSamples[idx].classId;
                for (int ci = 0; ci < numClasses; ++ci) {
                    if (classIds[ci] == cls) {
                        counts[ci]++;
                        break;
                    }
                }
            }
            double n = static_cast<double>(indices.size());
            double impurity = 1.0;
            for (int ci = 0; ci < numClasses; ++ci) {
                double p = counts[ci] / n;
                impurity -= p * p;
            }
            return impurity;
        };

        // Majority class for a set of samples
        auto majorityClass = [&classIds, numClasses](const std::vector<int>& indices,
                                                      const std::vector<Sample>& allSamples) -> int {
            std::vector<int> counts(numClasses, 0);
            for (int idx : indices) {
                int cls = allSamples[idx].classId;
                for (int ci = 0; ci < numClasses; ++ci) {
                    if (classIds[ci] == cls) {
                        counts[ci]++;
                        break;
                    }
                }
            }
            int bestIdx = 0;
            for (int ci = 1; ci < numClasses; ++ci) {
                if (counts[ci] > counts[bestIdx])
                    bestIdx = ci;
            }
            return classIds[bestIdx];
        };

        // Build tree recursively using a stack
        struct BuildTask {
            std::vector<int> indices;
            int nodeIdx;
            int depth;
        };

        // Create root node
        tree.push_back(TreeNode{});
        std::vector<int> rootIndices(samples.size());
        std::iota(rootIndices.begin(), rootIndices.end(), 0);

        std::vector<BuildTask> stack;
        stack.push_back({std::move(rootIndices), 0, 0});

        while (!stack.empty()) {
            BuildTask task = std::move(stack.back());
            stack.pop_back();

            int nIdx = static_cast<int>(task.indices.size());

            // Check termination conditions
            if (task.depth >= maxDepth || nIdx <= minSamples || gini(task.indices, samples) < 1e-10) {
                tree[task.nodeIdx].isLeaf = true;
                tree[task.nodeIdx].classId = majorityClass(task.indices, samples);
                continue;
            }

            // Find best split across all bands and thresholds
            double bestGain = -1.0;
            int bestBand = 0;
            double bestThreshold = 0.0;
            double parentGini = gini(task.indices, samples);

            for (int b = 0; b < numBands; ++b) {
                // Sort indices by feature value in this band
                std::vector<int> sorted = task.indices;
                std::sort(sorted.begin(), sorted.end(),
                    [&samples, b](int a, int bb) {
                        return samples[a].features[b] < samples[bb].features[b];
                    });

                // Try split points between consecutive distinct values
                for (int s = 1; s < nIdx; ++s) {
                    if (samples[sorted[s]].features[b] == samples[sorted[s - 1]].features[b])
                        continue;

                    double threshold = (samples[sorted[s - 1]].features[b] +
                                       samples[sorted[s]].features[b]) / 2.0;

                    std::vector<int> leftIdx(sorted.begin(), sorted.begin() + s);
                    std::vector<int> rightIdx(sorted.begin() + s, sorted.end());

                    double leftWeight = static_cast<double>(s) / nIdx;
                    double rightWeight = static_cast<double>(nIdx - s) / nIdx;
                    double gain = parentGini
                                  - leftWeight * gini(leftIdx, samples)
                                  - rightWeight * gini(rightIdx, samples);

                    if (gain > bestGain) {
                        bestGain = gain;
                        bestBand = b;
                        bestThreshold = threshold;
                    }
                }
            }

            if (bestGain <= 0) {
                tree[task.nodeIdx].isLeaf = true;
                tree[task.nodeIdx].classId = majorityClass(task.indices, samples);
                continue;
            }

            // Perform the split
            std::vector<int> leftIdx, rightIdx;
            for (int idx : task.indices) {
                if (samples[idx].features[bestBand] <= bestThreshold)
                    leftIdx.push_back(idx);
                else
                    rightIdx.push_back(idx);
            }

            tree[task.nodeIdx].splitBand = bestBand;
            tree[task.nodeIdx].splitThreshold = bestThreshold;

            int leftNode = static_cast<int>(tree.size());
            tree.push_back(TreeNode{});
            int rightNode = static_cast<int>(tree.size());
            tree.push_back(TreeNode{});

            tree[task.nodeIdx].left = leftNode;
            tree[task.nodeIdx].right = rightNode;

            stack.push_back({std::move(leftIdx), leftNode, task.depth + 1});
            stack.push_back({std::move(rightIdx), rightNode, task.depth + 1});
        }

        // --------------------------------------------------------------------
        // 5. Classify all pixels
        // --------------------------------------------------------------------
        reportProgress(0.6, "Classifying pixels...");

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

            // Traverse the tree
            int nodeIdx = 0;
            while (!tree[nodeIdx].isLeaf) {
                double val = (*bands[tree[nodeIdx].splitBand])[i];
                if (val <= tree[nodeIdx].splitThreshold)
                    nodeIdx = tree[nodeIdx].left;
                else
                    nodeIdx = tree[nodeIdx].right;
            }
            out[i] = tree[nodeIdx].classId;

            if (i % 1000000 == 0)
                reportProgress(0.6 + 0.35 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(CtaModule)

} // namespace aplaceholder
