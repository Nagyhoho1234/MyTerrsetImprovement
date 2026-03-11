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

class DecisionForestModule : public Module {
public:
    QString name() const override { return "DECISIONFOREST"; }
    QString description() const override {
        return "Random Forest classifier. Ensemble of decision trees with random feature subsets "
               "and bootstrap sampling. Classification by majority vote across all trees.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster where pixel values indicate class IDs (0 = unclassified)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::integer("num_trees", "Number of trees", 50, 1, 1000,
                "Number of decision trees in the forest"),
            ParameterDef::integer("max_depth", "Maximum tree depth", 10, 1, 50,
                "Maximum depth of each decision tree"),
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

        int numTrees = parameter("num_trees").toInt();
        int maxDepth = parameter("max_depth").toInt();

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

        // Number of features to consider at each split: sqrt(numBands)
        int featuresPerSplit = std::max(1, static_cast<int>(std::sqrt(numBands)));

        reportProgress(0.05, "Building random forest...");

        // --------------------------------------------------------------------
        // 4. Build tree data structures
        // --------------------------------------------------------------------
        struct TreeNode {
            bool isLeaf = false;
            int classId = 0;
            int splitBand = -1;
            double splitThreshold = 0;
            int left = -1;
            int right = -1;
        };

        // All trees stored in a single vector of vectors
        std::vector<std::vector<TreeNode>> forest(numTrees);

        std::mt19937 rng(42);

        // Helper lambdas
        auto giniForIndices = [&classIds, numClasses, &samples](const std::vector<int>& indices) -> double {
            if (indices.empty()) return 0.0;
            std::vector<int> counts(numClasses, 0);
            for (int idx : indices) {
                for (int ci = 0; ci < numClasses; ++ci) {
                    if (classIds[ci] == samples[idx].classId) {
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

        auto majorityClassForIndices = [&classIds, numClasses, &samples](const std::vector<int>& indices) -> int {
            std::vector<int> counts(numClasses, 0);
            for (int idx : indices) {
                for (int ci = 0; ci < numClasses; ++ci) {
                    if (classIds[ci] == samples[idx].classId) {
                        counts[ci]++;
                        break;
                    }
                }
            }
            int bestIdx = 0;
            for (int ci = 1; ci < numClasses; ++ci) {
                if (counts[ci] > counts[bestIdx]) bestIdx = ci;
            }
            return classIds[bestIdx];
        };

        // --------------------------------------------------------------------
        // 5. Build each tree
        // --------------------------------------------------------------------
        for (int t = 0; t < numTrees; ++t) {
            // Bootstrap sample
            std::vector<int> bootstrap(nSamples);
            std::uniform_int_distribution<int> sampleDist(0, nSamples - 1);
            for (int i = 0; i < nSamples; ++i)
                bootstrap[i] = sampleDist(rng);

            auto& tree = forest[t];
            tree.reserve(256);
            tree.push_back(TreeNode{});

            struct BuildTask {
                std::vector<int> indices;
                int nodeIdx;
                int depth;
            };

            std::vector<BuildTask> buildStack;
            buildStack.push_back({std::move(bootstrap), 0, 0});

            while (!buildStack.empty()) {
                BuildTask task = std::move(buildStack.back());
                buildStack.pop_back();

                int nIdx = static_cast<int>(task.indices.size());

                if (task.depth >= maxDepth || nIdx <= 2 ||
                    giniForIndices(task.indices) < 1e-10) {
                    tree[task.nodeIdx].isLeaf = true;
                    tree[task.nodeIdx].classId = majorityClassForIndices(task.indices);
                    continue;
                }

                // Select random subset of bands
                std::vector<int> bandSubset(numBands);
                std::iota(bandSubset.begin(), bandSubset.end(), 0);
                std::shuffle(bandSubset.begin(), bandSubset.end(), rng);
                bandSubset.resize(featuresPerSplit);

                double bestGain = -1.0;
                int bestBand = bandSubset[0];
                double bestThreshold = 0.0;
                double parentGini = giniForIndices(task.indices);

                for (int bi : bandSubset) {
                    std::vector<int> sorted = task.indices;
                    std::sort(sorted.begin(), sorted.end(),
                        [&samples, bi](int a, int b) {
                            return samples[a].features[bi] < samples[b].features[bi];
                        });

                    // Try a subset of split points (every step-th) for speed
                    int step = std::max(1, nIdx / 20);
                    for (int s = step; s < nIdx; s += step) {
                        if (samples[sorted[s]].features[bi] ==
                            samples[sorted[s - 1]].features[bi])
                            continue;

                        double threshold = (samples[sorted[s - 1]].features[bi] +
                                           samples[sorted[s]].features[bi]) / 2.0;

                        std::vector<int> leftIdx(sorted.begin(), sorted.begin() + s);
                        std::vector<int> rightIdx(sorted.begin() + s, sorted.end());

                        double leftWeight = static_cast<double>(s) / nIdx;
                        double rightWeight = static_cast<double>(nIdx - s) / nIdx;
                        double gain = parentGini
                                      - leftWeight * giniForIndices(leftIdx)
                                      - rightWeight * giniForIndices(rightIdx);

                        if (gain > bestGain) {
                            bestGain = gain;
                            bestBand = bi;
                            bestThreshold = threshold;
                        }
                    }
                }

                if (bestGain <= 0) {
                    tree[task.nodeIdx].isLeaf = true;
                    tree[task.nodeIdx].classId = majorityClassForIndices(task.indices);
                    continue;
                }

                std::vector<int> leftIdx, rightIdx;
                for (int idx : task.indices) {
                    if (samples[idx].features[bestBand] <= bestThreshold)
                        leftIdx.push_back(idx);
                    else
                        rightIdx.push_back(idx);
                }

                if (leftIdx.empty() || rightIdx.empty()) {
                    tree[task.nodeIdx].isLeaf = true;
                    tree[task.nodeIdx].classId = majorityClassForIndices(task.indices);
                    continue;
                }

                tree[task.nodeIdx].splitBand = bestBand;
                tree[task.nodeIdx].splitThreshold = bestThreshold;

                int leftNode = static_cast<int>(tree.size());
                tree.push_back(TreeNode{});
                int rightNode = static_cast<int>(tree.size());
                tree.push_back(TreeNode{});

                tree[task.nodeIdx].left = leftNode;
                tree[task.nodeIdx].right = rightNode;

                buildStack.push_back({std::move(leftIdx), leftNode, task.depth + 1});
                buildStack.push_back({std::move(rightIdx), rightNode, task.depth + 1});
            }

            if (t % 5 == 0)
                reportProgress(0.05 + 0.55 * static_cast<double>(t + 1) / numTrees,
                               QString("Built tree %1 of %2").arg(t + 1).arg(numTrees));
        }

        // --------------------------------------------------------------------
        // 6. Classify all pixels by majority vote
        // --------------------------------------------------------------------
        reportProgress(0.6, "Classifying pixels...");

        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);
        std::fill(out.begin(), out.end(), 0.0);

        // Build class id to vote index map
        std::vector<int> classToVoteIdx(classIds.back() + 1, -1);
        for (int ci = 0; ci < numClasses; ++ci)
            classToVoteIdx[classIds[ci]] = ci;

        std::vector<int> votes(numClasses);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                out[i] = 0;
                continue;
            }

            std::fill(votes.begin(), votes.end(), 0);

            for (int t = 0; t < numTrees; ++t) {
                const auto& tree = forest[t];
                int nodeIdx = 0;
                while (!tree[nodeIdx].isLeaf) {
                    double val = (*bands[tree[nodeIdx].splitBand])[i];
                    if (val <= tree[nodeIdx].splitThreshold)
                        nodeIdx = tree[nodeIdx].left;
                    else
                        nodeIdx = tree[nodeIdx].right;
                }
                int cls = tree[nodeIdx].classId;
                if (cls > 0 && cls < static_cast<int>(classToVoteIdx.size()) &&
                    classToVoteIdx[cls] >= 0) {
                    votes[classToVoteIdx[cls]]++;
                }
            }

            int bestIdx = 0;
            for (int ci = 1; ci < numClasses; ++ci) {
                if (votes[ci] > votes[bestIdx]) bestIdx = ci;
            }
            out[i] = classIds[bestIdx];

            if (i % 1000000 == 0)
                reportProgress(0.6 + 0.35 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(DecisionForestModule)

} // namespace aplaceholder
