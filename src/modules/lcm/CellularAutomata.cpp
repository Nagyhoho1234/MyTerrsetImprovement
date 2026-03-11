#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class CellularAutomataModule : public Module {
public:
    QString name() const override { return "CELLULAR_AUTOMATA"; }
    QString description() const override {
        return "CA-Markov land change simulation. Combines cellular automata "
               "spatial contiguity filters with Markov chain transition probabilities "
               "and suitability maps to simulate land cover change.";
    }
    QString category() const override { return "Land Change Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("base_image", "Base land cover image"),
            ParameterDef::file("transition_matrix", "Transition probability matrix (CSV)"),
            ParameterDef::file("suitability_maps", "Suitability maps (comma-separated)"),
            ParameterDef::integer("iterations", "Number of iterations", 10, 1, 1000,
                "Number of cellular automata iterations to run"),
            ParameterDef::output("output", "Output simulated land cover"),
        };
    }

    bool execute() override {
        // Read base land cover
        auto baseRaster = GdalIO::read(parameter("base_image").toString());
        if (!baseRaster) {
            setError("Failed to read base land cover image");
            return false;
        }

        int cols = baseRaster->cols(), rows = baseRaster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = baseRaster->noDataValue();
        bool hasND = baseRaster->hasNoData();

        // Parse transition matrix CSV
        std::vector<int> classes;
        std::vector<std::vector<double>> transMatrix;
        if (!readTransitionMatrix(parameter("transition_matrix").toString(),
                                  classes, transMatrix)) {
            return false;
        }
        int n = static_cast<int>(classes.size());

        std::map<int, int> classIndex;
        for (int i = 0; i < n; ++i)
            classIndex[classes[i]] = i;

        // Parse and read suitability maps (comma-separated file paths)
        // One suitability map per class, ordered to match classes in the matrix
        QString suitPaths = parameter("suitability_maps").toString();
        QStringList suitFiles = suitPaths.split(",", Qt::SkipEmptyParts);

        std::vector<std::shared_ptr<Raster>> suitRasters;
        for (const auto& path : suitFiles) {
            auto sr = GdalIO::read(path.trimmed());
            if (!sr) {
                setError("Failed to read suitability map: " + path.trimmed());
                return false;
            }
            if (sr->cols() != cols || sr->rows() != rows) {
                setError("Suitability map dimensions do not match base image: " + path.trimmed());
                return false;
            }
            suitRasters.push_back(std::move(sr));
        }

        // If fewer suitability maps than classes, fill remaining with uniform suitability
        while (static_cast<int>(suitRasters.size()) < n) {
            auto uniform = std::make_shared<Raster>(cols, rows, 1, DataType::Float64);
            auto& ud = uniform->data(0);
            for (int64_t i = 0; i < total; ++i)
                ud[i] = 1.0;
            suitRasters.push_back(uniform);
        }

        int iterations = parameter("iterations").toInt();

        // Working copy of the land cover map
        std::vector<double> current(total);
        const auto& baseData = baseRaster->data(0);
        for (int64_t i = 0; i < total; ++i)
            current[i] = baseData[i];

        // Compute expected pixel counts from Markov matrix
        // Count current class areas
        std::vector<int64_t> initialCounts(n, 0);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && current[i] == noData) continue;
            auto it = classIndex.find(static_cast<int>(current[i]));
            if (it != classIndex.end())
                initialCounts[it->second]++;
        }

        // Target counts for each class after full transition
        std::vector<int64_t> targetCounts(n, 0);
        for (int ci = 0; ci < n; ++ci) {
            double expected = 0.0;
            for (int cj = 0; cj < n; ++cj)
                expected += initialCounts[cj] * transMatrix[cj][ci];
            targetCounts[ci] = static_cast<int64_t>(std::round(expected));
        }

        // Iterative CA-Markov simulation
        std::vector<double> next(total);
        for (int iter = 0; iter < iterations; ++iter) {
            double iterFrac = static_cast<double>(iter + 1) / iterations;
            reportProgress(iterFrac * 0.9, "Iteration " + QString::number(iter + 1));

            // Interpolated target counts for this iteration
            std::vector<int64_t> currentCounts(n, 0);
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && current[i] == noData) continue;
                auto it = classIndex.find(static_cast<int>(current[i]));
                if (it != classIndex.end())
                    currentCounts[it->second]++;
            }

            // Per-iteration incremental target: move fractionally toward final target
            std::vector<int64_t> iterTarget(n);
            for (int ci = 0; ci < n; ++ci) {
                double stepTarget = initialCounts[ci] +
                    iterFrac * (targetCounts[ci] - initialCounts[ci]);
                iterTarget[ci] = static_cast<int64_t>(std::round(stepTarget));
            }

            // Determine how many pixels each class needs to gain/lose this iteration
            std::vector<int64_t> deficit(n);
            for (int ci = 0; ci < n; ++ci)
                deficit[ci] = iterTarget[ci] - currentCounts[ci];

            // Compute allocation scores for each pixel
            // For pixels eligible to change: score = suitability * transition_probability
            // Competitive multi-objective allocation: classes needing more pixels claim first

            // Copy current state
            for (int64_t i = 0; i < total; ++i)
                next[i] = current[i];

            // Build allocation candidates per target class that needs to gain pixels
            for (int toIdx = 0; toIdx < n; ++toIdx) {
                if (deficit[toIdx] <= 0) continue;  // Class doesn't need to gain

                // Collect candidate pixels (currently not this class, with positive
                // transition probability from their current class to toIdx)
                struct Candidate {
                    int64_t pixelIdx;
                    double score;
                };
                std::vector<Candidate> candidates;

                for (int64_t i = 0; i < total; ++i) {
                    if (hasND && current[i] == noData) continue;
                    int curClass = static_cast<int>(current[i]);
                    auto it = classIndex.find(curClass);
                    if (it == classIndex.end()) continue;
                    int fromIdx = it->second;
                    if (fromIdx == toIdx) continue;  // Already this class

                    double prob = transMatrix[fromIdx][toIdx];
                    if (prob <= 0.0) continue;

                    double suit = suitRasters[toIdx]->data(0)[i];
                    double score = suit * prob;
                    if (score > 0.0)
                        candidates.push_back({i, score});
                }

                // Sort by score descending
                std::sort(candidates.begin(), candidates.end(),
                    [](const Candidate& a, const Candidate& b) {
                        return a.score > b.score;
                    });

                // Allocate top deficit[toIdx] pixels
                int64_t toAllocate = std::min(deficit[toIdx],
                    static_cast<int64_t>(candidates.size()));
                for (int64_t k = 0; k < toAllocate; ++k)
                    next[candidates[k].pixelIdx] = static_cast<double>(classes[toIdx]);
            }

            // Apply 5x5 contiguity filter to prevent scattered isolated changes
            applyContiguityFilter(next, current, cols, rows, noData, hasND, classIndex);

            // Update current state
            for (int64_t i = 0; i < total; ++i)
                current[i] = next[i];
        }

        // Write output raster
        reportProgress(0.95, "Writing output...");
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(baseRaster->geoTransform());
        output.setProjection(baseRaster->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);
        for (int64_t i = 0; i < total; ++i)
            outData[i] = current[i];

        if (!GdalIO::write(output, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0, "CA-Markov simulation complete.");
        return true;
    }

private:
    bool readTransitionMatrix(const QString& path,
                              std::vector<int>& classes,
                              std::vector<std::vector<double>>& matrix) {
        std::ifstream file(path.toStdString());
        if (!file.is_open()) {
            setError("Failed to open transition matrix file: " + path);
            return false;
        }

        std::string line;
        // Read header to get class labels
        if (!std::getline(file, line)) {
            setError("Empty transition matrix file");
            return false;
        }

        std::stringstream headerSS(line);
        std::string token;
        // Skip first column header (From\To)
        std::getline(headerSS, token, ',');
        while (std::getline(headerSS, token, ',')) {
            token.erase(0, token.find_first_not_of(" \t\r\n"));
            token.erase(token.find_last_not_of(" \t\r\n") + 1);
            if (!token.empty())
                classes.push_back(std::stoi(token));
        }

        int n = static_cast<int>(classes.size());
        matrix.resize(n, std::vector<double>(n, 0.0));

        int row = 0;
        while (std::getline(file, line) && row < n) {
            std::stringstream ss(line);
            // Skip row label
            std::getline(ss, token, ',');
            for (int col = 0; col < n && std::getline(ss, token, ','); ++col) {
                token.erase(0, token.find_first_not_of(" \t\r\n"));
                token.erase(token.find_last_not_of(" \t\r\n") + 1);
                if (!token.empty())
                    matrix[row][col] = std::stod(token);
            }
            ++row;
        }

        file.close();
        return true;
    }

    // 5x5 contiguity filter: revert changes where the new class does not have
    // enough neighbors of the same class in the 5x5 neighborhood
    void applyContiguityFilter(std::vector<double>& next,
                               const std::vector<double>& current,
                               int cols, int rows,
                               double noData, bool hasND,
                               const std::map<int, int>& classIndex) {
        int64_t total = static_cast<int64_t>(cols) * rows;
        std::vector<double> filtered(next.begin(), next.end());

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                if (hasND && next[idx] == noData) continue;

                // Only check pixels that changed this iteration
                if (next[idx] == current[idx]) continue;

                int newClass = static_cast<int>(next[idx]);

                // Count neighbors of same new class in 5x5 window
                int sameCount = 0;
                int validCount = 0;
                for (int dr = -2; dr <= 2; ++dr) {
                    for (int dc = -2; dc <= 2; ++dc) {
                        if (dr == 0 && dc == 0) continue;
                        int nr = r + dr, nc = c + dc;
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                        int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                        if (hasND && next[nIdx] == noData) continue;
                        validCount++;
                        if (static_cast<int>(next[nIdx]) == newClass)
                            sameCount++;
                    }
                }

                // Require at least some contiguity: revert if isolated
                // Threshold: at least 2 neighbors of same class in 5x5 window
                if (sameCount < 2)
                    filtered[idx] = current[idx];
            }
        }

        next = filtered;
    }
};

REGISTER_MODULE(CellularAutomataModule)

} // namespace aplaceholder
