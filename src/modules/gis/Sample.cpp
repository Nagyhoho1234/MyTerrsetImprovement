#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <set>

namespace aplaceholder {

class SampleModule : public Module {
public:
    QString name() const override { return "SAMPLE"; }
    QString description() const override {
        return "Generates sample point locations using random, systematic, or "
               "stratified random sampling. Creates a point raster with sampling "
               "locations for ground truth collection and error assessment.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("reference", "Reference raster (for extent/dimensions)"),
            ParameterDef::output("output", "Output sample point raster"),
            ParameterDef::combo("method", "Sampling method",
                {"Random", "Systematic", "Stratified Random"}, 0,
                "Sampling scheme to use"),
            ParameterDef::integer("num_points", "Number of sample points",
                100, 1, 999999,
                "Total number of sample points to generate"),
            ParameterDef::file("strata", "Strata raster (optional)",
                "Categorical raster defining strata for stratified sampling"),
        };
    }

    bool execute() override {
        auto ref = GdalIO::read(parameter("reference").toString());
        if (!ref) {
            setError("Failed to read reference raster");
            return false;
        }

        int cols = ref->cols(), rows = ref->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int numPoints = parameter("num_points").toInt();
        int method = parameter("method").toInt();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(ref->geoTransform());
        output.setProjection(ref->projection());
        output.setNoDataValue(0);
        auto& out = output.data(0);

        for (int64_t i = 0; i < total; ++i)
            out[i] = 0.0;

        reportProgress(0.1, "Generating sample points...");

        if (method == 0) {
            // Random sampling
            generateRandom(out, cols, rows, numPoints);
        } else if (method == 1) {
            // Systematic sampling
            generateSystematic(out, cols, rows, numPoints);
        } else {
            // Stratified random sampling
            QString strataPath = parameter("strata").toString();
            if (strataPath.isEmpty()) {
                // Without strata raster: spatially stratified by grid division
                generateStratifiedGrid(out, cols, rows, numPoints);
            } else {
                auto strata = GdalIO::read(strataPath);
                if (!strata) {
                    setError("Failed to read strata raster");
                    return false;
                }
                if (strata->cols() != cols || strata->rows() != rows) {
                    setError("Strata raster must have same dimensions as reference");
                    return false;
                }
                generateStratifiedByRaster(out, strata->data(0), cols, rows, numPoints,
                                           strata->noDataValue(), strata->hasNoData());
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }

private:
    void generateRandom(std::vector<double>& out, int cols, int rows, int numPoints) {
        int64_t total = static_cast<int64_t>(cols) * rows;
        std::mt19937_64 rng(std::random_device{}());
        std::uniform_int_distribution<int64_t> dist(0, total - 1);

        std::set<int64_t> selected;
        int pointId = 1;
        int attempts = 0;
        int maxAttempts = numPoints * 20;

        while (static_cast<int>(selected.size()) < numPoints && attempts < maxAttempts) {
            int64_t idx = dist(rng);
            if (selected.insert(idx).second) {
                out[idx] = static_cast<double>(pointId++);
            }
            attempts++;
        }
    }

    void generateSystematic(std::vector<double>& out, int cols, int rows, int numPoints) {
        // Calculate grid spacing to distribute points evenly
        double area = static_cast<double>(cols) * rows;
        double cellsPerPoint = area / numPoints;
        double spacing = std::sqrt(cellsPerPoint);

        int spacingC = std::max(1, static_cast<int>(std::round(spacing)));
        int spacingR = std::max(1, static_cast<int>(std::round(spacing)));

        // Offset to center the grid
        int offsetC = spacingC / 2;
        int offsetR = spacingR / 2;

        int pointId = 1;
        for (int r = offsetR; r < rows; r += spacingR) {
            for (int c = offsetC; c < cols; c += spacingC) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                out[idx] = static_cast<double>(pointId++);
                if (pointId > numPoints + 1) return;
            }
        }
    }

    void generateStratifiedGrid(std::vector<double>& out, int cols, int rows, int numPoints) {
        // Divide the area into approximately numPoints rectangular strata
        // and place one random point in each stratum
        double area = static_cast<double>(cols) * rows;
        double cellsPerStratum = area / numPoints;
        double stratumSide = std::sqrt(cellsPerStratum);

        int numStrataC = std::max(1, static_cast<int>(std::round(static_cast<double>(cols) / stratumSide)));
        int numStrataR = std::max(1, static_cast<int>(std::round(static_cast<double>(rows) / stratumSide)));

        double stratumW = static_cast<double>(cols) / numStrataC;
        double stratumH = static_cast<double>(rows) / numStrataR;

        std::mt19937_64 rng(std::random_device{}());
        int pointId = 1;

        for (int sr = 0; sr < numStrataR; ++sr) {
            for (int sc = 0; sc < numStrataC; ++sc) {
                int minC = static_cast<int>(sc * stratumW);
                int maxC = std::min(cols - 1, static_cast<int>((sc + 1) * stratumW - 1));
                int minR = static_cast<int>(sr * stratumH);
                int maxR = std::min(rows - 1, static_cast<int>((sr + 1) * stratumH - 1));

                if (maxC < minC) maxC = minC;
                if (maxR < minR) maxR = minR;

                std::uniform_int_distribution<int> distC(minC, maxC);
                std::uniform_int_distribution<int> distR(minR, maxR);

                int c = distC(rng);
                int r = distR(rng);
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                out[idx] = static_cast<double>(pointId++);
            }
        }
    }

    void generateStratifiedByRaster(std::vector<double>& out,
                                     const std::vector<double>& strata,
                                     int cols, int rows, int numPoints,
                                     double noData, bool hasND) {
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Group cells by stratum value
        std::map<int, std::vector<int64_t>> strataGroups;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && strata[i] == noData)
                continue;
            int stratum = static_cast<int>(strata[i]);
            strataGroups[stratum].push_back(i);
        }

        if (strataGroups.empty()) return;

        // Count total valid cells for proportional allocation
        int64_t totalValid = 0;
        for (auto& [key, cells] : strataGroups)
            totalValid += static_cast<int64_t>(cells.size());

        std::mt19937_64 rng(std::random_device{}());
        int pointId = 1;

        // Allocate points proportionally to each stratum
        for (auto& [key, cells] : strataGroups) {
            int strataPoints = std::max(1,
                static_cast<int>(std::round(
                    static_cast<double>(numPoints) * cells.size() / totalValid)));

            // Random sample within this stratum
            std::shuffle(cells.begin(), cells.end(), rng);
            int n = std::min(strataPoints, static_cast<int>(cells.size()));
            for (int i = 0; i < n; ++i) {
                out[cells[i]] = static_cast<double>(pointId++);
            }
        }
    }
};

REGISTER_MODULE(SampleModule)

} // namespace aplaceholder
