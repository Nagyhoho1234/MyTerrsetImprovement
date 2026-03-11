#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <vector>
#include <limits>

namespace aplaceholder {

class CostGrowModule : public Module {
public:
    QString name() const override { return "COSTGROW"; }
    QString description() const override {
        return "Cost surface growth/allocation. Grows regions from seed cells outward "
               "based on a cost surface, assigning each cell to the cheapest seed.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("seed_raster", "Seed raster (regions with IDs)"),
            ParameterDef::file("cost_raster", "Cost/friction surface"),
            ParameterDef::output("output", "Output allocation raster"),
            ParameterDef::real("max_cost", "Maximum cost (0 = unlimited)", 0, 0, 1e12,
                "Maximum accumulated cost for growth"),
        };
    }

    bool execute() override {
        auto seeds = GdalIO::read(parameter("seed_raster").toString());
        auto cost = GdalIO::read(parameter("cost_raster").toString());
        if (!seeds || !cost) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = seeds->cols(), rows = seeds->rows();
        if (cost->cols() != cols || cost->rows() != rows) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& seedData = seeds->data(0);
        const auto& costData = cost->data(0);
        double seedND = seeds->noDataValue();
        bool seedHasND = seeds->hasNoData();
        double costND = cost->noDataValue();
        bool costHasND = cost->hasNoData();
        double maxCost = parameter("max_cost").toDouble();
        if (maxCost <= 0) maxCost = std::numeric_limits<double>::max();

        double dx = std::abs(seeds->geoTransform().pixelWidth);
        double dy = std::abs(seeds->geoTransform().pixelHeight);
        double dd = std::sqrt(dx * dx + dy * dy);

        const int dr[] = {-1, -1, -1, 0, 0, 1, 1, 1};
        const int dc[] = {-1, 0, 1, -1, 1, -1, 0, 1};
        const double stepDist[] = {dd, dy, dd, dx, dx, dd, dy, dd};

        // Output allocation raster
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(seeds->geoTransform());
        output.setProjection(seeds->projection());
        output.setNoDataValue(0);
        auto& out = output.data(0);

        double inf = std::numeric_limits<double>::max();
        std::vector<double> dist(total, inf);

        // Priority queue: (accumulated_cost, linear_index)
        using PQEntry = std::pair<double, int64_t>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

        // Initialize seeds
        reportProgress(0.0, "Initializing seeds...");
        for (int64_t i = 0; i < total; ++i) {
            if (seedHasND && seedData[i] == seedND) continue;
            if (seedData[i] == 0) continue;

            dist[i] = 0.0;
            out[i] = seedData[i];
            pq.push({0.0, i});
        }

        if (pq.empty()) {
            setError("No seed cells found");
            return false;
        }

        // Dijkstra growth
        reportProgress(0.1, "Growing regions...");
        int64_t processed = 0;

        while (!pq.empty()) {
            auto [curCost, idx] = pq.top();
            pq.pop();

            if (curCost > dist[idx]) continue;

            int r = static_cast<int>(idx / cols);
            int c = static_cast<int>(idx % cols);

            for (int d = 0; d < 8; ++d) {
                int nr = r + dr[d];
                int nc = c + dc[d];
                if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;

                int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                if (costHasND && costData[nIdx] == costND) continue;
                if (costData[nIdx] <= 0) continue;

                double moveCost = costData[nIdx] * stepDist[d];
                double newCost = curCost + moveCost;

                if (newCost > maxCost) continue;

                if (newCost < dist[nIdx]) {
                    dist[nIdx] = newCost;
                    out[nIdx] = out[idx]; // assign same region as source
                    pq.push({newCost, nIdx});
                }
            }

            processed++;
            if (processed % 500000 == 0)
                reportProgress(0.1 + 0.85 * static_cast<double>(processed) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(CostGrowModule)

} // namespace aplaceholder
