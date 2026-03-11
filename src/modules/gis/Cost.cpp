#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <vector>
#include <limits>

namespace aplaceholder {

class CostModule : public Module {
public:
    QString name() const override { return "COST"; }
    QString description() const override {
        return "Cost distance surface over friction layers. "
               "Computes the accumulated cost of moving from source features "
               "across a friction surface using Dijkstra's algorithm.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_features", "Input source features"),
            ParameterDef::file("friction_surface", "Friction surface image"),
            ParameterDef::output("output", "Output cost distance surface"),
        };
    }

    bool execute() override {
        auto features = GdalIO::read(parameter("input_features").toString());
        auto friction = GdalIO::read(parameter("friction_surface").toString());
        if (!features || !friction) {
            setError("Failed to read input rasters");
            return false;
        }

        if (features->cols() != friction->cols() || features->rows() != friction->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = features->cols(), rows = features->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& feat = features->data(0);
        const auto& fric = friction->data(0);
        double featNoData = features->noDataValue();
        bool featHasND = features->hasNoData();
        double fricNoData = friction->noDataValue();
        bool fricHasND = friction->hasNoData();

        double dx = std::abs(features->geoTransform().pixelWidth);
        double dy = std::abs(features->geoTransform().pixelHeight);
        double dd = std::sqrt(dx * dx + dy * dy);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(features->geoTransform());
        output.setProjection(features->projection());
        output.setNoDataValue(-1);
        auto& dist = output.data(0);

        double inf = std::numeric_limits<double>::max();
        for (int64_t i = 0; i < total; ++i)
            dist[i] = inf;

        // Priority queue: (accumulated_cost, linear_index)
        using PQEntry = std::pair<double, int64_t>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

        // Initialize: all target cells (non-zero, non-nodata) have cost 0
        reportProgress(0.0, "Initializing sources...");
        for (int64_t i = 0; i < total; ++i) {
            bool isTarget = (feat[i] != 0 && (!featHasND || feat[i] != featNoData));
            if (isTarget) {
                dist[i] = 0.0;
                pq.push({0.0, i});
            }
        }

        if (pq.empty()) {
            setError("No target features found in input");
            return false;
        }

        // 8-directional neighbors: (row_offset, col_offset, distance_multiplier)
        const int dr[] = {-1, -1, -1, 0, 0, 1, 1, 1};
        const int dc[] = {-1, 0, 1, -1, 1, -1, 0, 1};
        const double stepDist[] = {dd, dy, dd, dx, dx, dd, dy, dd};

        // Dijkstra's algorithm
        reportProgress(0.1, "Running Dijkstra cost distance...");
        int64_t processed = 0;
        int64_t totalSources = static_cast<int64_t>(pq.size());

        while (!pq.empty()) {
            auto [cost, idx] = pq.top();
            pq.pop();

            // Skip if we already found a shorter path
            if (cost > dist[idx]) continue;

            int r = static_cast<int>(idx / cols);
            int c = static_cast<int>(idx % cols);

            for (int d = 0; d < 8; ++d) {
                int nr = r + dr[d];
                int nc = c + dc[d];
                if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
                    continue;

                int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;

                // Skip nodata in friction surface
                if (fricHasND && fric[nIdx] == fricNoData)
                    continue;

                // Skip negative or zero friction (impassable)
                if (fric[nIdx] <= 0)
                    continue;

                // Cost = friction * distance to neighbor
                double moveCost = fric[nIdx] * stepDist[d];
                double newCost = cost + moveCost;

                if (newCost < dist[nIdx]) {
                    dist[nIdx] = newCost;
                    pq.push({newCost, nIdx});
                }
            }

            processed++;
            if (processed % 500000 == 0)
                reportProgress(0.1 + 0.85 * static_cast<double>(processed) / total);
        }

        // Set unreachable cells to nodata
        for (int64_t i = 0; i < total; ++i) {
            if (dist[i] == inf)
                dist[i] = -1;
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(CostModule)

} // namespace aplaceholder
