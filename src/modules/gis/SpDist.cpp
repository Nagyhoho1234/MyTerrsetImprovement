#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <vector>
#include <limits>

namespace aplaceholder {

class SpDistModule : public Module {
public:
    QString name() const override { return "SPDIST"; }
    QString description() const override {
        return "Shortest path distance on a cost surface using Dijkstra's algorithm. "
               "Outputs the actual path as a binary raster from source to destination.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("cost_raster", "Cost/friction surface"),
            ParameterDef::real("source_x", "Source X coordinate", 0, -1e12, 1e12),
            ParameterDef::real("source_y", "Source Y coordinate", 0, -1e12, 1e12),
            ParameterDef::real("dest_x", "Destination X coordinate", 0, -1e12, 1e12),
            ParameterDef::real("dest_y", "Destination Y coordinate", 0, -1e12, 1e12),
            ParameterDef::output("output", "Output path raster"),
        };
    }

    bool execute() override {
        auto cost = GdalIO::read(parameter("cost_raster").toString());
        if (!cost) { setError("Failed to read cost raster"); return false; }

        int cols = cost->cols(), rows = cost->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& costData = cost->data(0);
        double costND = cost->noDataValue();
        bool costHasND = cost->hasNoData();

        double srcX = parameter("source_x").toDouble();
        double srcY = parameter("source_y").toDouble();
        double dstX = parameter("dest_x").toDouble();
        double dstY = parameter("dest_y").toDouble();

        int srcCol, srcRow, dstCol, dstRow;
        cost->xyToColRow(srcX, srcY, srcCol, srcRow);
        cost->xyToColRow(dstX, dstY, dstCol, dstRow);

        if (srcCol < 0 || srcCol >= cols || srcRow < 0 || srcRow >= rows) {
            setError("Source point is outside the raster extent");
            return false;
        }
        if (dstCol < 0 || dstCol >= cols || dstRow < 0 || dstRow >= rows) {
            setError("Destination point is outside the raster extent");
            return false;
        }

        double dx = std::abs(cost->geoTransform().pixelWidth);
        double dy = std::abs(cost->geoTransform().pixelHeight);
        double dd = std::sqrt(dx * dx + dy * dy);

        const int dr[] = {-1, -1, -1, 0, 0, 1, 1, 1};
        const int dc[] = {-1, 0, 1, -1, 1, -1, 0, 1};
        const double stepDist[] = {dd, dy, dd, dx, dx, dd, dy, dd};

        double inf = std::numeric_limits<double>::max();
        std::vector<double> dist(total, inf);
        std::vector<int64_t> prev(total, -1);

        // Dijkstra from source
        reportProgress(0.0, "Running Dijkstra shortest path...");
        using PQEntry = std::pair<double, int64_t>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

        int64_t srcIdx = static_cast<int64_t>(srcRow) * cols + srcCol;
        int64_t dstIdx = static_cast<int64_t>(dstRow) * cols + dstCol;
        dist[srcIdx] = 0.0;
        pq.push({0.0, srcIdx});

        int64_t processed = 0;
        bool found = false;

        while (!pq.empty()) {
            auto [curCost, idx] = pq.top();
            pq.pop();

            if (curCost > dist[idx]) continue;

            if (idx == dstIdx) {
                found = true;
                break;
            }

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

                if (newCost < dist[nIdx]) {
                    dist[nIdx] = newCost;
                    prev[nIdx] = idx;
                    pq.push({newCost, nIdx});
                }
            }

            processed++;
            if (processed % 500000 == 0)
                reportProgress(0.1 + 0.7 * static_cast<double>(processed) / total);
        }

        if (!found) {
            setError("No path found between source and destination");
            return false;
        }

        // Trace back path
        reportProgress(0.85, "Tracing path...");
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(cost->geoTransform());
        output.setProjection(cost->projection());
        output.setNoDataValue(0);
        auto& out = output.data(0);

        int64_t cur = dstIdx;
        while (cur >= 0) {
            out[cur] = 1.0;
            cur = prev[cur];
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SpDistModule)

} // namespace aplaceholder
