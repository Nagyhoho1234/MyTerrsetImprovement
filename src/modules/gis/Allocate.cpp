#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <vector>
#include <limits>

namespace aplaceholder {

class AllocateModule : public Module {
public:
    QString name() const override { return "ALLOCATE"; }
    QString description() const override {
        return "Zone allocation by cost distance. Assigns each pixel to its "
               "nearest source feature based on accumulated cost through a "
               "friction surface, producing cost-based Thiessen-like polygons.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("source", "Source features image"),
            ParameterDef::file("friction", "Friction surface image"),
            ParameterDef::output("output_allocation", "Output allocation image"),
            ParameterDef::output("output_cost", "Output cost distance image"),
        };
    }

    bool execute() override {
        auto source = GdalIO::read(parameter("source").toString());
        auto friction = GdalIO::read(parameter("friction").toString());
        if (!source || !friction) {
            setError("Failed to read input rasters");
            return false;
        }

        if (source->cols() != friction->cols() || source->rows() != friction->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = source->cols(), rows = source->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& src = source->data(0);
        const auto& fric = friction->data(0);
        double srcNoData = source->noDataValue();
        bool srcHasND = source->hasNoData();
        double fricNoData = friction->noDataValue();
        bool fricHasND = friction->hasNoData();

        double dx = std::abs(source->geoTransform().pixelWidth);
        double dy = std::abs(source->geoTransform().pixelHeight);
        double dd = std::sqrt(dx * dx + dy * dy);

        // Cost distance output
        Raster costOut(cols, rows, 1, DataType::Float64);
        costOut.setGeoTransform(source->geoTransform());
        costOut.setProjection(source->projection());
        costOut.setNoDataValue(-1);
        auto& dist = costOut.data(0);

        // Allocation output
        Raster allocOut(cols, rows, 1, DataType::Float64);
        allocOut.setGeoTransform(source->geoTransform());
        allocOut.setProjection(source->projection());
        allocOut.setNoDataValue(0);
        auto& alloc = allocOut.data(0);

        double inf = std::numeric_limits<double>::max();
        for (int64_t i = 0; i < total; ++i) {
            dist[i] = inf;
            alloc[i] = 0;
        }

        // Priority queue: (accumulated_cost, linear_index)
        using PQEntry = std::pair<double, int64_t>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

        // Initialize: all source cells (non-zero, non-nodata) have cost 0
        reportProgress(0.0, "Initializing sources...");
        for (int64_t i = 0; i < total; ++i) {
            bool isSource = (src[i] != 0 && (!srcHasND || src[i] != srcNoData));
            if (isSource) {
                dist[i] = 0.0;
                alloc[i] = src[i];  // Track which source class this cell belongs to
                pq.push({0.0, i});
            }
        }

        if (pq.empty()) {
            setError("No source features found in input");
            return false;
        }

        // 8-directional neighbors
        const int dr[] = {-1, -1, -1, 0, 0, 1, 1, 1};
        const int dc[] = {-1, 0, 1, -1, 1, -1, 0, 1};
        const double stepDist[] = {dd, dy, dd, dx, dx, dd, dy, dd};

        // Dijkstra's algorithm tracking allocation
        reportProgress(0.1, "Running allocation by cost distance...");
        int64_t processed = 0;

        while (!pq.empty()) {
            auto [cost, idx] = pq.top();
            pq.pop();

            if (cost > dist[idx]) continue;

            int r = static_cast<int>(idx / cols);
            int c = static_cast<int>(idx % cols);

            for (int d = 0; d < 8; ++d) {
                int nr = r + dr[d];
                int nc = c + dc[d];
                if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
                    continue;

                int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;

                if (fricHasND && fric[nIdx] == fricNoData)
                    continue;

                if (fric[nIdx] <= 0)
                    continue;

                double moveCost = fric[nIdx] * stepDist[d];
                double newCost = cost + moveCost;

                if (newCost < dist[nIdx]) {
                    dist[nIdx] = newCost;
                    alloc[nIdx] = alloc[idx];  // Inherit source ID from parent
                    pq.push({newCost, nIdx});
                }
            }

            processed++;
            if (processed % 500000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(processed) / total);
        }

        // Set unreachable cells to nodata
        for (int64_t i = 0; i < total; ++i) {
            if (dist[i] == inf) {
                dist[i] = -1;
                alloc[i] = 0;
            }
        }

        reportProgress(0.95, "Writing outputs...");
        if (!GdalIO::write(allocOut, parameter("output_allocation").toString())) {
            setError("Failed to write allocation output");
            return false;
        }
        if (!GdalIO::write(costOut, parameter("output_cost").toString())) {
            setError("Failed to write cost distance output");
            return false;
        }

        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(AllocateModule)

} // namespace aplaceholder
