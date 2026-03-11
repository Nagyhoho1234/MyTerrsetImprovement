#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <queue>

namespace aplaceholder {

class HinterlandModule : public Module {
public:
    QString name() const override { return "HINTERLAND"; }
    QString description() const override {
        return "Compute hinterland (service area) regions around features. "
               "Each non-feature cell is assigned to the nearest feature based "
               "on Euclidean distance, producing a Voronoi-like allocation map.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input feature raster (non-zero = features)"),
            ParameterDef::output("output", "Output hinterland allocation image"),
            ParameterDef::real("max_distance", "Maximum distance (cells)", 0.0, 0.0, 999999.0,
                "Maximum allocation distance in cells (0 = unlimited)"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        double maxDist = parameter("max_distance").toDouble();

        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        // Distance grid
        std::vector<double> dist(total, 1e30);

        // Multi-source BFS from all feature cells
        struct Cell { int r, c; };
        std::queue<Cell> q;

        reportProgress(0.0, "Initializing features...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                if (hasND && data[idx] == noData) {
                    out[idx] = noData;
                    continue;
                }
                if (data[idx] != 0.0) {
                    // Feature cell
                    out[idx] = data[idx];
                    dist[idx] = 0.0;
                    q.push({r, c});
                } else {
                    out[idx] = 0.0;
                }
            }
        }

        reportProgress(0.1, "Computing hinterlands...");

        int dx[] = {1, -1, 0, 0, 1, 1, -1, -1};
        int dy[] = {0, 0, 1, -1, 1, -1, 1, -1};
        double dd[] = {1.0, 1.0, 1.0, 1.0, 1.414213562, 1.414213562, 1.414213562, 1.414213562};

        // Iterative relaxation (approximate Euclidean via chamfer distance)
        while (!q.empty()) {
            Cell cur = q.front(); q.pop();
            int64_t curIdx = static_cast<int64_t>(cur.r) * cols + cur.c;
            double curDist = dist[curIdx];

            for (int d = 0; d < 8; ++d) {
                int nr = cur.r + dy[d], nc = cur.c + dx[d];
                if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                if (hasND && data[nIdx] == noData) continue;

                double newDist = curDist + dd[d];
                if (maxDist > 0.0 && newDist > maxDist) continue;

                if (newDist < dist[nIdx]) {
                    dist[nIdx] = newDist;
                    out[nIdx] = out[curIdx];
                    q.push({nr, nc});
                }
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(HinterlandModule)

} // namespace aplaceholder
