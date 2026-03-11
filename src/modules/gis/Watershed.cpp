#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <vector>

namespace aplaceholder {

class WatershedModule : public Module {
public:
    QString name() const override { return "WATERSHED"; }
    QString description() const override {
        return "Delineates watersheds from a flow direction raster and pour points. "
               "Traces upstream from each pour point to find all contributing cells.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("flow_direction", "Flow direction raster (D8)"),
            ParameterDef::file("pour_points", "Pour points raster (point IDs)"),
            ParameterDef::output("output", "Output watershed raster"),
        };
    }

    bool execute() override {
        auto flowDir = GdalIO::read(parameter("flow_direction").toString());
        auto pourPts = GdalIO::read(parameter("pour_points").toString());
        if (!flowDir || !pourPts) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = flowDir->cols(), rows = flowDir->rows();
        if (pourPts->cols() != cols || pourPts->rows() != rows) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& dirData = flowDir->data(0);
        const auto& ppData = pourPts->data(0);
        double ppND = pourPts->noDataValue();
        bool ppHasND = pourPts->hasNoData();

        // D8 direction codes and offsets
        // Code:  1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE
        const int dr[] = {0, 1, 1, 1, 0, -1, -1, -1};
        const int dc[] = {1, 1, 0, -1, -1, -1, 0, 1};
        const int dirCode[] = {1, 2, 4, 8, 16, 32, 64, 128};
        // Opposite directions: if neighbor flows toward us with code X,
        // then neighbor's direction code is the opposite
        const int oppCode[] = {16, 32, 64, 128, 1, 2, 4, 8};

        // Build reverse flow graph: for each cell, find which neighbors flow INTO it
        reportProgress(0.0, "Building flow network...");
        std::vector<std::vector<int64_t>> inflowCells(total);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                int code = static_cast<int>(dirData[idx]);
                if (code <= 0) continue;

                for (int d = 0; d < 8; ++d) {
                    if (dirCode[d] == code) {
                        int nr = r + dr[d];
                        int nc = c + dc[d];
                        if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                            int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                            inflowCells[nIdx].push_back(idx);
                        }
                        break;
                    }
                }
            }
            if (r % 200 == 0)
                reportProgress(0.3 * static_cast<double>(r) / rows);
        }

        // Output raster
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(flowDir->geoTransform());
        output.setProjection(flowDir->projection());
        output.setNoDataValue(0);
        auto& out = output.data(0);

        // BFS upstream from each pour point
        reportProgress(0.3, "Delineating watersheds...");
        std::queue<int64_t> q;

        for (int64_t i = 0; i < total; ++i) {
            if (ppHasND && ppData[i] == ppND) continue;
            if (ppData[i] == 0) continue;

            double wsID = ppData[i];
            q.push(i);
            out[i] = wsID;

            while (!q.empty()) {
                int64_t cur = q.front();
                q.pop();

                // Visit all cells that flow into this cell
                for (int64_t upstream : inflowCells[cur]) {
                    if (out[upstream] == 0) {
                        out[upstream] = wsID;
                        q.push(upstream);
                    }
                }
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(WatershedModule)

} // namespace aplaceholder
