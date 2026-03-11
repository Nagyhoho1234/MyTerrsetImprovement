#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <queue>

namespace aplaceholder {

class GroupModule : public Module {
public:
    QString name() const override { return "GROUP"; }
    QString description() const override {
        return "Connected component labeling. Assigns unique IDs to connected "
               "groups of non-zero cells (4-connected or 8-connected).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::combo("connectivity", "Connectivity",
                {"4-connected", "8-connected"}, 0,
                "Neighborhood connectivity for grouping"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int connectivity = (parameter("connectivity").toInt() == 0) ? 4 : 8;
        int cols = raster->cols(), rows = raster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(0);

        const auto& src = raster->data(0);
        auto& dst = output.data(0);
        std::fill(dst.begin(), dst.end(), 0.0);

        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        // Direction offsets for 4-connectivity: right, down, left, up
        int dx4[] = {1, 0, -1, 0};
        int dy4[] = {0, 1, 0, -1};
        // Additional diagonals for 8-connectivity
        int dx8[] = {1, 0, -1, 0, 1, -1, 1, -1};
        int dy8[] = {0, 1, 0, -1, 1, 1, -1, -1};

        int* dx = (connectivity == 4) ? dx4 : dx8;
        int* dy = (connectivity == 4) ? dy4 : dy8;
        int numDirs = connectivity;

        int groupId = 0;

        reportProgress(0.0, "Labeling connected components...");

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                int64_t idx = static_cast<int64_t>(row) * cols + col;

                // Skip if already labeled, zero, or nodata
                if (dst[idx] != 0.0)
                    continue;
                if (hasND && src[idx] == noData)
                    continue;
                if (src[idx] == 0.0)
                    continue;

                // BFS flood fill
                ++groupId;
                std::queue<std::pair<int, int>> queue;
                queue.push({col, row});
                dst[idx] = static_cast<double>(groupId);

                while (!queue.empty()) {
                    auto [cx, cy] = queue.front();
                    queue.pop();

                    for (int d = 0; d < numDirs; ++d) {
                        int nx = cx + dx[d];
                        int ny = cy + dy[d];

                        if (nx < 0 || nx >= cols || ny < 0 || ny >= rows)
                            continue;

                        int64_t nIdx = static_cast<int64_t>(ny) * cols + nx;
                        if (dst[nIdx] != 0.0)
                            continue;
                        if (hasND && src[nIdx] == noData)
                            continue;
                        if (src[nIdx] == 0.0)
                            continue;

                        dst[nIdx] = static_cast<double>(groupId);
                        queue.push({nx, ny});
                    }
                }
            }

            if (row % 100 == 0)
                reportProgress(static_cast<double>(row) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(GroupModule)

} // namespace aplaceholder
