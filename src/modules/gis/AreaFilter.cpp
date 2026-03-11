#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <queue>
#include <map>

namespace aplaceholder {

class AreaFilterModule : public Module {
public:
    QString name() const override { return "AREAFILTER"; }
    QString description() const override {
        return "Remove small contiguous groups of pixels below a size threshold. "
               "Uses flood-fill to identify connected regions and replaces groups "
               "smaller than the threshold with the most common neighboring value.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input categorical raster image"),
            ParameterDef::output("output", "Output filtered image"),
            ParameterDef::integer("threshold", "Minimum area (pixels)", 10, 1, 999999,
                "Groups with fewer pixels than this will be removed"),
            ParameterDef::combo("connectivity", "Connectivity",
                {"4-connected", "8-connected"}, 1,
                "Pixel connectivity for defining groups"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        int threshold = parameter("threshold").toInt();
        int connectivity = parameter("connectivity").toInt();
        bool use8 = (connectivity == 1);

        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Copy data to output
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);
        for (int64_t i = 0; i < total; ++i) out[i] = data[i];

        // Label connected components via flood fill
        std::vector<int> labels(total, -1);
        int nextLabel = 0;

        // Direction offsets
        int dx4[] = {1, -1, 0, 0};
        int dy4[] = {0, 0, 1, -1};
        int dx8[] = {1, -1, 0, 0, 1, 1, -1, -1};
        int dy8[] = {0, 0, 1, -1, 1, -1, 1, -1};
        int numDirs = use8 ? 8 : 4;
        int* dx = use8 ? dx8 : dx4;
        int* dy = use8 ? dy8 : dy4;

        reportProgress(0.0, "Labeling connected regions...");

        struct GroupInfo {
            double value;
            int64_t count;
            std::vector<int64_t> cells;
        };
        std::vector<GroupInfo> groups;

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                if (labels[idx] >= 0) continue;
                if (hasND && data[idx] == noData) { labels[idx] = -2; continue; }

                // BFS flood fill
                int label = nextLabel++;
                double val = data[idx];
                GroupInfo gi;
                gi.value = val;
                gi.count = 0;

                std::queue<int64_t> q;
                q.push(idx);
                labels[idx] = label;

                while (!q.empty()) {
                    int64_t cur = q.front(); q.pop();
                    gi.cells.push_back(cur);
                    gi.count++;
                    int cr = static_cast<int>(cur / cols);
                    int cc = static_cast<int>(cur % cols);

                    for (int d = 0; d < numDirs; ++d) {
                        int nr = cr + dy[d], nc = cc + dx[d];
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                        int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                        if (labels[nIdx] >= 0 || labels[nIdx] == -2) continue;
                        if (hasND && data[nIdx] == noData) { labels[nIdx] = -2; continue; }
                        if (data[nIdx] != val) continue;
                        labels[nIdx] = label;
                        q.push(nIdx);
                    }
                }
                groups.push_back(std::move(gi));
            }
            if (r % 200 == 0) reportProgress(0.3 * static_cast<double>(r) / rows);
        }

        reportProgress(0.5, "Removing small groups...");

        // For small groups, replace with most common neighbor value
        for (size_t g = 0; g < groups.size(); ++g) {
            if (groups[g].count >= threshold) continue;

            // Find most common neighboring value
            std::map<double, int> neighborFreq;
            for (int64_t cell : groups[g].cells) {
                int cr = static_cast<int>(cell / cols);
                int cc = static_cast<int>(cell % cols);
                for (int d = 0; d < numDirs; ++d) {
                    int nr = cr + dy[d], nc = cc + dx[d];
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                    int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                    if (hasND && data[nIdx] == noData) continue;
                    if (data[nIdx] == groups[g].value) continue;
                    neighborFreq[data[nIdx]]++;
                }
            }

            double replaceVal = noData;
            int maxFreq = 0;
            for (const auto& kv : neighborFreq) {
                if (kv.second > maxFreq) {
                    maxFreq = kv.second;
                    replaceVal = kv.first;
                }
            }

            if (maxFreq > 0) {
                for (int64_t cell : groups[g].cells) {
                    out[cell] = replaceVal;
                }
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(AreaFilterModule)

} // namespace aplaceholder
