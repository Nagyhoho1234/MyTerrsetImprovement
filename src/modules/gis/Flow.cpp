#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <queue>

namespace aplaceholder {

class FlowModule : public Module {
public:
    QString name() const override { return "FLOW"; }
    QString description() const override {
        return "Computes D8 flow direction and flow accumulation from a DEM.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dem", "Input DEM"),
            ParameterDef::output("output_direction", "Output flow direction raster"),
            ParameterDef::output("output_accumulation", "Output flow accumulation raster"),
        };
    }

    bool execute() override {
        auto dem = GdalIO::read(parameter("dem").toString());
        if (!dem) { setError("Failed to read DEM"); return false; }

        int cols = dem->cols(), rows = dem->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& elev = dem->data(0);
        double noData = dem->noDataValue();
        bool hasND = dem->hasNoData();

        double dx = std::abs(dem->geoTransform().pixelWidth);
        double dy = std::abs(dem->geoTransform().pixelHeight);
        double dd = std::sqrt(dx * dx + dy * dy);

        // D8 direction encoding: 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE
        //                  dr: -1, -1, -1,  0, 0,  1, 1, 1
        //                  dc: -1,  0,  1, -1, 1, -1, 0, 1
        const int dr[] = {0, 1, 1, 1, 0, -1, -1, -1};
        const int dc[] = {1, 1, 0, -1, -1, -1, 0, 1};
        const int dirCode[] = {1, 2, 4, 8, 16, 32, 64, 128};
        const double stepDist[] = {dx, dd, dy, dd, dx, dd, dy, dd};

        // Flow direction
        Raster dirRaster(cols, rows, 1, DataType::Float64);
        dirRaster.setGeoTransform(dem->geoTransform());
        dirRaster.setProjection(dem->projection());
        dirRaster.setNoDataValue(-1);
        auto& dirData = dirRaster.data(0);

        // Flow accumulation
        Raster accRaster(cols, rows, 1, DataType::Float64);
        accRaster.setGeoTransform(dem->geoTransform());
        accRaster.setProjection(dem->projection());
        accRaster.setNoDataValue(-1);
        auto& accData = accRaster.data(0);

        // Compute flow directions using D8
        reportProgress(0.0, "Computing flow directions...");
        std::vector<int> inDegree(total, 0);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;

                if (hasND && elev[idx] == noData) {
                    dirData[idx] = -1;
                    continue;
                }

                double maxDrop = -1.0;
                int bestDir = -1;

                for (int d = 0; d < 8; ++d) {
                    int nr = r + dr[d];
                    int nc = c + dc[d];
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;

                    int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                    if (hasND && elev[nIdx] == noData) continue;

                    double drop = (elev[idx] - elev[nIdx]) / stepDist[d];
                    if (drop > maxDrop) {
                        maxDrop = drop;
                        bestDir = d;
                    }
                }

                if (bestDir >= 0 && maxDrop > 0) {
                    dirData[idx] = dirCode[bestDir];
                    int nr = r + dr[bestDir];
                    int nc = c + dc[bestDir];
                    int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                    inDegree[nIdx]++;
                } else {
                    dirData[idx] = 0; // pit or flat
                }
            }
            if (r % 100 == 0)
                reportProgress(0.4 * static_cast<double>(r) / rows);
        }

        // Compute flow accumulation using topological sort
        reportProgress(0.4, "Computing flow accumulation...");
        for (int64_t i = 0; i < total; ++i)
            accData[i] = 1.0; // each cell counts itself

        std::queue<int64_t> q;
        for (int64_t i = 0; i < total; ++i) {
            if (inDegree[i] == 0 && !(hasND && elev[i] == noData))
                q.push(i);
        }

        int64_t processed = 0;
        while (!q.empty()) {
            int64_t idx = q.front();
            q.pop();

            int dirVal = static_cast<int>(dirData[idx]);
            if (dirVal <= 0) { processed++; continue; }

            // Find which direction index this code corresponds to
            int bestDir = -1;
            for (int d = 0; d < 8; ++d) {
                if (dirCode[d] == dirVal) { bestDir = d; break; }
            }
            if (bestDir < 0) { processed++; continue; }

            int r = static_cast<int>(idx / cols);
            int c = static_cast<int>(idx % cols);
            int nr = r + dr[bestDir];
            int nc = c + dc[bestDir];
            if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                accData[nIdx] += accData[idx];
                inDegree[nIdx]--;
                if (inDegree[nIdx] == 0)
                    q.push(nIdx);
            }

            processed++;
            if (processed % 500000 == 0)
                reportProgress(0.4 + 0.55 * static_cast<double>(processed) / total);
        }

        // Set nodata cells
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && elev[i] == noData)
                accData[i] = -1;
        }

        reportProgress(0.95, "Writing outputs...");
        if (!GdalIO::write(dirRaster, parameter("output_direction").toString())) {
            setError("Failed to write flow direction raster");
            return false;
        }
        if (!GdalIO::write(accRaster, parameter("output_accumulation").toString())) {
            setError("Failed to write flow accumulation raster");
            return false;
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(FlowModule)

} // namespace aplaceholder
