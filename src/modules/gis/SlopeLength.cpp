#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <vector>

namespace aplaceholder {

class SlopeLengthModule : public Module {
public:
    QString name() const override { return "SLOPELENGTH"; }
    QString description() const override {
        return "Computes the LS factor for RUSLE from a DEM using slope and flow accumulation. "
               "LS = (flow_accum * cellsize / 22.13)^0.4 * (sin(slope) / 0.0896)^1.3.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dem", "Input DEM"),
            ParameterDef::output("output", "Output LS factor raster"),
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
        double cellsize = (dx + dy) / 2.0;

        // D8 directions
        const int dr[] = {0, 1, 1, 1, 0, -1, -1, -1};
        const int dc[] = {1, 1, 0, -1, -1, -1, 0, 1};
        const double stepDist[] = {dx, dd, dy, dd, dx, dd, dy, dd};

        // Step 1: Compute flow direction and accumulation
        reportProgress(0.0, "Computing flow directions...");
        std::vector<int> flowDir(total, -1);
        std::vector<int> inDegree(total, 0);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                if (hasND && elev[idx] == noData) continue;

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
                    flowDir[idx] = bestDir;
                    int nr = r + dr[bestDir];
                    int nc = c + dc[bestDir];
                    int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                    inDegree[nIdx]++;
                }
            }
        }

        // Topological sort for flow accumulation
        reportProgress(0.2, "Computing flow accumulation...");
        std::vector<double> accum(total, 1.0);
        std::queue<int64_t> q;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && elev[i] == noData) { accum[i] = 0; continue; }
            if (inDegree[i] == 0)
                q.push(i);
        }

        while (!q.empty()) {
            int64_t idx = q.front();
            q.pop();
            int d = flowDir[idx];
            if (d < 0) continue;

            int r = static_cast<int>(idx / cols);
            int c = static_cast<int>(idx % cols);
            int nr = r + dr[d];
            int nc = c + dc[d];
            if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                accum[nIdx] += accum[idx];
                inDegree[nIdx]--;
                if (inDegree[nIdx] == 0)
                    q.push(nIdx);
            }
        }

        // Step 2: Compute slope using Horn's method and LS factor
        reportProgress(0.5, "Computing LS factor...");
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(dem->geoTransform());
        output.setProjection(dem->projection());
        output.setNoDataValue(-9999);
        auto& out = output.data(0);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;

                if (hasND && elev[idx] == noData) {
                    out[idx] = -9999;
                    continue;
                }

                // Boundary cells: use simple slope estimate
                if (r == 0 || r == rows - 1 || c == 0 || c == cols - 1) {
                    out[idx] = 0.0;
                    continue;
                }

                // Horn's method for slope
                double z1 = elev[(r-1)*cols+(c-1)];
                double z2 = elev[(r-1)*cols+c];
                double z3 = elev[(r-1)*cols+(c+1)];
                double z4 = elev[r*cols+(c-1)];
                double z6 = elev[r*cols+(c+1)];
                double z7 = elev[(r+1)*cols+(c-1)];
                double z8 = elev[(r+1)*cols+c];
                double z9 = elev[(r+1)*cols+(c+1)];

                double dzdx = ((z3 + 2*z6 + z9) - (z1 + 2*z4 + z7)) / (8.0 * dx);
                double dzdy = ((z7 + 2*z8 + z9) - (z1 + 2*z2 + z3)) / (8.0 * dy);
                double slopeRad = std::atan(std::sqrt(dzdx*dzdx + dzdy*dzdy));
                double sinSlope = std::sin(slopeRad);

                double flowLength = accum[idx] * cellsize / 22.13;
                double slopePart = sinSlope / 0.0896;

                if (flowLength <= 0) flowLength = 1e-10;
                if (slopePart <= 0) slopePart = 1e-10;

                out[idx] = std::pow(flowLength, 0.4) * std::pow(slopePart, 1.3);
            }
            if (r % 100 == 0)
                reportProgress(0.5 + 0.45 * static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SlopeLengthModule)

} // namespace aplaceholder
