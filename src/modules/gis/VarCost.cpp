#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <vector>
#include <limits>

namespace aplaceholder {

class VarCostModule : public Module {
public:
    QString name() const override { return "VARCOST"; }
    QString description() const override {
        return "Anisotropic cost distance analysis for phenomena with their own "
               "motive force. Extends COST by moderating frictional effects based "
               "on direction of movement relative to a direction surface.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("source", "Source features image"),
            ParameterDef::file("friction", "Friction surface image"),
            ParameterDef::file("direction", "Direction image (degrees from north)",
                "Azimuth in degrees clockwise from north indicating direction of maximum friction"),
            ParameterDef::real("k_coefficient", "Anisotropy coefficient k",
                1.0, 0.0, 100.0,
                "Higher k = more direction-specific angular response"),
            ParameterDef::output("output", "Output cost distance surface"),
        };
    }

    bool execute() override {
        auto source = GdalIO::read(parameter("source").toString());
        auto friction = GdalIO::read(parameter("friction").toString());
        auto dirRaster = GdalIO::read(parameter("direction").toString());
        if (!source || !friction || !dirRaster) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = source->cols(), rows = source->rows();
        if (friction->cols() != cols || friction->rows() != rows ||
            dirRaster->cols() != cols || dirRaster->rows() != rows) {
            setError("All input rasters must have the same dimensions");
            return false;
        }

        double k = parameter("k_coefficient").toDouble();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& src = source->data(0);
        const auto& fric = friction->data(0);
        const auto& dir = dirRaster->data(0);
        double srcNoData = source->noDataValue();
        bool srcHasND = source->hasNoData();
        double fricNoData = friction->noDataValue();
        bool fricHasND = friction->hasNoData();

        double dx = std::abs(source->geoTransform().pixelWidth);
        double dy = std::abs(source->geoTransform().pixelHeight);
        double dd = std::sqrt(dx * dx + dy * dy);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(source->geoTransform());
        output.setProjection(source->projection());
        output.setNoDataValue(-1);
        auto& dist = output.data(0);

        double inf = std::numeric_limits<double>::max();
        for (int64_t i = 0; i < total; ++i)
            dist[i] = inf;

        using PQEntry = std::pair<double, int64_t>;
        std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

        reportProgress(0.0, "Initializing sources...");
        for (int64_t i = 0; i < total; ++i) {
            bool isSource = (src[i] != 0 && (!srcHasND || src[i] != srcNoData));
            if (isSource) {
                dist[i] = 0.0;
                pq.push({0.0, i});
            }
        }

        if (pq.empty()) {
            setError("No source features found in input");
            return false;
        }

        // 8-directional neighbors: dr, dc, step distance, and movement azimuth (degrees)
        // Movement azimuths: direction of travel from current cell to neighbor
        const int dr[] = {-1, -1, -1,  0, 0,  1,  1, 1};
        const int dc[] = {-1,  0,  1, -1, 1, -1,  0, 1};
        const double stepDist[] = {dd, dy, dd, dx, dx, dd, dy, dd};
        // Azimuths: N=0, NE=45, E=90, etc. Movement to neighbor (row-1=north)
        const double moveAzimuth[] = {315.0, 0.0, 45.0, 270.0, 90.0, 225.0, 180.0, 135.0};

        static constexpr double DEG2RAD = 3.14159265358979323846 / 180.0;

        reportProgress(0.1, "Running anisotropic cost distance...");
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

                double statedFriction = fric[nIdx];
                double effectiveFriction = statedFriction;

                // Apply anisotropic modification if direction is valid (not -1)
                if (dir[nIdx] >= 0) {
                    // alpha = angle between movement direction and max friction direction
                    double alpha = moveAzimuth[d] - dir[nIdx];
                    // Normalize to -180..180
                    while (alpha > 180.0) alpha -= 360.0;
                    while (alpha < -180.0) alpha += 360.0;
                    double alphaRad = alpha * DEG2RAD;

                    // VARCOST formula: effective_friction = stated_friction^f
                    // where f = cos^k(alpha)
                    double cosAlpha = std::cos(alphaRad);
                    double f;
                    if (cosAlpha > 0) {
                        f = std::pow(cosAlpha, k);
                    } else if (cosAlpha < 0) {
                        // Negative cosine: friction becomes force (< 1)
                        f = -std::pow(-cosAlpha, k);
                    } else {
                        f = 0.0;  // At 90 degrees: friction neutralized
                    }

                    // effective_friction = stated_friction^f
                    if (f > 0) {
                        effectiveFriction = std::pow(statedFriction, f);
                    } else if (f < 0) {
                        // Force: effective friction < 1 (1/friction^|f|)
                        effectiveFriction = 1.0 / std::pow(statedFriction, -f);
                    } else {
                        effectiveFriction = 1.0;  // Neutralized
                    }
                }

                double moveCost = effectiveFriction * stepDist[d];
                if (moveCost < 0) moveCost = 0;  // Forces can't produce negative cost
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

REGISTER_MODULE(VarCostModule)

} // namespace aplaceholder
