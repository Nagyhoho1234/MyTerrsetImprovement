#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <vector>
#include <limits>

namespace aplaceholder {

class DisperseModule : public Module {
public:
    QString name() const override { return "DISPERSE"; }
    QString description() const override {
        return "Dispersal/spread modeling for phenomena with no motive force. "
               "Simulates spread from source locations over multiple time steps "
               "with friction and distance decay, producing plume-like patterns.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("source", "Source features image"),
            ParameterDef::file("friction", "Friction/force magnitude image"),
            ParameterDef::output("output", "Output dispersion surface"),
            ParameterDef::integer("num_steps", "Number of time steps",
                10, 1, 10000,
                "Number of iterative spread steps"),
            ParameterDef::real("decay_rate", "Distance decay rate",
                0.1, 0.0, 10.0,
                "Rate at which spread intensity decays with distance"),
        };
    }

    bool execute() override {
        auto source = GdalIO::read(parameter("source").toString());
        auto friction = GdalIO::read(parameter("friction").toString());
        if (!source || !friction) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = source->cols(), rows = source->rows();
        if (friction->cols() != cols || friction->rows() != rows) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int numSteps = parameter("num_steps").toInt();
        double decayRate = parameter("decay_rate").toDouble();
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

        // Intensity surface: tracks accumulated spread intensity at each cell
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(source->geoTransform());
        output.setProjection(source->projection());
        output.setNoDataValue(-1);
        auto& intensity = output.data(0);

        for (int64_t i = 0; i < total; ++i)
            intensity[i] = 0.0;

        // Current wavefront: cells that can spread in this step
        // Each entry: (intensity_at_cell, linear_index)
        std::vector<double> current(total, 0.0);
        std::vector<double> next(total, 0.0);

        // Initialize sources with intensity 1.0
        reportProgress(0.0, "Initializing sources...");
        for (int64_t i = 0; i < total; ++i) {
            bool isSource = (src[i] != 0 && (!srcHasND || src[i] != srcNoData));
            if (isSource) {
                current[i] = 1.0;
                intensity[i] = 1.0;
            }
        }

        // 8-directional neighbors
        const int dr[] = {-1, -1, -1, 0, 0, 1, 1, 1};
        const int dc[] = {-1, 0, 1, -1, 1, -1, 0, 1};
        const double stepDist[] = {dd, dy, dd, dx, dx, dd, dy, dd};

        // Iterative spread over time steps
        for (int step = 0; step < numSteps; ++step) {
            reportProgress(0.05 + 0.9 * static_cast<double>(step) / numSteps,
                           QString("Step %1 of %2").arg(step + 1).arg(numSteps));

            for (int64_t i = 0; i < total; ++i)
                next[i] = 0.0;

            bool anySpread = false;

            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    int64_t idx = static_cast<int64_t>(r) * cols + c;
                    if (current[idx] <= 0.0)
                        continue;

                    for (int d = 0; d < 8; ++d) {
                        int nr = r + dr[d];
                        int nc = c + dc[d];
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
                            continue;

                        int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;

                        // Skip nodata/impassable friction
                        if (fricHasND && fric[nIdx] == fricNoData)
                            continue;
                        if (fric[nIdx] <= 0)
                            continue;

                        // Spread intensity decays with distance and friction
                        double distFactor = stepDist[d] / dx;  // Normalize by cell size
                        double decay = std::exp(-decayRate * distFactor);
                        double spreadIntensity = current[idx] * decay / fric[nIdx];

                        if (spreadIntensity > next[nIdx]) {
                            next[nIdx] = spreadIntensity;
                            anySpread = true;
                        }
                    }
                }
            }

            // Update intensity: accumulate maximum spread intensity
            for (int64_t i = 0; i < total; ++i) {
                if (next[i] > intensity[i])
                    intensity[i] = next[i];
            }

            // Swap current and next for next iteration
            std::swap(current, next);

            if (!anySpread) break;  // No more spread possible
        }

        // Mark cells never reached as nodata
        for (int64_t i = 0; i < total; ++i) {
            if (intensity[i] <= 0.0) {
                // Check if the cell is impassable
                if (fricHasND && fric[i] == fricNoData)
                    intensity[i] = -1;
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(DisperseModule)

} // namespace aplaceholder
