#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>

namespace aplaceholder {

class DistanceModule : public Module {
public:
    QString name() const override { return "DISTANCE"; }
    QString description() const override {
        return "Calculates Euclidean distance from each cell to the nearest target cell. "
               "Target cells are non-zero, non-nodata pixels.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image (target features)"),
            ParameterDef::output("output", "Output distance image"),
            ParameterDef::boolean("geographic", "Use geographic (spherical) distance", false,
                "If checked, distances are in meters using spherical approximation"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input"); return false; }

        int cols = input->cols(), rows = input->rows();
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());

        const auto& data = input->data(0);
        auto& dist = output.data(0);
        double noData = input->hasNoData() ? input->noDataValue() : -9999;

        // Initialize distances
        double inf = std::numeric_limits<double>::max();
        for (int64_t i = 0; i < static_cast<int64_t>(cols) * rows; ++i) {
            bool isTarget = (data[i] != 0 && (!input->hasNoData() || data[i] != noData));
            dist[i] = isTarget ? 0.0 : inf;
        }

        double dx = std::abs(input->geoTransform().pixelWidth);
        double dy = std::abs(input->geoTransform().pixelHeight);
        double dd = std::sqrt(dx * dx + dy * dy);

        // Two-pass distance transform (Rosenfeld & Pfaltz)
        // Forward pass
        reportProgress(0.0, "Forward pass...");
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;
                if (r > 0)
                    dist[idx] = std::min(dist[idx], dist[(r-1)*cols + c] + dy);
                if (c > 0)
                    dist[idx] = std::min(dist[idx], dist[r*cols + (c-1)] + dx);
                if (r > 0 && c > 0)
                    dist[idx] = std::min(dist[idx], dist[(r-1)*cols + (c-1)] + dd);
                if (r > 0 && c < cols - 1)
                    dist[idx] = std::min(dist[idx], dist[(r-1)*cols + (c+1)] + dd);
            }
            if (r % 100 == 0)
                reportProgress(0.5 * r / rows);
        }

        // Backward pass
        reportProgress(0.5, "Backward pass...");
        for (int r = rows - 1; r >= 0; --r) {
            for (int c = cols - 1; c >= 0; --c) {
                size_t idx = static_cast<size_t>(r) * cols + c;
                if (r < rows - 1)
                    dist[idx] = std::min(dist[idx], dist[(r+1)*cols + c] + dy);
                if (c < cols - 1)
                    dist[idx] = std::min(dist[idx], dist[r*cols + (c+1)] + dx);
                if (r < rows - 1 && c < cols - 1)
                    dist[idx] = std::min(dist[idx], dist[(r+1)*cols + (c+1)] + dd);
                if (r < rows - 1 && c > 0)
                    dist[idx] = std::min(dist[idx], dist[(r+1)*cols + (c-1)] + dd);
            }
            if (r % 100 == 0)
                reportProgress(0.5 + 0.5 * (rows - r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(DistanceModule)

} // namespace aplaceholder
