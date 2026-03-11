#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class PathwayModule : public Module {
public:
    QString name() const override { return "PATHWAY"; }
    QString description() const override {
        return "Finds the least-cost path between a destination point and the nearest "
               "source by tracing back through a cost distance and backlink surface.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("cost_distance_raster", "Cost distance raster"),
            ParameterDef::file("backlink_raster", "Backlink (direction) raster"),
            ParameterDef::real("destination_x", "Destination X coordinate", 0, -1e12, 1e12),
            ParameterDef::real("destination_y", "Destination Y coordinate", 0, -1e12, 1e12),
            ParameterDef::output("output", "Output path raster"),
        };
    }

    bool execute() override {
        auto costDist = GdalIO::read(parameter("cost_distance_raster").toString());
        auto backlink = GdalIO::read(parameter("backlink_raster").toString());
        if (!costDist || !backlink) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = costDist->cols(), rows = costDist->rows();
        if (backlink->cols() != cols || backlink->rows() != rows) {
            setError("Cost distance and backlink rasters must have the same dimensions");
            return false;
        }

        double destX = parameter("destination_x").toDouble();
        double destY = parameter("destination_y").toDouble();

        int destCol, destRow;
        costDist->xyToColRow(destX, destY, destCol, destRow);

        if (destCol < 0 || destCol >= cols || destRow < 0 || destRow >= rows) {
            setError("Destination point is outside the raster extent");
            return false;
        }

        const auto& costData = costDist->data(0);
        const auto& blData = backlink->data(0);

        // Backlink encoding (reverse direction to trace back):
        // 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE
        // To trace back, we go in the opposite direction of the backlink
        // Backlink tells us where the flow came FROM, so we follow it directly
        const int blDr[] = {0, 1, 1, 1, 0, -1, -1, -1};
        const int blDc[] = {1, 1, 0, -1, -1, -1, 0, 1};
        const int blCode[] = {1, 2, 4, 8, 16, 32, 64, 128};

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(costDist->geoTransform());
        output.setProjection(costDist->projection());
        output.setNoDataValue(0);
        auto& out = output.data(0);

        // Trace back from destination to source
        reportProgress(0.0, "Tracing least-cost path...");
        int curRow = destRow, curCol = destCol;
        int maxSteps = cols * rows; // safety limit

        for (int step = 0; step < maxSteps; ++step) {
            int64_t idx = static_cast<int64_t>(curRow) * cols + curCol;
            out[idx] = 1.0;

            // If cost is 0 or near 0, we reached a source
            if (costData[idx] <= 0.0)
                break;

            int bl = static_cast<int>(blData[idx]);
            if (bl <= 0) break;

            // Find direction from backlink code
            int nextRow = curRow, nextCol = curCol;
            bool found = false;
            for (int d = 0; d < 8; ++d) {
                if (blCode[d] == bl) {
                    // Backlink points toward the source, so we move in that direction
                    nextRow = curRow + blDr[d];
                    nextCol = curCol + blDc[d];
                    found = true;
                    break;
                }
            }

            if (!found) break;
            if (nextRow < 0 || nextRow >= rows || nextCol < 0 || nextCol >= cols)
                break;

            curRow = nextRow;
            curCol = nextCol;
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PathwayModule)

} // namespace aplaceholder
