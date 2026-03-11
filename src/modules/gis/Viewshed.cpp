#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class ViewshedModule : public Module {
public:
    QString name() const override { return "VIEWSHED"; }
    QString description() const override {
        return "Viewshed analysis. Determines which cells are visible from an observer "
               "point using line-of-sight across a DEM.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dem", "Input DEM"),
            ParameterDef::output("output", "Output viewshed image"),
            ParameterDef::real("observer_x", "Observer X coordinate", 0.0, -999999, 999999,
                               "X coordinate of observer in map units"),
            ParameterDef::real("observer_y", "Observer Y coordinate", 0.0, -999999, 999999,
                               "Y coordinate of observer in map units"),
            ParameterDef::real("observer_height", "Observer height", 1.7, 0.0, 999999,
                               "Height of observer above ground (map units)"),
            ParameterDef::real("max_distance", "Maximum distance (0=unlimited)", 0.0, 0.0, 999999,
                               "Maximum viewshed distance in map units (0 for unlimited)"),
        };
    }

    bool execute() override {
        auto dem = GdalIO::read(parameter("dem").toString());
        if (!dem) {
            setError("Failed to read DEM raster");
            return false;
        }

        double obsX = parameter("observer_x").toDouble();
        double obsY = parameter("observer_y").toDouble();
        double obsH = parameter("observer_height").toDouble();
        double maxDist = parameter("max_distance").toDouble();

        int cols = dem->cols(), rows = dem->rows();
        double noData = dem->noDataValue();
        bool hasND = dem->hasNoData();
        const auto& gt = dem->geoTransform();
        double cellW = std::abs(gt.pixelWidth);
        double cellH = std::abs(gt.pixelHeight);

        // Find observer col/row
        int obsCol, obsRow;
        dem->xyToColRow(obsX, obsY, obsCol, obsRow);

        if (obsCol < 0 || obsCol >= cols || obsRow < 0 || obsRow >= rows) {
            setError("Observer location is outside the DEM extent");
            return false;
        }

        double obsElev = dem->value(obsCol, obsRow);
        if (hasND && obsElev == noData) {
            setError("Observer location has NoData elevation");
            return false;
        }
        obsElev += obsH;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(dem->geoTransform());
        output.setProjection(dem->projection());
        output.setNoDataValue(noData);

        auto& out = output.data(0);
        int64_t total = static_cast<int64_t>(cols) * rows;
        for (int64_t i = 0; i < total; ++i)
            out[i] = 0.0;

        // Observer cell is always visible
        output.setValue(obsCol, obsRow, 1.0);

        // For each cell, trace line-of-sight from observer
        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                if (col == obsCol && row == obsRow) continue;

                double targetVal = dem->value(col, row);
                if (hasND && targetVal == noData) {
                    output.setValue(col, row, noData);
                    continue;
                }

                // Check distance
                double dx = (col - obsCol) * cellW;
                double dy = (row - obsRow) * cellH;
                double dist = std::sqrt(dx * dx + dy * dy);

                if (maxDist > 0.0 && dist > maxDist) {
                    output.setValue(col, row, 0.0);
                    continue;
                }

                // Line-of-sight check using Bresenham-like stepping
                bool visible = isVisible(*dem, obsCol, obsRow, obsElev, col, row, noData, hasND);
                output.setValue(col, row, visible ? 1.0 : 0.0);
            }

            if (row % 10 == 0)
                reportProgress(static_cast<double>(row) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }

private:
    bool isVisible(const Raster& dem, int obsCol, int obsRow, double obsElev,
                   int tgtCol, int tgtRow, double noData, bool hasND) const {
        int dc = tgtCol - obsCol;
        int dr = tgtRow - obsRow;
        int steps = std::max(std::abs(dc), std::abs(dr));
        if (steps == 0) return true;

        double stepC = static_cast<double>(dc) / steps;
        double stepR = static_cast<double>(dr) / steps;

        double tgtElev = dem.value(tgtCol, tgtRow);
        double totalDist = std::sqrt(static_cast<double>(dc * dc + dr * dr));

        // Maximum slope angle from observer to target
        double tgtAngle = (tgtElev - obsElev) / totalDist;

        // Check each intermediate cell
        double maxAngle = -1e30;
        for (int s = 1; s < steps; ++s) {
            int c = obsCol + static_cast<int>(std::round(stepC * s));
            int r = obsRow + static_cast<int>(std::round(stepR * s));

            if (c < 0 || c >= dem.cols() || r < 0 || r >= dem.rows())
                continue;

            double elev = dem.value(c, r);
            if (hasND && elev == noData) continue;

            double d = std::sqrt(
                static_cast<double>((c - obsCol) * (c - obsCol) +
                                    (r - obsRow) * (r - obsRow)));
            if (d == 0.0) continue;

            double angle = (elev - obsElev) / d;
            if (angle > maxAngle) maxAngle = angle;
        }

        return tgtAngle >= maxAngle;
    }
};

REGISTER_MODULE(ViewshedModule)

} // namespace aplaceholder
