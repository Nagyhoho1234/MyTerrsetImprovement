#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace aplaceholder {

class PitRemovalModule : public Module {
public:
    QString name() const override { return "PITREMOVAL"; }
    QString description() const override {
        return "Fill sinks/pits in a DEM. Iteratively raises pit cells to the level "
               "of their lowest neighbor until no pits remain.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_dem", "Input DEM"),
            ParameterDef::output("output", "Output filled DEM"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input_dem").toString());
        if (!r1) {
            setError("Failed to read input DEM");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();

        // Create output as a copy of input
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        if (hasND) output.setNoDataValue(noData);

        const auto& src = r1->data(0);
        auto& out = output.data(0);
        int64_t total = static_cast<int64_t>(cols) * rows;
        for (int64_t i = 0; i < total; ++i)
            out[i] = src[i];

        // 8-connected neighbor offsets
        const int dx[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
        const int dy[8] = {-1, -1, -1, 0, 0, 1, 1, 1};

        // Iteratively fill pits
        bool changed = true;
        int iteration = 0;
        const int maxIterations = 10000;

        while (changed && iteration < maxIterations) {
            changed = false;
            ++iteration;

            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    double val = output.value(col, row);
                    if (hasND && val == noData) continue;

                    // Check if this cell is a pit (lower than all neighbors)
                    double minNeighbor = std::numeric_limits<double>::max();
                    bool isPit = true;
                    bool hasValidNeighbor = false;

                    for (int n = 0; n < 8; ++n) {
                        int nc = col + dx[n];
                        int nr = row + dy[n];
                        if (nc < 0 || nc >= cols || nr < 0 || nr >= rows) continue;

                        double nv = output.value(nc, nr);
                        if (hasND && nv == noData) continue;

                        hasValidNeighbor = true;
                        if (nv <= val) {
                            isPit = false;
                            break;
                        }
                        if (nv < minNeighbor) {
                            minNeighbor = nv;
                        }
                    }

                    // If it is a pit (lower than all valid neighbors), raise it
                    if (isPit && hasValidNeighbor && minNeighbor < std::numeric_limits<double>::max()) {
                        output.setValue(col, row, minNeighbor);
                        changed = true;
                    }
                }
            }

            if (iteration % 10 == 0)
                reportProgress(0.0, QString("Iteration %1...").arg(iteration));
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PitRemovalModule)

} // namespace aplaceholder
