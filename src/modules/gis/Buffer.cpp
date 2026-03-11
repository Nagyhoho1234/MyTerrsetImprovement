#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <queue>
#include <limits>

namespace aplaceholder {

class BufferModule : public Module {
public:
    QString name() const override { return "BUFFER"; }
    QString description() const override {
        return "Buffer around non-zero features. Computes Euclidean distance from features "
               "and outputs 1 where distance <= buffer_distance, else 0.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input feature image"),
            ParameterDef::output("output", "Output image"),
            ParameterDef::real("buffer_distance", "Buffer distance (map units)", 100.0, 0.0, 999999,
                               "Maximum distance from features to buffer"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input").toString());
        if (!r1) {
            setError("Failed to read input raster");
            return false;
        }

        double bufDist = parameter("buffer_distance").toDouble();
        int cols = r1->cols(), rows = r1->rows();
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();
        const auto& gt = r1->geoTransform();
        double cellW = std::abs(gt.pixelWidth);
        double cellH = std::abs(gt.pixelHeight);

        // Convert buffer distance to pixel radius for search window
        int searchR = static_cast<int>(std::ceil(bufDist / std::min(cellW, cellH)));

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(noData);

        // Compute distance via brute-force within search window
        // For large rasters a proper distance transform would be faster,
        // but this is correct and straightforward.
        const auto& d1 = r1->data(0);
        auto& out = output.data(0);
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Build list of feature cell locations for efficiency
        std::vector<std::pair<int, int>> featureCells;
        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                double val = r1->value(col, row);
                if (hasND && val == noData) continue;
                if (val != 0.0) {
                    featureCells.emplace_back(col, row);
                }
            }
        }

        // For each cell, check distance to nearest feature within buffer
        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                double val = r1->value(col, row);
                if (hasND && val == noData) {
                    output.setValue(col, row, noData);
                    continue;
                }

                // If cell is already a feature, it's within the buffer
                if (val != 0.0) {
                    output.setValue(col, row, 1.0);
                    continue;
                }

                // Search within window for nearest feature
                double minDist = std::numeric_limits<double>::max();
                int rMin = std::max(0, row - searchR);
                int rMax = std::min(rows - 1, row + searchR);
                int cMin = std::max(0, col - searchR);
                int cMax = std::min(cols - 1, col + searchR);

                for (int sr = rMin; sr <= rMax; ++sr) {
                    for (int sc = cMin; sc <= cMax; ++sc) {
                        double sv = r1->value(sc, sr);
                        if (hasND && sv == noData) continue;
                        if (sv != 0.0) {
                            double dx = (sc - col) * cellW;
                            double dy = (sr - row) * cellH;
                            double dist = std::sqrt(dx * dx + dy * dy);
                            if (dist < minDist) minDist = dist;
                        }
                    }
                }

                output.setValue(col, row, (minDist <= bufDist) ? 1.0 : 0.0);
            }

            if (row % 10 == 0)
                reportProgress(static_cast<double>(row) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(BufferModule)

} // namespace aplaceholder
