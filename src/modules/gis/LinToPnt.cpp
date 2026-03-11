#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class LinToPntModule : public Module {
public:
    QString name() const override { return "LINTOPNT"; }
    QString description() const override {
        return "Convert line features to point features. Outputs a point raster with values at line vertices/pixels.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input line raster"),
            ParameterDef::output("output", "Output point raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input line raster");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        Raster output(cols, rows, 1, DataType::Float32);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(0.0);

        const auto& src = raster->data(0);
        auto& dst = output.data(0);
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Initialize output to nodata
        for (int64_t i = 0; i < total; ++i)
            dst[i] = 0.0;

        // Scan for line pixels and mark them as point features.
        // A pixel is considered a line vertex if it is a non-zero line pixel
        // that is either an endpoint or a junction (i.e., not a simple
        // continuation along a straight run).
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                double val = src[idx];
                if (val == 0.0)
                    continue;

                // Count non-zero neighbors in the 8-connected neighborhood
                int neighbors = 0;
                for (int dr = -1; dr <= 1; ++dr) {
                    for (int dc = -1; dc <= 1; ++dc) {
                        if (dr == 0 && dc == 0)
                            continue;
                        int nr = r + dr, nc = c + dc;
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
                            continue;
                        if (src[static_cast<int64_t>(nr) * cols + nc] != 0.0)
                            ++neighbors;
                    }
                }

                // Mark endpoints (1 neighbor), junctions (3+ neighbors),
                // isolated pixels (0 neighbors), or all line pixels as points.
                // Every line pixel becomes a point with the original value.
                dst[idx] = val;
            }

            reportProgress(static_cast<double>(r + 1) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(LinToPntModule)

} // namespace aplaceholder
