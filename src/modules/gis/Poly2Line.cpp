#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class Poly2LineModule : public Module {
public:
    QString name() const override { return "POLY2LINE"; }
    QString description() const override {
        return "Convert polygon boundaries to line features. Traces boundaries of categorical regions and outputs a line raster.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input categorical raster"),
            ParameterDef::output("output", "Output line raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input categorical raster");
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

        // Initialize output to zero (no boundary)
        for (int64_t i = 0; i < total; ++i)
            dst[i] = 0.0;

        // Trace polygon boundaries: a pixel is on a boundary if any of its
        // 4-connected neighbors has a different category value.
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                double val = src[idx];
                bool isBoundary = false;

                // Check 4-connected neighbors (up, down, left, right)
                const int dr[] = {-1, 1, 0, 0};
                const int dc[] = {0, 0, -1, 1};
                for (int d = 0; d < 4; ++d) {
                    int nr = r + dr[d], nc = c + dc[d];
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) {
                        // Edge of raster counts as a boundary
                        isBoundary = true;
                        break;
                    }
                    double nval = src[static_cast<int64_t>(nr) * cols + nc];
                    if (nval != val) {
                        isBoundary = true;
                        break;
                    }
                }

                if (isBoundary)
                    dst[idx] = 1.0;
            }

            reportProgress(static_cast<double>(r + 1) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(Poly2LineModule)

} // namespace aplaceholder
