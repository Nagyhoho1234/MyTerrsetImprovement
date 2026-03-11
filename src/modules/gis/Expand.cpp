#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ExpandModule : public Module {
public:
    QString name() const override { return "EXPAND"; }
    QString description() const override {
        return "Expand (dilate) non-zero features by N pixels using a circular structuring element.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output image"),
            ParameterDef::integer("radius", "Radius (pixels)", 1, 1, 100,
                                  "Number of pixels to expand features"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input").toString());
        if (!r1) {
            setError("Failed to read input raster");
            return false;
        }

        int radius = parameter("radius").toInt();
        int cols = r1->cols(), rows = r1->rows();
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(noData);

        auto& out = output.data(0);
        double radiusSq = static_cast<double>(radius) * radius;

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                double val = r1->value(col, row);
                if (hasND && val == noData) {
                    output.setValue(col, row, noData);
                    continue;
                }

                // If already non-zero, keep it
                if (val != 0.0) {
                    output.setValue(col, row, val);
                    continue;
                }

                // Check if any non-zero neighbor within radius
                bool found = false;
                int rMin = std::max(0, row - radius);
                int rMax = std::min(rows - 1, row + radius);
                int cMin = std::max(0, col - radius);
                int cMax = std::min(cols - 1, col + radius);

                for (int sr = rMin; sr <= rMax && !found; ++sr) {
                    for (int sc = cMin; sc <= cMax && !found; ++sc) {
                        double dx = sc - col;
                        double dy = sr - row;
                        if (dx * dx + dy * dy <= radiusSq) {
                            double sv = r1->value(sc, sr);
                            if (!(hasND && sv == noData) && sv != 0.0) {
                                found = true;
                            }
                        }
                    }
                }

                output.setValue(col, row, found ? 1.0 : 0.0);
            }

            if (row % 100 == 0)
                reportProgress(static_cast<double>(row) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ExpandModule)

} // namespace aplaceholder
