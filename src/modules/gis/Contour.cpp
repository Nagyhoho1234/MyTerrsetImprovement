#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ContourModule : public Module {
public:
    QString name() const override { return "CONTOUR"; }
    QString description() const override {
        return "Extract contour lines as a raster. Marks pixels where values cross "
               "contour interval boundaries by checking 4-connected neighbors.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input DEM"),
            ParameterDef::output("output", "Output image"),
            ParameterDef::real("contour_interval", "Contour interval", 100.0, 0.001, 999999,
                               "Elevation interval between contour lines"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input").toString());
        if (!r1) {
            setError("Failed to read input raster");
            return false;
        }

        double interval = parameter("contour_interval").toDouble();
        if (interval <= 0.0) {
            setError("Contour interval must be positive");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(noData);

        auto& out = output.data(0);
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Initialize output to 0
        for (int64_t i = 0; i < total; ++i)
            out[i] = 0.0;

        // For each pixel, check if a contour line crosses between it and its neighbors
        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                double val = r1->value(col, row);
                if (hasND && val == noData) {
                    output.setValue(col, row, noData);
                    continue;
                }

                int contourLevel = static_cast<int>(std::floor(val / interval));

                // Check right neighbor
                if (col + 1 < cols) {
                    double nv = r1->value(col + 1, row);
                    if (!(hasND && nv == noData)) {
                        int nLevel = static_cast<int>(std::floor(nv / interval));
                        if (contourLevel != nLevel) {
                            output.setValue(col, row, 1.0);
                        }
                    }
                }

                // Check bottom neighbor
                if (row + 1 < rows) {
                    double nv = r1->value(col, row + 1);
                    if (!(hasND && nv == noData)) {
                        int nLevel = static_cast<int>(std::floor(nv / interval));
                        if (contourLevel != nLevel) {
                            output.setValue(col, row, 1.0);
                        }
                    }
                }
            }

            if (row % 100 == 0)
                reportProgress(static_cast<double>(row) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ContourModule)

} // namespace aplaceholder
