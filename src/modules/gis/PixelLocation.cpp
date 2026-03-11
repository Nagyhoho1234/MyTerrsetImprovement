#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class PixelLocationModule : public Module {
public:
    QString name() const override { return "PIXELLOCATION"; }
    QString description() const override {
        return "Create coordinate images from a reference raster. "
               "Generates X (easting) or Y (northing) coordinate values "
               "for each pixel based on the raster's geo-referencing.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input reference raster"),
            ParameterDef::output("output", "Output coordinate image"),
            ParameterDef::combo("coordinate", "Coordinate type",
                {"X (Easting/Longitude)", "Y (Northing/Latitude)"}, 0,
                "Which coordinate axis to output"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        int coordType = parameter("coordinate").toInt();

        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        const auto& data = input->data(0);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        reportProgress(0.0, "Generating coordinate image...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;
                if (hasND && data[idx] == noData) {
                    out[idx] = noData;
                    continue;
                }

                double x, y;
                input->colRowToXY(c, r, x, y);
                out[idx] = (coordType == 0) ? x : y;
            }
            if (r % 200 == 0) reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PixelLocationModule)

} // namespace aplaceholder
