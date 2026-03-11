#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class FlipModule : public Module {
public:
    QString name() const override { return "FLIP"; }
    QString description() const override {
        return "Flip raster horizontally, vertically, or both.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output image"),
            ParameterDef::combo("direction", "Flip direction",
                {"horizontal", "vertical", "both"}, 0,
                "Direction of flip"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input").toString());
        if (!r1) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        int direction = parameter("direction").toInt();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        if (r1->hasNoData())
            output.setNoDataValue(r1->noDataValue());

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                int srcCol = col, srcRow = row;

                switch (direction) {
                    case 0: // horizontal
                        srcCol = cols - 1 - col;
                        break;
                    case 1: // vertical
                        srcRow = rows - 1 - row;
                        break;
                    case 2: // both
                        srcCol = cols - 1 - col;
                        srcRow = rows - 1 - row;
                        break;
                }

                output.setValue(col, row, r1->value(srcCol, srcRow));
            }

            if (row % 100 == 0)
                reportProgress(static_cast<double>(row) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(FlipModule)

} // namespace aplaceholder
