#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class DecompModule : public Module {
public:
    QString name() const override { return "DECOMP"; }
    QString description() const override {
        return "Decomposes magnitude/direction rasters into X and Y components, or "
               "recomposes X/Y components back to magnitude/direction. Essential for "
               "vector interpolation where direction values cannot be directly interpolated.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::combo("mode", "Operation mode",
                {"Decompose (mag/dir -> X/Y)", "Recompose (X/Y -> mag/dir)"}, 0,
                "Decompose splits magnitude/direction into X/Y; Recompose does the reverse"),
            ParameterDef::file("input1", "Input 1 (Magnitude or X component)",
                "Decompose: magnitude raster. Recompose: X component raster"),
            ParameterDef::file("input2", "Input 2 (Direction or Y component)",
                "Decompose: direction raster (degrees). Recompose: Y component raster"),
            ParameterDef::output("output1", "Output 1 (X component or Magnitude)"),
            ParameterDef::output("output2", "Output 2 (Y component or Direction)"),
        };
    }

    bool execute() override {
        int mode = parameter("mode").toInt();
        QString in1Path = parameter("input1").toString();
        QString in2Path = parameter("input2").toString();
        QString out1Path = parameter("output1").toString();
        QString out2Path = parameter("output2").toString();

        auto raster1 = GdalIO::read(in1Path);
        auto raster2 = GdalIO::read(in2Path);
        if (!raster1) { setError("Failed to read input 1: " + in1Path); return false; }
        if (!raster2) { setError("Failed to read input 2: " + in2Path); return false; }

        int cols = raster1->cols();
        int rows = raster1->rows();
        if (cols != raster2->cols() || rows != raster2->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        const auto& data1 = raster1->data(0);
        const auto& data2 = raster2->data(0);
        double noData1 = raster1->hasNoData() ? raster1->noDataValue() : -9999.0;
        double noData2 = raster2->hasNoData() ? raster2->noDataValue() : -9999.0;

        Raster out1(cols, rows, 1, DataType::Float64);
        Raster out2(cols, rows, 1, DataType::Float64);
        out1.setGeoTransform(raster1->geoTransform());
        out1.setProjection(raster1->projection());
        out1.setNoDataValue(-9999.0);
        out2.setGeoTransform(raster1->geoTransform());
        out2.setProjection(raster1->projection());
        out2.setNoDataValue(-9999.0);

        auto& outData1 = out1.data(0);
        auto& outData2 = out2.data(0);
        double outNoData = -9999.0;

        int64_t totalCells = static_cast<int64_t>(cols) * rows;
        static constexpr double DEG2RAD = 3.14159265358979323846 / 180.0;
        static constexpr double RAD2DEG = 180.0 / 3.14159265358979323846;

        if (mode == 0) {
            // Decompose: magnitude + direction -> X, Y
            reportProgress(0.0, "Decomposing to X/Y components...");
            for (int64_t i = 0; i < totalCells; ++i) {
                double mag = data1[i];
                double dir = data2[i];

                if ((raster1->hasNoData() && mag == noData1) ||
                    (raster2->hasNoData() && dir == noData2)) {
                    outData1[i] = outNoData;
                    outData2[i] = outNoData;
                    continue;
                }

                double dirRad = dir * DEG2RAD;
                outData1[i] = mag * std::cos(dirRad);
                outData2[i] = mag * std::sin(dirRad);
            }
        } else {
            // Recompose: X + Y -> magnitude, direction
            reportProgress(0.0, "Recomposing to magnitude/direction...");
            for (int64_t i = 0; i < totalCells; ++i) {
                double x = data1[i];
                double y = data2[i];

                if ((raster1->hasNoData() && x == noData1) ||
                    (raster2->hasNoData() && y == noData2)) {
                    outData1[i] = outNoData;
                    outData2[i] = outNoData;
                    continue;
                }

                outData1[i] = std::sqrt(x * x + y * y);
                double dir = std::atan2(y, x) * RAD2DEG;
                if (dir < 0.0) dir += 360.0;
                outData2[i] = dir;
            }
        }

        reportProgress(0.8, "Writing outputs...");
        if (!GdalIO::write(out1, out1Path)) {
            setError("Failed to write output 1");
            return false;
        }
        if (!GdalIO::write(out2, out2Path)) {
            setError("Failed to write output 2");
            return false;
        }

        reportProgress(1.0, "Done");
        return true;
    }
};

REGISTER_MODULE(DecompModule)

} // namespace aplaceholder
