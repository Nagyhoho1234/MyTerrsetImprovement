#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>

namespace aplaceholder {

class RasterToPointModule : public Module {
public:
    QString name() const override { return "RASTERTOPOINT"; }
    QString description() const override {
        return "Convert raster to point CSV. For each non-NoData pixel, "
               "outputs x, y, value.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output_file", "Output CSV file"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        QFile outFile(parameter("output_file").toString());
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file for writing");
            return false;
        }

        QTextStream stream(&outFile);
        stream << "x,y,value\n";

        int cols = raster->cols();
        int rows = raster->rows();
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();
        const auto& d = raster->data(0);
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                int64_t idx = static_cast<int64_t>(row) * cols + col;
                double val = d[idx];

                if (hasND && val == noData)
                    continue;

                double x, y;
                raster->colRowToXY(col, row, x, y);
                stream << x << "," << y << "," << val << "\n";
            }

            if (row % 100 == 0)
                reportProgress(static_cast<double>(row) / rows);
        }

        reportProgress(1.0, "Done");
        return true;
    }
};

REGISTER_MODULE(RasterToPointModule)

} // namespace aplaceholder
