#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <cmath>
#include <limits>

namespace aplaceholder {

class PointToRasterModule : public Module {
public:
    QString name() const override { return "POINTTORASTER"; }
    QString description() const override {
        return "Rasterize point CSV (x,y,value) onto a reference grid. "
               "Handles multiple points per cell using mean, sum, max, min, or count.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("points_file", "Input points CSV (x,y,value)"),
            ParameterDef::file("reference_raster", "Reference raster (defines grid)"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::combo("method", "Aggregation method",
                {"Mean", "Sum", "Max", "Min", "Count"}, 0,
                "Method for combining multiple points per cell"),
        };
    }

    bool execute() override {
        auto ref = GdalIO::read(parameter("reference_raster").toString());
        if (!ref) {
            setError("Failed to read reference raster");
            return false;
        }

        int cols = ref->cols(), rows = ref->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int method = parameter("method").toInt();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(ref->geoTransform());
        output.setProjection(ref->projection());
        double noData = -9999.0;
        output.setNoDataValue(noData);

        auto& out = output.data(0);
        std::vector<double> sums(total, 0.0);
        std::vector<int> counts(total, 0);
        std::vector<double> extremes(total, 0.0);

        // Initialize extremes
        if (method == 2) // Max
            std::fill(extremes.begin(), extremes.end(), -std::numeric_limits<double>::max());
        else if (method == 3) // Min
            std::fill(extremes.begin(), extremes.end(), std::numeric_limits<double>::max());

        // Read CSV
        QFile csvFile(parameter("points_file").toString());
        if (!csvFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open points CSV file");
            return false;
        }

        QTextStream stream(&csvFile);
        QString header = stream.readLine(); // skip header

        reportProgress(0.0, "Reading points...");
        int lineCount = 0;

        while (!stream.atEnd()) {
            QString line = stream.readLine().trimmed();
            if (line.isEmpty()) continue;

            QStringList parts = line.split(',');
            if (parts.size() < 3) continue;

            double x = parts[0].toDouble();
            double y = parts[1].toDouble();
            double val = parts[2].toDouble();

            int col, row;
            ref->xyToColRow(x, y, col, row);

            if (col < 0 || col >= cols || row < 0 || row >= rows)
                continue;

            int64_t idx = static_cast<int64_t>(row) * cols + col;
            counts[idx]++;

            switch (method) {
                case 0: // Mean
                case 1: // Sum
                    sums[idx] += val;
                    break;
                case 2: // Max
                    if (val > extremes[idx]) extremes[idx] = val;
                    break;
                case 3: // Min
                    if (val < extremes[idx]) extremes[idx] = val;
                    break;
                case 4: // Count
                    break;
            }

            if (++lineCount % 100000 == 0)
                reportProgress(0.5, QString("Read %1 points...").arg(lineCount));
        }

        reportProgress(0.7, "Building output raster...");

        for (int64_t i = 0; i < total; ++i) {
            if (counts[i] == 0) {
                out[i] = noData;
            } else {
                switch (method) {
                    case 0: out[i] = sums[i] / counts[i]; break; // Mean
                    case 1: out[i] = sums[i]; break;             // Sum
                    case 2: out[i] = extremes[i]; break;         // Max
                    case 3: out[i] = extremes[i]; break;         // Min
                    case 4: out[i] = static_cast<double>(counts[i]); break; // Count
                }
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PointToRasterModule)

} // namespace aplaceholder
