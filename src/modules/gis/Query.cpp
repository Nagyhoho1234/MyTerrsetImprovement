#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <QStringList>

namespace aplaceholder {

class QueryModule : public Module {
public:
    QString name() const override { return "QUERY"; }
    QString description() const override {
        return "Extract attribute values at point locations from a raster image. "
               "Reads a CSV of point coordinates and outputs a CSV with extracted values.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::file("points_file", "Points file (CSV with x,y)",
                "CSV file with x,y coordinates"),
            ParameterDef::output("output_file", "Output CSV file"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        // Read points CSV
        QFile pointsFile(parameter("points_file").toString());
        if (!pointsFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open points file: " + parameter("points_file").toString());
            return false;
        }

        struct Point { double x, y; };
        std::vector<Point> points;
        QTextStream in(&pointsFile);

        // Skip header if present
        QString header = in.readLine();
        bool hasHeader = false;
        if (!header.isEmpty()) {
            bool okX, okY;
            QStringList parts = header.split(",");
            if (parts.size() >= 2) {
                parts[0].trimmed().toDouble(&okX);
                parts[1].trimmed().toDouble(&okY);
                if (okX && okY) {
                    points.push_back({parts[0].trimmed().toDouble(),
                                      parts[1].trimmed().toDouble()});
                } else {
                    hasHeader = true;
                }
            }
        }

        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty()) continue;
            QStringList parts = line.split(",");
            if (parts.size() < 2) continue;
            bool okX, okY;
            double x = parts[0].trimmed().toDouble(&okX);
            double y = parts[1].trimmed().toDouble(&okY);
            if (okX && okY)
                points.push_back({x, y});
        }
        pointsFile.close();

        if (points.empty()) {
            setError("No valid points found in the points file");
            return false;
        }

        // Write output CSV
        QFile outFile(parameter("output_file").toString());
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file for writing");
            return false;
        }

        QTextStream out(&outFile);
        out << "x,y,value\n";

        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();
        int64_t total = static_cast<int64_t>(points.size());

        for (int64_t i = 0; i < total; ++i) {
            int col, row;
            raster->xyToColRow(points[i].x, points[i].y, col, row);

            double value;
            if (col >= 0 && col < raster->cols() && row >= 0 && row < raster->rows()) {
                value = raster->value(col, row);
                if (hasND && value == noData) {
                    out << points[i].x << "," << points[i].y << ",NoData\n";
                } else {
                    out << points[i].x << "," << points[i].y << "," << value << "\n";
                }
            } else {
                out << points[i].x << "," << points[i].y << ",OutOfBounds\n";
            }

            if (i % 10000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        outFile.close();
        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(QueryModule)

} // namespace aplaceholder
