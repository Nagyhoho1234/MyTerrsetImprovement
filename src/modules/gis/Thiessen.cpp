#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

namespace aplaceholder {

class ThiessenModule : public Module {
public:
    QString name() const override { return "THIESSEN"; }
    QString description() const override {
        return "Creates a Thiessen (Voronoi) tessellation from point data. Each pixel "
               "is assigned the ID of the nearest data point based on Euclidean distance.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("points_file", "Input point data (CSV: x,y,id)",
                "CSV file with columns x, y, id"),
            ParameterDef::file("reference_raster", "Reference raster (extent/dimensions)",
                "Raster defining the output extent and resolution"),
            ParameterDef::output("output", "Output tessellation raster"),
        };
    }

    bool execute() override {
        QString pointsPath = parameter("points_file").toString();
        QString refPath = parameter("reference_raster").toString();
        QString outPath = parameter("output").toString();

        // Read points
        struct Point { double x, y; int id; };
        std::vector<Point> points;
        {
            std::ifstream ifs(pointsPath.toStdString());
            if (!ifs.is_open()) {
                setError("Failed to open points file: " + pointsPath);
                return false;
            }
            std::string line;
            std::getline(ifs, line); // skip header
            while (std::getline(ifs, line)) {
                if (line.empty()) continue;
                std::istringstream ss(line);
                Point pt;
                char sep;
                if (ss >> pt.x >> sep >> pt.y >> sep >> pt.id) {
                    points.push_back(pt);
                }
            }
        }

        if (points.empty()) {
            setError("No valid points read from CSV file");
            return false;
        }

        // Read reference raster metadata
        auto ref = GdalIO::readMetadata(refPath);
        if (!ref) {
            setError("Failed to read reference raster: " + refPath);
            return false;
        }

        int cols = ref->cols();
        int rows = ref->rows();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(ref->geoTransform());
        output.setProjection(ref->projection());
        output.setNoDataValue(-9999.0);

        auto& outData = output.data(0);

        reportProgress(0.0, "Computing Thiessen polygons...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                double px, py;
                output.colRowToXY(c, r, px, py);

                double minDistSq = std::numeric_limits<double>::max();
                int nearestId = -1;

                for (const auto& pt : points) {
                    double dx = px - pt.x;
                    double dy = py - pt.y;
                    double dSq = dx * dx + dy * dy;
                    if (dSq < minDistSq) {
                        minDistSq = dSq;
                        nearestId = pt.id;
                    }
                }

                outData[static_cast<size_t>(r) * cols + c] = static_cast<double>(nearestId);
            }

            if (r % 50 == 0)
                reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outPath);
    }
};

REGISTER_MODULE(ThiessenModule)

} // namespace aplaceholder
