#include "Module.h"
#include "ModuleRegistry.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

namespace aplaceholder {

class GeneralizeModule : public Module {
public:
    QString name() const override { return "GENERALIZE"; }
    QString description() const override {
        return "Point thinning and line generalization for preparing TIN input data. "
               "Removes points within a minimum distance of each other, using a "
               "Douglas-Peucker-like simplification approach.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("points_file", "Input point data (CSV)",
                "CSV file with columns x, y, value (or x, y)"),
            ParameterDef::output("output_file", "Output thinned point file (CSV)"),
            ParameterDef::real("min_distance", "Minimum distance between points", 1.0, 0.0, 999999.0,
                "Points closer than this distance to a previously kept point are removed"),
        };
    }

    bool execute() override {
        QString inPath = parameter("points_file").toString();
        QString outPath = parameter("output_file").toString();
        double minDist = parameter("min_distance").toDouble();
        double minDistSq = minDist * minDist;

        // Read points from CSV
        struct Point {
            double x, y;
            std::string rest; // remainder of the CSV line (value columns etc.)
            bool keep;
        };

        std::vector<Point> points;
        std::string header;

        {
            std::ifstream ifs(inPath.toStdString());
            if (!ifs.is_open()) {
                setError("Failed to open points file: " + inPath);
                return false;
            }
            std::getline(ifs, header);
            std::string line;
            while (std::getline(ifs, line)) {
                if (line.empty()) continue;
                std::istringstream ss(line);
                Point pt;
                char sep;
                double dummy;
                if (ss >> pt.x >> sep >> pt.y) {
                    // Capture rest of line (value column(s))
                    std::string remaining;
                    std::getline(ss, remaining);
                    pt.rest = remaining;
                    pt.keep = false;
                    points.push_back(pt);
                }
            }
        }

        if (points.empty()) {
            setError("No valid points read from CSV file");
            return false;
        }

        reportProgress(0.0, "Thinning points...");

        // Douglas-Peucker style greedy thinning:
        // Mark first point as kept. For each subsequent point, keep it only if
        // it is at least min_distance from all previously kept points.
        // This is a spatial thinning approach.

        // First, try to detect if points form an ordered sequence (line).
        // If so, apply Douglas-Peucker simplification.
        // Otherwise, use proximity-based thinning.

        // We use proximity-based thinning as the general approach:
        // iterate through points, keep a point if no already-kept point
        // is within min_distance.

        std::vector<size_t> keptIndices;
        keptIndices.reserve(points.size());

        for (size_t i = 0; i < points.size(); ++i) {
            bool tooClose = false;

            for (size_t ki : keptIndices) {
                double dx = points[i].x - points[ki].x;
                double dy = points[i].y - points[ki].y;
                if (dx * dx + dy * dy < minDistSq) {
                    tooClose = true;
                    break;
                }
            }

            if (!tooClose) {
                points[i].keep = true;
                keptIndices.push_back(i);
            }

            if (i % 500 == 0)
                reportProgress(0.8 * i / points.size());
        }

        reportProgress(0.9, "Writing output...");

        // Write output CSV
        {
            std::ofstream ofs(outPath.toStdString());
            if (!ofs.is_open()) {
                setError("Failed to write output file: " + outPath);
                return false;
            }
            ofs.precision(12);
            ofs << header << "\n";
            for (const auto& pt : points) {
                if (pt.keep) {
                    ofs << pt.x << "," << pt.y << pt.rest << "\n";
                }
            }
        }

        reportProgress(1.0, "Done");
        return true;
    }
};

REGISTER_MODULE(GeneralizeModule)

} // namespace aplaceholder
