#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

namespace aplaceholder {

struct IDWPoint {
    double x, y, value;
};

class InterpolModule : public Module {
public:
    QString name() const override { return "INTERPOL"; }
    QString description() const override {
        return "Interpolates a continuous raster surface from point data using "
               "Inverse Distance Weighting (IDW). Exact interpolator that preserves "
               "original sample values at point locations.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("points_file", "Input point data (CSV: x,y,value)",
                "CSV file with columns x, y, value"),
            ParameterDef::file("reference_raster", "Reference raster (extent/dimensions)",
                "Raster defining the output extent and resolution"),
            ParameterDef::output("output", "Output interpolated surface"),
            ParameterDef::real("power", "Distance weighting exponent", 2.0, 0.1, 10.0,
                "Higher values give more weight to nearest points"),
            ParameterDef::real("search_radius", "Search radius (0 = all points)", 0.0, 0.0, 999999.0,
                "Maximum distance to search for points; 0 uses all points"),
            ParameterDef::integer("num_neighbors", "Maximum number of neighbors", 12, 1, 99999,
                "Maximum number of nearest points to use"),
        };
    }

    bool execute() override {
        QString pointsPath = parameter("points_file").toString();
        QString refPath = parameter("reference_raster").toString();
        QString outPath = parameter("output").toString();
        double power = parameter("power").toDouble();
        double searchRadius = parameter("search_radius").toDouble();
        int numNeighbors = parameter("num_neighbors").toInt();

        // Read points from CSV
        std::vector<IDWPoint> points;
        {
            std::ifstream ifs(pointsPath.toStdString());
            if (!ifs.is_open()) {
                setError("Failed to open points file: " + pointsPath);
                return false;
            }
            std::string line;
            // Skip header
            std::getline(ifs, line);
            while (std::getline(ifs, line)) {
                if (line.empty()) continue;
                std::istringstream ss(line);
                IDWPoint pt;
                char sep;
                if (ss >> pt.x >> sep >> pt.y >> sep >> pt.value) {
                    points.push_back(pt);
                }
            }
        }

        if (points.empty()) {
            setError("No valid points read from CSV file");
            return false;
        }

        // Read reference raster for extent and dimensions
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
        double noData = output.noDataValue();
        bool useRadius = (searchRadius > 0.0);
        double radiusSq = searchRadius * searchRadius;

        // Structure for sorting by distance
        struct DistVal {
            double distSq;
            double value;
            bool operator<(const DistVal& o) const { return distSq < o.distSq; }
        };

        reportProgress(0.0, "Interpolating surface...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                double px, py;
                output.colRowToXY(c, r, px, py);

                // Gather candidate points and their distances
                std::vector<DistVal> candidates;
                candidates.reserve(points.size());

                for (const auto& pt : points) {
                    double dx = px - pt.x;
                    double dy = py - pt.y;
                    double dSq = dx * dx + dy * dy;

                    if (useRadius && dSq > radiusSq)
                        continue;

                    candidates.push_back({dSq, pt.value});
                }

                if (candidates.empty()) {
                    outData[static_cast<size_t>(r) * cols + c] = noData;
                    continue;
                }

                // Sort by distance and take nearest neighbors
                std::sort(candidates.begin(), candidates.end());
                int count = std::min(static_cast<int>(candidates.size()), numNeighbors);

                // Check if pixel coincides with a point (distance ~ 0)
                if (candidates[0].distSq < 1e-20) {
                    outData[static_cast<size_t>(r) * cols + c] = candidates[0].value;
                    continue;
                }

                // Compute IDW
                double weightSum = 0.0;
                double valueSum = 0.0;
                for (int i = 0; i < count; ++i) {
                    double dist = std::sqrt(candidates[i].distSq);
                    double w = 1.0 / std::pow(dist, power);
                    weightSum += w;
                    valueSum += w * candidates[i].value;
                }

                outData[static_cast<size_t>(r) * cols + c] = valueSum / weightSum;
            }

            if (r % 50 == 0)
                reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outPath);
    }
};

REGISTER_MODULE(InterpolModule)

} // namespace aplaceholder
