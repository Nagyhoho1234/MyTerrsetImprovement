#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <map>
#include <fstream>

namespace aplaceholder {

class PerimModule : public Module {
public:
    QString name() const override { return "PERIM"; }
    QString description() const override {
        return "Perimeter calculation for categorical rasters. Counts edge pixels "
               "(adjacent to different class) for each class and outputs a CSV with "
               "class, area, and perimeter.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input categorical raster"),
            ParameterDef::output("output_file", "Output CSV file"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input raster"); return false; }

        int cols = raster->cols(), rows = raster->rows();
        const auto& data = raster->data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        double dx = std::abs(raster->geoTransform().pixelWidth);
        double dy = std::abs(raster->geoTransform().pixelHeight);
        double cellArea = dx * dy;

        // 4-connected neighbors for perimeter calculation
        const int dr[] = {-1, 0, 0, 1};
        const int dc[] = {0, -1, 1, 0};

        // Count area and perimeter edges per class
        std::map<int, int64_t> areaCount;   // class -> pixel count
        std::map<int, int64_t> perimCount;   // class -> edge pixel count

        reportProgress(0.0, "Computing perimeter...");
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                double val = data[idx];

                if (hasND && val == noData) continue;

                int classVal = static_cast<int>(val);
                areaCount[classVal]++;

                // Check if this pixel is on an edge (adjacent to different class or boundary)
                bool isEdge = false;
                for (int d = 0; d < 4; ++d) {
                    int nr = r + dr[d];
                    int nc = c + dc[d];

                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) {
                        isEdge = true;
                        break;
                    }

                    int64_t nIdx = static_cast<int64_t>(nr) * cols + nc;
                    double nVal = data[nIdx];

                    if (hasND && nVal == noData) {
                        isEdge = true;
                        break;
                    }

                    if (static_cast<int>(nVal) != classVal) {
                        isEdge = true;
                        break;
                    }
                }

                if (isEdge) {
                    perimCount[classVal]++;
                }
            }
            if (r % 100 == 0)
                reportProgress(0.9 * static_cast<double>(r) / rows);
        }

        // Write CSV output
        reportProgress(0.9, "Writing CSV...");
        QString outPath = parameter("output_file").toString();
        std::ofstream ofs(outPath.toStdString());
        if (!ofs.is_open()) {
            setError("Failed to open output file: " + outPath);
            return false;
        }

        ofs << "Class,Area,Perimeter\n";
        for (auto& [cls, area] : areaCount) {
            double areaVal = area * cellArea;
            // Perimeter: count of edge pixels * cell edge length (approximate)
            double perimVal = perimCount[cls] * ((dx + dy) / 2.0);
            ofs << cls << "," << areaVal << "," << perimVal << "\n";
        }

        ofs.close();
        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(PerimModule)

} // namespace aplaceholder
