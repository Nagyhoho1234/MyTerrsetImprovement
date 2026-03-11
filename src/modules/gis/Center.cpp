#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <map>
#include <fstream>

namespace aplaceholder {

class CenterModule : public Module {
public:
    QString name() const override { return "CENTER"; }
    QString description() const override {
        return "Compute the spatial centroid (center of gravity) for each category "
               "in a categorical raster. Outputs a text file with X/Y coordinates "
               "for each category's mean center and optionally a point raster.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input categorical raster image"),
            ParameterDef::output("output", "Output centroid text file"),
            ParameterDef::output("output_raster", "Output centroid point raster (optional)"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();

        reportProgress(0.0, "Computing category centroids...");

        struct CentroidAccum {
            double sumX = 0, sumY = 0;
            int64_t count = 0;
        };
        std::map<int, CentroidAccum> accum;

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;
                if (hasND && data[idx] == noData) continue;

                int cat = static_cast<int>(std::round(data[idx]));
                double x, y;
                input->colRowToXY(c, r, x, y);
                accum[cat].sumX += x;
                accum[cat].sumY += y;
                accum[cat].count++;
            }
            if (r % 200 == 0) reportProgress(0.5 * static_cast<double>(r) / rows);
        }

        reportProgress(0.6, "Writing centroid results...");

        // Write text output
        QString outputPath = parameter("output").toString();
        std::ofstream outFile(outputPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output file: " + outputPath);
            return false;
        }

        outFile << "Category Centroids\n";
        outFile << "==================\n\n";
        outFile << "Category\tCenterX\tCenterY\tPixelCount\n";

        for (const auto& kv : accum) {
            double cx = kv.second.sumX / kv.second.count;
            double cy = kv.second.sumY / kv.second.count;
            outFile << kv.first << "\t" << cx << "\t" << cy << "\t" << kv.second.count << "\n";
        }
        outFile.close();

        // Optionally create point raster
        QString rasterOut = parameter("output_raster").toString();
        if (!rasterOut.isEmpty()) {
            Raster outRaster(cols, rows, 1, DataType::Float64);
            outRaster.setGeoTransform(input->geoTransform());
            outRaster.setProjection(input->projection());
            outRaster.setNoDataValue(0.0);
            auto& rOut = outRaster.data(0);

            for (const auto& kv : accum) {
                double cx = kv.second.sumX / kv.second.count;
                double cy = kv.second.sumY / kv.second.count;
                int col, row;
                input->xyToColRow(cx, cy, col, row);
                if (col >= 0 && col < cols && row >= 0 && row < rows) {
                    rOut[static_cast<size_t>(row) * cols + col] = kv.first;
                }
            }

            GdalIO::write(outRaster, rasterOut);
        }

        reportProgress(1.0,
            QString("Computed centroids for %1 categories").arg(accum.size()));
        return true;
    }
};

REGISTER_MODULE(CenterModule)

} // namespace aplaceholder
