#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <cmath>

namespace aplaceholder {

class ProfileModule : public Module {
public:
    QString name() const override { return "PROFILE"; }
    QString description() const override {
        return "Extract values along a transect line from a raster image. "
               "Samples at pixel intervals and outputs CSV with distance and value columns.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::real("start_x", "Start X coordinate", 0, -999999, 999999),
            ParameterDef::real("start_y", "Start Y coordinate", 0, -999999, 999999),
            ParameterDef::real("end_x", "End X coordinate", 0, -999999, 999999),
            ParameterDef::real("end_y", "End Y coordinate", 0, -999999, 999999),
            ParameterDef::output("output_file", "Output CSV file"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        double startX = parameter("start_x").toDouble();
        double startY = parameter("start_y").toDouble();
        double endX = parameter("end_x").toDouble();
        double endY = parameter("end_y").toDouble();

        // Compute the sampling interval as the pixel size
        double pixW = std::abs(raster->geoTransform().pixelWidth);
        double pixH = std::abs(raster->geoTransform().pixelHeight);
        double pixelSize = std::min(pixW, pixH);

        double lineLen = std::sqrt((endX - startX) * (endX - startX) +
                                   (endY - startY) * (endY - startY));
        if (lineLen < 1e-12) {
            setError("Start and end points are the same");
            return false;
        }

        int numSamples = static_cast<int>(std::ceil(lineLen / pixelSize)) + 1;
        double dx = (endX - startX) / (numSamples - 1);
        double dy = (endY - startY) / (numSamples - 1);

        // Write output CSV
        QFile outFile(parameter("output_file").toString());
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file for writing");
            return false;
        }

        QTextStream out(&outFile);
        out << "distance,value\n";

        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        for (int s = 0; s < numSamples; ++s) {
            double x = startX + s * dx;
            double y = startY + s * dy;
            double distance = s * pixelSize;
            // Adjust last sample to exact line length
            if (s == numSamples - 1)
                distance = lineLen;

            int col, row;
            raster->xyToColRow(x, y, col, row);

            if (col >= 0 && col < raster->cols() && row >= 0 && row < raster->rows()) {
                double value = raster->value(col, row);
                if (hasND && value == noData) {
                    out << distance << ",NoData\n";
                } else {
                    out << distance << "," << value << "\n";
                }
            } else {
                out << distance << ",OutOfBounds\n";
            }

            if (s % 10000 == 0)
                reportProgress(static_cast<double>(s) / numSamples);
        }

        outFile.close();
        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(ProfileModule)

} // namespace aplaceholder
