#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class OrthoModule : public Module {
public:
    QString name() const override { return "ORTHO"; }
    QString description() const override {
        return "Simple orthorectification. Removes relief displacement from an image "
               "using a DEM, sensor height, and focal length. Corrects for parallax "
               "caused by terrain elevation relative to nadir.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::file("dem", "Digital Elevation Model"),
            ParameterDef::output("output", "Output orthorectified image"),
            ParameterDef::real("sensor_height", "Sensor height (meters)", 700000.0,
                1.0, 999999999.0, "Height of the sensor above the reference datum"),
            ParameterDef::real("focal_length", "Focal length (meters)", 0.15,
                0.001, 999999.0, "Camera/sensor focal length in meters"),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString demPath = parameter("dem").toString();
        QString outputPath = parameter("output").toString();
        double sensorH = parameter("sensor_height").toDouble();
        double focalLen = parameter("focal_length").toDouble();

        auto input = GdalIO::read(inputPath);
        auto dem = GdalIO::read(demPath);
        if (!input || !dem) {
            setError("Failed to read input image or DEM.");
            return false;
        }

        int cols = input->cols();
        int rows = input->rows();
        int numBands = input->bands();

        if (dem->cols() != cols || dem->rows() != rows) {
            setError("Input image and DEM must have the same dimensions.");
            return false;
        }

        bool hasND = input->hasNoData();
        double inputND = input->noDataValue();
        double outNoData = -9999.0;

        const auto& demData = dem->data(0);
        const GeoTransform& gt = input->geoTransform();

        // Image center (nadir point)
        double centerCol = cols / 2.0;
        double centerRow = rows / 2.0;

        Raster output(cols, rows, numBands, DataType::Float32);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(outNoData);

        reportProgress(0.0, "Orthorectifying...");

        // For each output pixel, compute the source pixel by removing relief displacement
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                double elev = demData[idx];

                // Relief displacement: dr = h * r_dist / H
                // where h = terrain elevation, r_dist = radial distance from nadir, H = sensor height
                double dx = c - centerCol;
                double dy = r - centerRow;
                double rDist = std::sqrt(dx * dx + dy * dy);

                double displacement = 0.0;
                if (rDist > 0.0 && sensorH > elev) {
                    displacement = elev * rDist / sensorH;
                }

                // Source pixel (remove displacement in radial direction)
                double srcCol, srcRow;
                if (rDist > 0.0) {
                    double scale = displacement / rDist;
                    srcCol = c + dx * scale;
                    srcRow = r + dy * scale;
                } else {
                    srcCol = c;
                    srcRow = r;
                }

                // Nearest neighbor resampling
                int sc = static_cast<int>(std::round(srcCol));
                int sr = static_cast<int>(std::round(srcRow));

                for (int b = 0; b < numBands; ++b) {
                    auto& outData = output.data(b);
                    if (sc < 0 || sc >= cols || sr < 0 || sr >= rows) {
                        outData[idx] = outNoData;
                    } else {
                        int64_t srcIdx = static_cast<int64_t>(sr) * cols + sc;
                        double val = input->data(b)[srcIdx];
                        if (hasND && val == inputND)
                            outData[idx] = outNoData;
                        else
                            outData[idx] = val;
                    }
                }
            }

            if (r % 500 == 0)
                reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
    }
};

REGISTER_MODULE(OrthoModule)

} // namespace aplaceholder
