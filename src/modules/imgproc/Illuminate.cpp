#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class IlluminateModule : public Module {
public:
    QString name() const override { return "ILLUMINATE"; }
    QString description() const override {
        return "Topographic illumination correction using the Minnaert method. "
               "Corrects for sun angle effects on mountainous terrain by computing "
               "the illumination angle from a DEM and applying the Minnaert "
               "constant to normalize reflectance.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::file("dem", "Digital Elevation Model"),
            ParameterDef::output("output", "Output corrected image"),
            ParameterDef::real("sun_elevation", "Sun elevation angle (degrees)", 45.0,
                0.0, 90.0, "Sun elevation angle in degrees above the horizon"),
            ParameterDef::real("sun_azimuth", "Sun azimuth angle (degrees)", 180.0,
                0.0, 360.0, "Sun azimuth angle in degrees clockwise from north"),
            ParameterDef::real("k_coefficient", "Minnaert k coefficient", 0.5,
                0.0, 1.0, "Minnaert constant controlling the correction strength"),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString demPath = parameter("dem").toString();
        QString outputPath = parameter("output").toString();
        double sunElev = parameter("sun_elevation").toDouble();
        double sunAz = parameter("sun_azimuth").toDouble();
        double kCoeff = parameter("k_coefficient").toDouble();

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

        // Convert sun angles to radians
        double sunZenith = (90.0 - sunElev) * M_PI / 180.0;
        double sunAzRad = sunAz * M_PI / 180.0;
        double cosZ = std::cos(sunZenith);
        double sinZ = std::sin(sunZenith);

        const auto& demData = dem->data(0);
        const GeoTransform& gt = dem->geoTransform();
        double cellSizeX = std::abs(gt.pixelWidth);
        double cellSizeY = std::abs(gt.pixelHeight);

        // Compute slope and aspect from DEM
        reportProgress(0.0, "Computing slope and aspect...");
        std::vector<double> cosIllum(static_cast<size_t>(cols) * rows, cosZ);

        for (int r = 1; r < rows - 1; ++r) {
            for (int c = 1; c < cols - 1; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;

                // 3x3 Sobel-like gradient
                double z1 = demData[(r - 1) * cols + (c - 1)];
                double z2 = demData[(r - 1) * cols + c];
                double z3 = demData[(r - 1) * cols + (c + 1)];
                double z4 = demData[r * cols + (c - 1)];
                double z6 = demData[r * cols + (c + 1)];
                double z7 = demData[(r + 1) * cols + (c - 1)];
                double z8 = demData[(r + 1) * cols + c];
                double z9 = demData[(r + 1) * cols + (c + 1)];

                double dzdx = ((z3 + 2.0 * z6 + z9) - (z1 + 2.0 * z4 + z7)) / (8.0 * cellSizeX);
                double dzdy = ((z7 + 2.0 * z8 + z9) - (z1 + 2.0 * z2 + z3)) / (8.0 * cellSizeY);

                double slope = std::atan(std::sqrt(dzdx * dzdx + dzdy * dzdy));
                double aspect = std::atan2(dzdy, -dzdx);

                // Cosine of illumination angle
                cosIllum[idx] = std::cos(slope) * cosZ +
                                std::sin(slope) * sinZ * std::cos(sunAzRad - aspect);
            }
        }

        // Apply Minnaert correction per band
        Raster output(cols, rows, numBands, DataType::Float32);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(outNoData);

        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int b = 0; b < numBands; ++b) {
            reportProgress(0.2 + 0.7 * static_cast<double>(b) / numBands,
                           QString("Correcting band %1...").arg(b + 1));

            const auto& inData = input->data(b);
            auto& outData = output.data(b);

            for (int64_t i = 0; i < total; ++i) {
                double val = inData[i];
                if (hasND && val == inputND) {
                    outData[i] = outNoData;
                    continue;
                }

                double ci = cosIllum[i];
                if (ci <= 0.0) {
                    outData[i] = outNoData;
                    continue;
                }

                // Minnaert correction: L_corrected = L * (cosZ / cosIllum)^k
                double correction = std::pow(cosZ / ci, kCoeff);
                outData[i] = val * correction;
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
    }
};

REGISTER_MODULE(IlluminateModule)

} // namespace aplaceholder
