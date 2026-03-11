#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ResultantModule : public Module {
public:
    QString name() const override { return "RESULTANT"; }
    QString description() const override {
        return "Combines two force/friction vector image pairs into a single resultant "
               "vector pair using vector addition. Directions are in degrees from north.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("mag1", "Magnitude image 1"),
            ParameterDef::file("dir1", "Direction image 1",
                "Direction in degrees from north"),
            ParameterDef::file("mag2", "Magnitude image 2"),
            ParameterDef::file("dir2", "Direction image 2",
                "Direction in degrees from north"),
            ParameterDef::output("output_magnitude", "Output resultant magnitude image"),
            ParameterDef::output("output_direction", "Output resultant direction image"),
        };
    }

    bool execute() override {
        auto rMag1 = GdalIO::read(parameter("mag1").toString());
        auto rDir1 = GdalIO::read(parameter("dir1").toString());
        auto rMag2 = GdalIO::read(parameter("mag2").toString());
        auto rDir2 = GdalIO::read(parameter("dir2").toString());

        if (!rMag1 || !rDir1 || !rMag2 || !rDir2) {
            setError("Failed to read one or more input rasters");
            return false;
        }

        int cols = rMag1->cols(), rows = rMag1->rows();
        if (rDir1->cols() != cols || rDir1->rows() != rows ||
            rMag2->cols() != cols || rMag2->rows() != rows ||
            rDir2->cols() != cols || rDir2->rows() != rows) {
            setError("All input rasters must have the same dimensions");
            return false;
        }

        Raster outMag(cols, rows, 1, DataType::Float64);
        Raster outDir(cols, rows, 1, DataType::Float64);
        outMag.setGeoTransform(rMag1->geoTransform());
        outMag.setProjection(rMag1->projection());
        outMag.setNoDataValue(rMag1->noDataValue());
        outDir.setGeoTransform(rMag1->geoTransform());
        outDir.setProjection(rMag1->projection());
        outDir.setNoDataValue(rMag1->noDataValue());

        const auto& dM1 = rMag1->data(0);
        const auto& dD1 = rDir1->data(0);
        const auto& dM2 = rMag2->data(0);
        const auto& dD2 = rDir2->data(0);
        auto& oMag = outMag.data(0);
        auto& oDir = outDir.data(0);
        double noData = rMag1->noDataValue();
        bool hasND = rMag1->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        constexpr double DEG2RAD = M_PI / 180.0;
        constexpr double RAD2DEG = 180.0 / M_PI;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (dM1[i] == noData || dD1[i] == noData ||
                          dM2[i] == noData || dD2[i] == noData)) {
                oMag[i] = noData;
                oDir[i] = noData;
                continue;
            }

            // Convert direction (degrees from north, clockwise) to math angle (radians)
            double angle1 = (90.0 - dD1[i]) * DEG2RAD;
            double angle2 = (90.0 - dD2[i]) * DEG2RAD;

            // Decompose to x,y components
            double x1 = dM1[i] * std::cos(angle1);
            double y1 = dM1[i] * std::sin(angle1);
            double x2 = dM2[i] * std::cos(angle2);
            double y2 = dM2[i] * std::sin(angle2);

            // Add components
            double rx = x1 + x2;
            double ry = y1 + y2;

            // Convert back to magnitude and direction
            oMag[i] = std::sqrt(rx * rx + ry * ry);

            // Convert math angle back to azimuth (degrees from north, clockwise)
            double azimuth = 90.0 - std::atan2(ry, rx) * RAD2DEG;
            if (azimuth < 0.0) azimuth += 360.0;
            if (azimuth >= 360.0) azimuth -= 360.0;
            oDir[i] = azimuth;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");
        bool ok = GdalIO::write(outMag, parameter("output_magnitude").toString());
        ok = ok && GdalIO::write(outDir, parameter("output_direction").toString());

        if (!ok) {
            setError("Failed to write one or more output rasters");
            return false;
        }

        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(ResultantModule)

} // namespace aplaceholder
