#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class CVAModule : public Module {
public:
    QString name() const override { return "CVA"; }
    QString description() const override {
        return "Change Vector Analysis. Computes magnitude and direction of change "
               "between two multi-band images. Magnitude is the Euclidean distance "
               "in spectral space; direction is computed from the first two bands.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time1_bands", "Time 1 band rasters (comma-separated paths)"),
            ParameterDef::file("time2_bands", "Time 2 band rasters (comma-separated paths)"),
            ParameterDef::output("output_magnitude", "Output change magnitude raster"),
            ParameterDef::output("output_direction", "Output change direction raster (degrees)"),
        };
    }

    bool execute() override {
        QStringList t1Paths = parameter("time1_bands").toString().split(",", Qt::SkipEmptyParts);
        QStringList t2Paths = parameter("time2_bands").toString().split(",", Qt::SkipEmptyParts);

        if (t1Paths.size() != t2Paths.size()) {
            setError("Number of Time 1 bands must equal number of Time 2 bands");
            return false;
        }

        int numBands = t1Paths.size();
        if (numBands < 1) {
            setError("At least one band pair is required");
            return false;
        }

        reportProgress(0.0, "Loading band rasters...");

        // Load all band rasters
        std::vector<std::unique_ptr<Raster>> t1(numBands), t2(numBands);
        for (int b = 0; b < numBands; ++b) {
            t1[b] = GdalIO::read(t1Paths[b].trimmed());
            t2[b] = GdalIO::read(t2Paths[b].trimmed());
            if (!t1[b] || !t2[b]) {
                setError("Failed to read band raster at index " + QString::number(b));
                return false;
            }
        }

        // Validate dimensions
        int cols = t1[0]->cols(), rows = t1[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        for (int b = 0; b < numBands; ++b) {
            if (t1[b]->cols() != cols || t1[b]->rows() != rows ||
                t2[b]->cols() != cols || t2[b]->rows() != rows) {
                setError("All band rasters must have the same dimensions");
                return false;
            }
        }

        reportProgress(0.2, "Computing change vectors...");

        double outNoData = -9999.0;

        Raster magnitude(cols, rows, 1, DataType::Float64);
        magnitude.setGeoTransform(t1[0]->geoTransform());
        magnitude.setProjection(t1[0]->projection());
        magnitude.setNoDataValue(outNoData);
        auto& magData = magnitude.data(0);

        Raster direction(cols, rows, 1, DataType::Float64);
        direction.setGeoTransform(t1[0]->geoTransform());
        direction.setProjection(t1[0]->projection());
        direction.setNoDataValue(outNoData);
        auto& dirData = direction.data(0);

        for (int64_t i = 0; i < total; ++i) {
            bool valid = true;
            double sumSqDiff = 0.0;
            double diff0 = 0.0, diff1 = 0.0;

            for (int b = 0; b < numBands; ++b) {
                const auto& d1 = t1[b]->data(0);
                const auto& d2 = t2[b]->data(0);
                bool nd1 = t1[b]->hasNoData() && d1[i] == t1[b]->noDataValue();
                bool nd2 = t2[b]->hasNoData() && d2[i] == t2[b]->noDataValue();
                if (nd1 || nd2) { valid = false; break; }

                double diff = d2[i] - d1[i];
                sumSqDiff += diff * diff;
                if (b == 0) diff0 = diff;
                if (b == 1) diff1 = diff;
            }

            if (!valid) {
                magData[i] = outNoData;
                dirData[i] = outNoData;
            } else {
                magData[i] = std::sqrt(sumSqDiff);
                // Direction from first two bands (or single band)
                if (numBands >= 2) {
                    double angle = std::atan2(diff1, diff0) * 180.0 / M_PI;
                    if (angle < 0.0) angle += 360.0;
                    dirData[i] = angle;
                } else {
                    // Single band: direction is sign (0 or 180)
                    dirData[i] = (diff0 >= 0.0) ? 0.0 : 180.0;
                }
            }

            if (i % 1000000 == 0)
                reportProgress(0.2 + 0.6 * static_cast<double>(i) / total);
        }

        reportProgress(0.85, "Writing output rasters...");

        if (!GdalIO::write(magnitude, parameter("output_magnitude").toString())) {
            setError("Failed to write magnitude raster");
            return false;
        }

        if (!GdalIO::write(direction, parameter("output_direction").toString())) {
            setError("Failed to write direction raster");
            return false;
        }

        reportProgress(1.0,
            QString("CVA complete: %1 bands, %2 pixels")
                .arg(numBands)
                .arg(total));

        return true;
    }
};

REGISTER_MODULE(CVAModule)

} // namespace aplaceholder
