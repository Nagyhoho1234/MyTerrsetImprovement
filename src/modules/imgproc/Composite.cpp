#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class CompositeModule : public Module {
public:
    QString name() const override { return "COMPOSITE"; }
    QString description() const override {
        return "Color composite generation. Combines three single-band images into "
               "an RGB color composite for visual interpretation.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("red_band", "Red band image"),
            ParameterDef::file("green_band", "Green band image"),
            ParameterDef::file("blue_band", "Blue band image"),
            ParameterDef::output("output", "Output composite image"),
        };
    }

    bool execute() override {
        auto rRaster = GdalIO::read(parameter("red_band").toString());
        auto gRaster = GdalIO::read(parameter("green_band").toString());
        auto bRaster = GdalIO::read(parameter("blue_band").toString());

        if (!rRaster || !gRaster || !bRaster) {
            setError("Failed to read one or more input band rasters");
            return false;
        }

        int cols = rRaster->cols(), rows = rRaster->rows();
        if (gRaster->cols() != cols || gRaster->rows() != rows ||
            bRaster->cols() != cols || bRaster->rows() != rows) {
            setError("All input band rasters must have the same dimensions");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;

        // Create 3-band output raster
        Raster output(cols, rows, 3, DataType::Float64);
        output.setGeoTransform(rRaster->geoTransform());
        output.setProjection(rRaster->projection());
        output.setNoDataValue(rRaster->noDataValue());

        // Compute stats for auto min-max stretch per band
        auto rStats = rRaster->computeStats(0);
        auto gStats = gRaster->computeStats(0);
        auto bStats = bRaster->computeStats(0);

        const auto& rSrc = rRaster->data(0);
        const auto& gSrc = gRaster->data(0);
        const auto& bSrc = bRaster->data(0);
        auto& rDst = output.data(0);
        auto& gDst = output.data(1);
        auto& bDst = output.data(2);

        double noData = rRaster->noDataValue();
        bool hasND = rRaster->hasNoData();

        reportProgress(0.0, "Creating composite...");

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (rSrc[i] == noData || gSrc[i] == noData || bSrc[i] == noData)) {
                rDst[i] = noData;
                gDst[i] = noData;
                bDst[i] = noData;
                continue;
            }

            // Auto min-max linear stretch to [0, 255] per band
            double rRange = rStats.max - rStats.min;
            double gRange = gStats.max - gStats.min;
            double bRange = bStats.max - bStats.min;

            rDst[i] = (rRange != 0.0)
                ? std::clamp((rSrc[i] - rStats.min) / rRange * 255.0, 0.0, 255.0)
                : 127.0;
            gDst[i] = (gRange != 0.0)
                ? std::clamp((gSrc[i] - gStats.min) / gRange * 255.0, 0.0, 255.0)
                : 127.0;
            bDst[i] = (bRange != 0.0)
                ? std::clamp((bSrc[i] - bStats.min) / bRange * 255.0, 0.0, 255.0)
                : 127.0;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(CompositeModule)

} // namespace aplaceholder
