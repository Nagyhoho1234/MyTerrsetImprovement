#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class VegIndexModule : public Module {
public:
    QString name() const override { return "VEGINDEX"; }
    QString description() const override {
        return "Vegetation index computation. Supports SAVI, EVI, MSAVI, GNDVI, and NDWI "
               "indices from multispectral band inputs.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("nir_band", "NIR band image"),
            ParameterDef::file("red_band", "Red band image"),
            ParameterDef::file("blue_band", "Blue band image (required for EVI)",
                "Blue band image; only required for EVI computation"),
            ParameterDef::file("green_band", "Green band image (required for GNDVI/NDWI)",
                "Green band image; required for GNDVI and NDWI"),
            ParameterDef::output("output", "Output vegetation index image"),
            ParameterDef::combo("index_type", "Vegetation index",
                {"SAVI", "EVI", "MSAVI", "GNDVI", "NDWI"}, 0,
                "Type of vegetation index to compute"),
            ParameterDef::real("savi_l", "SAVI L factor", 0.5, 0.0, 1.0,
                "Soil brightness correction factor for SAVI (0=NDVI, 0.5=intermediate, 1.0=low veg)"),
        };
    }

    bool execute() override {
        int indexType = parameter("index_type").toInt();

        // ------------------------------------------------------------------
        // 1. Read required bands
        // ------------------------------------------------------------------
        auto nirRaster = GdalIO::read(parameter("nir_band").toString());
        if (!nirRaster) {
            setError("Failed to read NIR band image");
            return false;
        }

        auto redRaster = GdalIO::read(parameter("red_band").toString());
        if (!redRaster) {
            setError("Failed to read Red band image");
            return false;
        }

        int cols = nirRaster->cols(), rows = nirRaster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = nirRaster->hasNoData();
        double noData = nirRaster->noDataValue();

        if (redRaster->cols() != cols || redRaster->rows() != rows) {
            setError("NIR and Red bands must have the same dimensions");
            return false;
        }

        // Optional bands
        std::unique_ptr<Raster> blueRaster;
        std::unique_ptr<Raster> greenRaster;

        if (indexType == 1) { // EVI
            blueRaster = GdalIO::read(parameter("blue_band").toString());
            if (!blueRaster) {
                setError("Failed to read Blue band image (required for EVI)");
                return false;
            }
            if (blueRaster->cols() != cols || blueRaster->rows() != rows) {
                setError("Blue band must have the same dimensions as NIR/Red");
                return false;
            }
        }

        if (indexType == 3 || indexType == 4) { // GNDVI or NDWI
            greenRaster = GdalIO::read(parameter("green_band").toString());
            if (!greenRaster) {
                setError("Failed to read Green band image (required for GNDVI/NDWI)");
                return false;
            }
            if (greenRaster->cols() != cols || greenRaster->rows() != rows) {
                setError("Green band must have the same dimensions as NIR/Red");
                return false;
            }
        }

        double saviL = parameter("savi_l").toDouble();

        reportProgress(0.1, "Computing vegetation index...");

        // ------------------------------------------------------------------
        // 2. Create output
        // ------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(nirRaster->geoTransform());
        output.setProjection(nirRaster->projection());
        output.setNoDataValue(noData);

        const auto& nir = nirRaster->data(0);
        const auto& red = redRaster->data(0);
        auto& out = output.data(0);

        const std::vector<double>* blue = blueRaster ? &blueRaster->data(0) : nullptr;
        const std::vector<double>* green = greenRaster ? &greenRaster->data(0) : nullptr;

        // ------------------------------------------------------------------
        // 3. Compute index per pixel
        // ------------------------------------------------------------------
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && nir[i] == noData) {
                out[i] = noData;
                continue;
            }

            double n = nir[i];
            double r = red[i];

            switch (indexType) {
            case 0: {
                // SAVI = ((NIR - Red) / (NIR + Red + L)) * (1 + L)
                double denom = n + r + saviL;
                if (std::abs(denom) < 1e-10)
                    out[i] = 0.0;
                else
                    out[i] = ((n - r) / denom) * (1.0 + saviL);
                break;
            }
            case 1: {
                // EVI = G * (NIR - Red) / (NIR + C1*Red - C2*Blue + L)
                // Default MODIS coefficients: G=2.5, C1=6, C2=7.5, L=1
                double b = (*blue)[i];
                double denom = n + 6.0 * r - 7.5 * b + 1.0;
                if (std::abs(denom) < 1e-10)
                    out[i] = 0.0;
                else
                    out[i] = 2.5 * (n - r) / denom;
                break;
            }
            case 2: {
                // MSAVI2 = (2*NIR + 1 - sqrt((2*NIR+1)^2 - 8*(NIR-Red))) / 2
                double inner = (2.0 * n + 1.0);
                double disc = inner * inner - 8.0 * (n - r);
                if (disc < 0.0)
                    out[i] = 0.0;
                else
                    out[i] = (inner - std::sqrt(disc)) / 2.0;
                break;
            }
            case 3: {
                // GNDVI = (NIR - Green) / (NIR + Green)
                double g = (*green)[i];
                double denom = n + g;
                if (std::abs(denom) < 1e-10)
                    out[i] = 0.0;
                else
                    out[i] = (n - g) / denom;
                break;
            }
            case 4: {
                // NDWI = (Green - NIR) / (Green + NIR)
                double g = (*green)[i];
                double denom = g + n;
                if (std::abs(denom) < 1e-10)
                    out[i] = 0.0;
                else
                    out[i] = (g - n) / denom;
                break;
            }
            default:
                out[i] = 0.0;
                break;
            }

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(VegIndexModule)

} // namespace aplaceholder
