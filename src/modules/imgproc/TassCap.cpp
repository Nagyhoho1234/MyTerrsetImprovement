#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <array>

namespace aplaceholder {

class TassCapModule : public Module {
public:
    QString name() const override { return "TASSCAP"; }
    QString description() const override {
        return "Tasseled Cap Transformation. Applies the Kauth-Thomas transformation "
               "to Landsat multispectral imagery, producing Brightness, Greenness, "
               "and Wetness components.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated, 6 bands)",
                "Comma-separated list of 6 Landsat band raster file paths "
                "(TM/ETM+: bands 1-5,7; OLI: bands 2-7)"),
            ParameterDef::output("prefix", "Output filename prefix",
                "Output rasters will be named prefix_brightness, prefix_greenness, prefix_wetness"),
            ParameterDef::combo("sensor", "Sensor type",
                {"TM", "ETM+", "OLI"}, 0,
                "Landsat sensor determining the transformation coefficients"),
        };
    }

    bool execute() override {
        // Tasseled Cap coefficients for 6 bands: [component][band]
        // Brightness, Greenness, Wetness

        // Landsat TM coefficients (Crist 1985)
        static const double TM_COEFFS[3][6] = {
            { 0.2043,  0.4158,  0.5524,  0.5741,  0.3124,  0.2303},  // Brightness
            {-0.1603, -0.2819, -0.4934,  0.7940, -0.0002, -0.1446},  // Greenness
            { 0.0315,  0.2021,  0.3102,  0.1594, -0.6806, -0.6109},  // Wetness
        };

        // Landsat ETM+ coefficients (Huang et al. 2002)
        static const double ETM_COEFFS[3][6] = {
            { 0.3561,  0.3972,  0.3904,  0.6966,  0.2286,  0.1596},  // Brightness
            {-0.3344, -0.3544, -0.4556,  0.6966, -0.0242, -0.2630},  // Greenness
            { 0.2626,  0.2141,  0.0926,  0.0656, -0.7629, -0.5388},  // Wetness
        };

        // Landsat OLI coefficients (Baig et al. 2014)
        static const double OLI_COEFFS[3][6] = {
            { 0.3029,  0.2786,  0.4733,  0.5599,  0.5080,  0.1872},  // Brightness
            {-0.2941, -0.2430, -0.5424,  0.7276,  0.0713, -0.1608},  // Greenness
            { 0.1511,  0.1973,  0.3283,  0.3407, -0.7117, -0.4559},  // Wetness
        };

        static const char* COMPONENT_NAMES[3] = {"brightness", "greenness", "wetness"};

        // Parse band file paths
        QString bandsParam = parameter("bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths)
            p = p.trimmed();

        if (bandPaths.size() != 6) {
            setError(QString("Tasseled Cap requires exactly 6 input bands, got %1")
                     .arg(bandPaths.size()));
            return false;
        }

        // Select coefficients based on sensor type
        int sensorIdx = parameter("sensor").toInt();
        const double (*coeffs)[6] = nullptr;
        QString sensorName;

        switch (sensorIdx) {
        case 0:
            coeffs = TM_COEFFS;
            sensorName = "TM";
            break;
        case 1:
            coeffs = ETM_COEFFS;
            sensorName = "ETM+";
            break;
        case 2:
            coeffs = OLI_COEFFS;
            sensorName = "OLI";
            break;
        default:
            setError("Invalid sensor type selection");
            return false;
        }

        // Read band rasters
        reportProgress(0.0, "Reading input bands...");
        std::vector<std::unique_ptr<Raster>> bandRasters(6);
        int cols = 0, rows = 0;

        for (int b = 0; b < 6; ++b) {
            bandRasters[b] = GdalIO::read(bandPaths[b]);
            if (!bandRasters[b]) {
                setError("Failed to read band image: " + bandPaths[b]);
                return false;
            }
            if (b == 0) {
                cols = bandRasters[b]->cols();
                rows = bandRasters[b]->rows();
            } else {
                if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                    setError("All band images must have the same dimensions. "
                             "Mismatch in: " + bandPaths[b]);
                    return false;
                }
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;

        // Collect band data pointers
        std::vector<const std::vector<double>*> bands(6);
        for (int b = 0; b < 6; ++b)
            bands[b] = &bandRasters[b]->data(0);

        bool hasNoData = bandRasters[0]->hasNoData();
        double noData = bandRasters[0]->noDataValue();

        QString prefix = parameter("prefix").toString();

        // Compute and write each component
        for (int comp = 0; comp < 3; ++comp) {
            reportProgress(0.1 + 0.25 * comp,
                           QString("Computing %1...").arg(COMPONENT_NAMES[comp]));

            Raster output(cols, rows, 1, DataType::Float32);
            output.setGeoTransform(bandRasters[0]->geoTransform());
            output.setProjection(bandRasters[0]->projection());
            double outNoData = -9999.0;
            output.setNoDataValue(outNoData);
            auto& out = output.data(0);

            for (int64_t i = 0; i < total; ++i) {
                // Check nodata on first band as proxy
                if (hasNoData && (*bands[0])[i] == noData) {
                    out[i] = outNoData;
                    continue;
                }

                double result = 0.0;
                for (int b = 0; b < 6; ++b)
                    result += coeffs[comp][b] * (*bands[b])[i];

                out[i] = result;

                if (i % 1000000 == 0)
                    reportProgress(0.1 + 0.25 * comp +
                                   0.25 * static_cast<double>(i) / total);
            }

            // Write output raster
            QString outPath = prefix + "_" + COMPONENT_NAMES[comp];
            if (!GdalIO::write(output, outPath)) {
                setError("Failed to write output raster: " + outPath);
                return false;
            }
        }

        reportProgress(1.0, QString("Tasseled Cap transformation complete (%1 sensor). "
                                     "Output: %2_brightness, %2_greenness, %2_wetness")
                        .arg(sensorName).arg(prefix));
        return true;
    }
};

REGISTER_MODULE(TassCapModule)

} // namespace aplaceholder
