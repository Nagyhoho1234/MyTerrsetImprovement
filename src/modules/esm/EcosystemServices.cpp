#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <QMap>
#include <cmath>

namespace aplaceholder {

class EcosystemServicesModule : public Module {
public:
    QString name() const override { return "ECOSYSTEM_SERVICES"; }
    QString description() const override {
        return "Ecosystem services valuation based on land cover classification. "
               "Maps each land cover class to an economic value ($/ha/year) using a "
               "Costanza et al. style valuation table (CSV), then computes per-pixel "
               "ecosystem service value scaled by pixel area. Produces a continuous "
               "raster of monetary values and reports total study-area value.";
    }
    QString category() const override { return "Ecosystem Services"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("land_cover", "Land cover classification raster",
                "Integer raster where each pixel value is a land cover class ID"),
            ParameterDef::file("valuation_table", "Valuation table (CSV)",
                "CSV with columns: class_id, class_name, value_per_ha. "
                "Values in $/ha/year (Costanza et al. approach)"),
            ParameterDef::output("output", "Output ecosystem services value raster",
                "Continuous raster where each pixel = ecosystem service value in $/pixel/year"),
            ParameterDef::real("pixel_area_ha", "Pixel area override (ha)", 0.0, 0.0, 1e12,
                "Override pixel area in hectares. If 0, area is computed from "
                "the raster geo-transform (pixel width * pixel height converted to ha)"),
        };
    }

    bool execute() override {
        // ----------------------------------------------------------
        // 1. Read the land cover raster
        // ----------------------------------------------------------
        QString lcPath = parameter("land_cover").toString();
        auto landCover = GdalIO::read(lcPath);
        if (!landCover) {
            setError("Failed to read land cover raster: " + lcPath);
            return false;
        }

        reportProgress(0.05, "Land cover raster loaded.");

        // ----------------------------------------------------------
        // 2. Parse the valuation CSV
        //    Expected format: class_id,class_name,value_per_ha
        //    First row is a header and is skipped.
        // ----------------------------------------------------------
        QString csvPath = parameter("valuation_table").toString();
        QMap<int, double> valueTable;          // class_id -> $/ha/year

        {
            QFile csvFile(csvPath);
            if (!csvFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
                setError("Failed to open valuation table: " + csvPath);
                return false;
            }

            QTextStream in(&csvFile);
            bool headerSkipped = false;
            int lineNum = 0;

            while (!in.atEnd()) {
                QString line = in.readLine().trimmed();
                ++lineNum;

                if (line.isEmpty() || line.startsWith('#'))
                    continue;

                if (!headerSkipped) {
                    headerSkipped = true;
                    continue;
                }

                QStringList fields = line.split(',');
                if (fields.size() < 3) {
                    setError(QString("Valuation CSV line %1: expected at least 3 "
                                     "comma-separated fields (class_id, class_name, "
                                     "value_per_ha)").arg(lineNum));
                    return false;
                }

                bool okId = false, okVal = false;
                int classId = fields[0].trimmed().toInt(&okId);
                double valPerHa = fields[2].trimmed().toDouble(&okVal);

                if (!okId || !okVal) {
                    setError(QString("Valuation CSV line %1: could not parse "
                                     "class_id or value_per_ha").arg(lineNum));
                    return false;
                }

                valueTable[classId] = valPerHa;
            }
        }

        if (valueTable.isEmpty()) {
            setError("Valuation table is empty or could not be parsed.");
            return false;
        }

        reportProgress(0.10, QString("Valuation table loaded: %1 classes.")
                        .arg(valueTable.size()));

        // ----------------------------------------------------------
        // 3. Determine pixel area in hectares
        // ----------------------------------------------------------
        double pixelAreaHa = parameter("pixel_area_ha").toDouble();

        if (pixelAreaHa <= 0.0) {
            // Derive from geo-transform. Pixel dimensions are in map units;
            // if the CRS is geographic (degrees) this will be approximate.
            // For projected CRS in metres, 1 ha = 10 000 m^2.
            const GeoTransform& gt = landCover->geoTransform();
            double pw = std::abs(gt.pixelWidth);
            double ph = std::abs(gt.pixelHeight);
            double pixelAreaMapUnits = pw * ph;

            // Heuristic: if both pixel dimensions are < 1 the units are
            // likely degrees (geographic CRS).  Convert assuming equatorial
            // approximation: 1 degree ~ 111 320 m.
            if (pw < 1.0 && ph < 1.0) {
                double metersPerDegree = 111320.0;
                double pixelAreaM2 = pixelAreaMapUnits
                                     * metersPerDegree * metersPerDegree;
                pixelAreaHa = pixelAreaM2 / 10000.0;
            } else {
                // Assume map units are metres
                pixelAreaHa = pixelAreaMapUnits / 10000.0;
            }

            if (pixelAreaHa <= 0.0) {
                setError("Could not determine pixel area from geo-transform. "
                         "Please provide a pixel_area_ha override.");
                return false;
            }
        }

        reportProgress(0.15, QString("Pixel area: %1 ha").arg(pixelAreaHa));

        // ----------------------------------------------------------
        // 4. Create output raster and compute per-pixel values
        // ----------------------------------------------------------
        int cols = landCover->cols();
        int rows = landCover->rows();
        int64_t total = landCover->cellCount();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(landCover->geoTransform());
        output.setProjection(landCover->projection());
        output.setNoDataValue(landCover->noDataValue());

        const auto& lcData = landCover->data(0);
        auto& outData = output.data(0);

        double noData = landCover->noDataValue();
        bool hasND = landCover->hasNoData();

        double totalValue = 0.0;
        int64_t valuedPixels = 0;
        int64_t unmappedPixels = 0;

        for (int64_t i = 0; i < total; ++i) {
            double raw = lcData[i];

            // Skip NoData pixels
            if (hasND && raw == noData) {
                outData[i] = noData;
                continue;
            }

            int classId = static_cast<int>(std::round(raw));

            auto it = valueTable.find(classId);
            if (it == valueTable.end()) {
                // Class not in valuation table — assign 0 value
                outData[i] = 0.0;
                ++unmappedPixels;
            } else {
                double pixelValue = it.value() * pixelAreaHa;
                outData[i] = pixelValue;
                totalValue += pixelValue;
                ++valuedPixels;
            }

            if (i % 500000 == 0)
                reportProgress(0.15 + 0.75 * static_cast<double>(i) / total);
        }

        reportProgress(0.90, "Computation complete. Writing output...");

        // ----------------------------------------------------------
        // 5. Write output raster
        // ----------------------------------------------------------
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        reportProgress(1.0,
            QString("Done. Total ecosystem service value: $%1/year over %2 "
                    "valued pixels (%3 pixels had unmapped classes).")
                .arg(totalValue, 0, 'f', 2)
                .arg(valuedPixels)
                .arg(unmappedPixels));

        return true;
    }
};

REGISTER_MODULE(EcosystemServicesModule)

} // namespace aplaceholder
