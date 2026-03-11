#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <QMap>
#include <cmath>

namespace aplaceholder {

class CarbonModule : public Module {
public:
    QString name() const override { return "CARBON_STOCK"; }
    QString description() const override {
        return "Carbon stock estimation module. Maps land cover classes to carbon "
               "density values (tons C/ha) from a lookup CSV, then multiplies by "
               "pixel area to produce a per-pixel carbon stock raster (tons C).";
    }
    QString category() const override { return "Ecosystem Services"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("landcover", "Land cover classification raster",
                "Integer raster where each pixel value is a land cover class ID"),
            ParameterDef::file("carbon_table", "Carbon density table (CSV)",
                "CSV with columns: class_id, carbon_tons_per_ha"),
            ParameterDef::output("output", "Output carbon stock raster",
                "Continuous raster where each pixel = carbon stock in tons C"),
            ParameterDef::real("pixel_area_ha", "Pixel area override (ha)", 0.0, 0.0, 1e12,
                "Override pixel area in hectares. If 0, area is computed from "
                "the raster geo-transform"),
        };
    }

    bool execute() override {
        // 1. Read land cover raster
        QString lcPath = parameter("landcover").toString();
        auto landCover = GdalIO::read(lcPath);
        if (!landCover) {
            setError("Failed to read land cover raster: " + lcPath);
            return false;
        }

        reportProgress(0.05, "Land cover raster loaded.");

        // 2. Parse the carbon density CSV
        //    Expected format: class_id, carbon_tons_per_ha
        //    First row is a header and is skipped.
        QString csvPath = parameter("carbon_table").toString();
        QMap<int, double> carbonTable;

        {
            QFile csvFile(csvPath);
            if (!csvFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
                setError("Failed to open carbon table: " + csvPath);
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
                if (fields.size() < 2) {
                    setError(QString("Carbon CSV line %1: expected at least 2 "
                                     "comma-separated fields (class_id, carbon_tons_per_ha)")
                                 .arg(lineNum));
                    return false;
                }

                bool okId = false, okVal = false;
                int classId = fields[0].trimmed().toInt(&okId);
                double carbonPerHa = fields[1].trimmed().toDouble(&okVal);

                if (!okId || !okVal) {
                    setError(QString("Carbon CSV line %1: could not parse "
                                     "class_id or carbon_tons_per_ha").arg(lineNum));
                    return false;
                }

                carbonTable[classId] = carbonPerHa;
            }
        }

        if (carbonTable.isEmpty()) {
            setError("Carbon table is empty or could not be parsed.");
            return false;
        }

        reportProgress(0.10, QString("Carbon table loaded: %1 classes.")
                        .arg(carbonTable.size()));

        // 3. Determine pixel area in hectares
        double pixelAreaHa = parameter("pixel_area_ha").toDouble();

        if (pixelAreaHa <= 0.0) {
            const GeoTransform& gt = landCover->geoTransform();
            double pw = std::abs(gt.pixelWidth);
            double ph = std::abs(gt.pixelHeight);
            double pixelAreaMapUnits = pw * ph;

            if (pw < 1.0 && ph < 1.0) {
                double metersPerDegree = 111320.0;
                double pixelAreaM2 = pixelAreaMapUnits
                                     * metersPerDegree * metersPerDegree;
                pixelAreaHa = pixelAreaM2 / 10000.0;
            } else {
                pixelAreaHa = pixelAreaMapUnits / 10000.0;
            }

            if (pixelAreaHa <= 0.0) {
                setError("Could not determine pixel area from geo-transform. "
                         "Please provide a pixel_area_ha override.");
                return false;
            }
        }

        reportProgress(0.15, QString("Pixel area: %1 ha").arg(pixelAreaHa));

        // 4. Create output raster and compute per-pixel carbon stock
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

        double totalCarbon = 0.0;
        int64_t valuedPixels = 0;
        int64_t unmappedPixels = 0;

        for (int64_t i = 0; i < total; ++i) {
            double raw = lcData[i];

            if (hasND && raw == noData) {
                outData[i] = noData;
                continue;
            }

            int classId = static_cast<int>(std::round(raw));

            auto it = carbonTable.find(classId);
            if (it == carbonTable.end()) {
                outData[i] = 0.0;
                ++unmappedPixels;
            } else {
                double carbonStock = it.value() * pixelAreaHa;
                outData[i] = carbonStock;
                totalCarbon += carbonStock;
                ++valuedPixels;
            }

            if (i % 500000 == 0)
                reportProgress(0.15 + 0.75 * static_cast<double>(i) / total);
        }

        reportProgress(0.90, "Computation complete. Writing output...");

        // 5. Write output raster
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        reportProgress(1.0,
            QString("Done. Total carbon stock: %1 tons C over %2 "
                    "valued pixels (%3 pixels had unmapped classes).")
                .arg(totalCarbon, 0, 'f', 2)
                .arg(valuedPixels)
                .arg(unmappedPixels));

        return true;
    }
};

REGISTER_MODULE(CarbonModule)

} // namespace aplaceholder
