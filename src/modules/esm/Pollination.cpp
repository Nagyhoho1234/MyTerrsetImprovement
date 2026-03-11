#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <QMap>
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class PollinationModule : public Module {
public:
    QString name() const override { return "POLLINATION"; }
    QString description() const override {
        return "Pollination service modeling (InVEST-like). Models pollinator abundance "
               "based on nesting suitability and floral resources within a foraging range. "
               "For each pixel, nesting suitability determines potential pollinator supply, "
               "and floral resources within the foraging distance are aggregated to estimate "
               "pollinator abundance visiting each cell.";
    }
    QString category() const override { return "Ecosystem Services"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("landcover", "Land cover classification raster",
                "Integer raster where each pixel value is a land cover class ID"),
            ParameterDef::file("nesting_table", "Nesting suitability table (CSV)",
                "CSV with columns: class_id, nesting_suitability (0-1)"),
            ParameterDef::file("floral_table", "Floral resources table (CSV)",
                "CSV with columns: class_id, floral_resources (0-1)"),
            ParameterDef::output("output", "Output pollinator abundance raster",
                "Continuous raster of modeled pollinator abundance"),
            ParameterDef::real("foraging_distance", "Foraging distance (pixels)", 10.0, 1.0, 500.0,
                "Maximum foraging distance in pixel units. Determines the radius "
                "over which floral resources are aggregated with distance decay."),
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

        int cols = landCover->cols();
        int rows = landCover->rows();
        int64_t total = landCover->cellCount();
        double noData = landCover->noDataValue();
        bool hasND = landCover->hasNoData();

        // 2. Parse nesting suitability CSV
        QMap<int, double> nestingTable;
        if (!parseLookupCSV(parameter("nesting_table").toString(), nestingTable)) {
            return false;
        }
        reportProgress(0.08, QString("Nesting table loaded: %1 classes.").arg(nestingTable.size()));

        // 3. Parse floral resources CSV
        QMap<int, double> floralTable;
        if (!parseLookupCSV(parameter("floral_table").toString(), floralTable)) {
            return false;
        }
        reportProgress(0.10, QString("Floral table loaded: %1 classes.").arg(floralTable.size()));

        // 4. Build nesting suitability and floral resource rasters from land cover
        const auto& lcData = landCover->data(0);
        std::vector<double> nestingMap(total, 0.0);
        std::vector<double> floralMap(total, 0.0);

        for (int64_t i = 0; i < total; ++i) {
            double raw = lcData[i];
            if (hasND && raw == noData)
                continue;

            int classId = static_cast<int>(std::round(raw));

            auto nIt = nestingTable.find(classId);
            if (nIt != nestingTable.end())
                nestingMap[i] = nIt.value();

            auto fIt = floralTable.find(classId);
            if (fIt != floralTable.end())
                floralMap[i] = fIt.value();
        }

        reportProgress(0.20, "Nesting and floral maps built.");

        // 5. Compute pollinator abundance
        //    For each pixel with nesting suitability > 0, aggregate floral resources
        //    within foraging distance using exponential distance decay, then
        //    multiply by nesting suitability.
        //    Pollinator supply at source j: PS_j = nesting_j * sum_i(floral_i * exp(-d/alpha))
        //    Then sum contributions at each destination pixel.
        double forageDist = parameter("foraging_distance").toDouble();
        int radius = static_cast<int>(std::ceil(forageDist));
        double alpha = forageDist / 2.0; // decay parameter

        // First pass: compute pollinator supply at each nesting site
        std::vector<double> pollinatorSupply(total, 0.0);

        for (int r = 0; r < rows; ++r) {
            if (r % 50 == 0)
                reportProgress(0.20 + 0.30 * static_cast<double>(r) / rows,
                               "Computing pollinator supply...");

            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;

                if (hasND && lcData[idx] == noData)
                    continue;

                if (nestingMap[idx] <= 0.0)
                    continue;

                // Aggregate floral resources within foraging distance
                double floralSum = 0.0;
                double weightSum = 0.0;

                int rStart = std::max(0, r - radius);
                int rEnd = std::min(rows - 1, r + radius);
                int cStart = std::max(0, c - radius);
                int cEnd = std::min(cols - 1, c + radius);

                for (int wr = rStart; wr <= rEnd; ++wr) {
                    for (int wc = cStart; wc <= cEnd; ++wc) {
                        double dist = std::sqrt(static_cast<double>((wr - r) * (wr - r) +
                                                                     (wc - c) * (wc - c)));
                        if (dist > forageDist)
                            continue;

                        int64_t wIdx = static_cast<int64_t>(wr) * cols + wc;
                        if (hasND && lcData[wIdx] == noData)
                            continue;

                        double weight = std::exp(-dist / alpha);
                        floralSum += floralMap[wIdx] * weight;
                        weightSum += weight;
                    }
                }

                if (weightSum > 0.0)
                    pollinatorSupply[idx] = nestingMap[idx] * (floralSum / weightSum);
            }
        }

        reportProgress(0.55, "Pollinator supply computed. Distributing to landscape...");

        // Second pass: distribute pollinator visits from nesting sites to surrounding cells
        std::vector<double> pollinatorAbundance(total, 0.0);

        for (int r = 0; r < rows; ++r) {
            if (r % 50 == 0)
                reportProgress(0.55 + 0.30 * static_cast<double>(r) / rows,
                               "Distributing pollinator visits...");

            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;

                if (hasND && lcData[idx] == noData)
                    continue;

                if (pollinatorSupply[idx] <= 0.0)
                    continue;

                // Distribute supply to surrounding cells based on floral attractiveness
                int rStart = std::max(0, r - radius);
                int rEnd = std::min(rows - 1, r + radius);
                int cStart = std::max(0, c - radius);
                int cEnd = std::min(cols - 1, c + radius);

                double totalAttractiveness = 0.0;
                for (int wr = rStart; wr <= rEnd; ++wr) {
                    for (int wc = cStart; wc <= cEnd; ++wc) {
                        double dist = std::sqrt(static_cast<double>((wr - r) * (wr - r) +
                                                                     (wc - c) * (wc - c)));
                        if (dist > forageDist)
                            continue;

                        int64_t wIdx = static_cast<int64_t>(wr) * cols + wc;
                        if (hasND && lcData[wIdx] == noData)
                            continue;

                        totalAttractiveness += floralMap[wIdx] * std::exp(-dist / alpha);
                    }
                }

                if (totalAttractiveness <= 0.0)
                    continue;

                for (int wr = rStart; wr <= rEnd; ++wr) {
                    for (int wc = cStart; wc <= cEnd; ++wc) {
                        double dist = std::sqrt(static_cast<double>((wr - r) * (wr - r) +
                                                                     (wc - c) * (wc - c)));
                        if (dist > forageDist)
                            continue;

                        int64_t wIdx = static_cast<int64_t>(wr) * cols + wc;
                        if (hasND && lcData[wIdx] == noData)
                            continue;

                        double fraction = (floralMap[wIdx] * std::exp(-dist / alpha))
                                          / totalAttractiveness;
                        pollinatorAbundance[wIdx] += pollinatorSupply[idx] * fraction;
                    }
                }
            }
        }

        reportProgress(0.88, "Building output raster...");

        // 6. Create and write output
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(landCover->geoTransform());
        output.setProjection(landCover->projection());
        output.setNoDataValue(noData);

        auto& outData = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && lcData[i] == noData) {
                outData[i] = noData;
            } else {
                outData[i] = pollinatorAbundance[i];
            }
        }

        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        auto stats = output.computeStats(0);
        reportProgress(1.0,
            QString("Done. Pollinator abundance — min: %1, max: %2, mean: %3")
                .arg(stats.min, 0, 'f', 4)
                .arg(stats.max, 0, 'f', 4)
                .arg(stats.mean, 0, 'f', 4));

        return true;
    }

private:
    bool parseLookupCSV(const QString& csvPath, QMap<int, double>& table) {
        QFile csvFile(csvPath);
        if (!csvFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open CSV: " + csvPath);
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
                setError(QString("CSV line %1: expected at least 2 fields").arg(lineNum));
                return false;
            }

            bool okId = false, okVal = false;
            int classId = fields[0].trimmed().toInt(&okId);
            double value = fields[1].trimmed().toDouble(&okVal);

            if (!okId || !okVal) {
                setError(QString("CSV line %1: could not parse values").arg(lineNum));
                return false;
            }

            table[classId] = value;
        }

        return true;
    }
};

REGISTER_MODULE(PollinationModule)

} // namespace aplaceholder
