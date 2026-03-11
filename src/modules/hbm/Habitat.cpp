#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <QMap>
#include <cmath>
#include <vector>
#include <memory>
#include <algorithm>

namespace aplaceholder {

class HabitatModule : public Module {
public:
    QString name() const override { return "HABITAT_QUALITY"; }
    QString description() const override {
        return "Habitat quality modeling (InVEST-like). Evaluates habitat quality based on "
               "proximity to threat sources and habitat sensitivity. For each pixel, computes "
               "a degradation score from distance-weighted threat intensities modulated by "
               "habitat sensitivity, then derives quality as Q = H * (1 - D^z / (D^z + k^z)) "
               "where H is habitat suitability, D is total degradation, and k,z are parameters.";
    }
    QString category() const override { return "Habitat & Biodiversity"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("landcover", "Land cover classification raster",
                "Integer raster of land cover class IDs"),
            ParameterDef::file("threats", "Threat rasters (comma-separated)",
                "Comma-separated paths to binary/continuous threat rasters. "
                "Order must match threat_id in threat_table."),
            ParameterDef::file("threat_table", "Threat attributes table (CSV)",
                "CSV with columns: threat_id, max_dist, weight, decay. "
                "threat_id is 0-based index matching the threats list order. "
                "decay: linear or exponential."),
            ParameterDef::file("sensitivity_table", "Habitat sensitivity table (CSV)",
                "CSV with columns: class_id, habitat (0-1 suitability), "
                "followed by sensitivity to each threat (sens_0, sens_1, ...)."),
            ParameterDef::output("output", "Output habitat quality raster",
                "Continuous raster [0-1] where 1 = highest quality habitat"),
        };
    }

    bool execute() override {
        // 1. Read land cover
        auto landCover = GdalIO::read(parameter("landcover").toString());
        if (!landCover) {
            setError("Failed to read land cover raster.");
            return false;
        }

        int cols = landCover->cols();
        int rows = landCover->rows();
        int64_t total = landCover->cellCount();
        double noData = landCover->noDataValue();
        bool hasND = landCover->hasNoData();

        reportProgress(0.02, "Land cover loaded.");

        // 2. Load threat rasters
        QString threatsStr = parameter("threats").toString();
        QStringList threatPaths = threatsStr.split(',', Qt::SkipEmptyParts);
        for (int i = 0; i < threatPaths.size(); ++i)
            threatPaths[i] = threatPaths[i].trimmed();

        int numThreats = threatPaths.size();
        if (numThreats == 0) {
            setError("No threat rasters provided.");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> threatRasters;
        for (int t = 0; t < numThreats; ++t) {
            auto tr = GdalIO::read(threatPaths[t]);
            if (!tr) {
                setError("Failed to read threat raster: " + threatPaths[t]);
                return false;
            }
            if (tr->cols() != cols || tr->rows() != rows) {
                setError(QString("Threat raster %1 dimension mismatch").arg(t));
                return false;
            }
            threatRasters.push_back(std::move(tr));
        }

        reportProgress(0.08, QString("%1 threat rasters loaded.").arg(numThreats));

        // 3. Parse threat table CSV: threat_id, max_dist, weight, decay
        struct ThreatInfo {
            int maxDist;        // in pixels
            double weight;
            bool exponential;   // true = exponential decay, false = linear
        };
        std::vector<ThreatInfo> threatInfos(numThreats);

        {
            QFile csvFile(parameter("threat_table").toString());
            if (!csvFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
                setError("Failed to open threat table.");
                return false;
            }
            QTextStream in(&csvFile);
            bool headerSkipped = false;
            while (!in.atEnd()) {
                QString line = in.readLine().trimmed();
                if (line.isEmpty() || line.startsWith('#')) continue;
                if (!headerSkipped) { headerSkipped = true; continue; }

                QStringList f = line.split(',');
                if (f.size() < 4) continue;

                int tid = f[0].trimmed().toInt();
                if (tid < 0 || tid >= numThreats) continue;

                threatInfos[tid].maxDist = f[1].trimmed().toInt();
                threatInfos[tid].weight = f[2].trimmed().toDouble();
                threatInfos[tid].exponential =
                    (f[3].trimmed().toLower() == "exponential");
            }
        }

        reportProgress(0.10, "Threat table parsed.");

        // 4. Parse sensitivity table CSV: class_id, habitat, sens_0, sens_1, ...
        struct SensInfo {
            double habitat;
            std::vector<double> sensitivity;  // per threat
        };
        QMap<int, SensInfo> sensTable;

        {
            QFile csvFile(parameter("sensitivity_table").toString());
            if (!csvFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
                setError("Failed to open sensitivity table.");
                return false;
            }
            QTextStream in(&csvFile);
            bool headerSkipped = false;
            while (!in.atEnd()) {
                QString line = in.readLine().trimmed();
                if (line.isEmpty() || line.startsWith('#')) continue;
                if (!headerSkipped) { headerSkipped = true; continue; }

                QStringList f = line.split(',');
                if (f.size() < 2 + numThreats) continue;

                int classId = f[0].trimmed().toInt();
                SensInfo si;
                si.habitat = f[1].trimmed().toDouble();
                for (int t = 0; t < numThreats; ++t)
                    si.sensitivity.push_back(f[2 + t].trimmed().toDouble());
                sensTable[classId] = si;
            }
        }

        reportProgress(0.15, QString("Sensitivity table: %1 classes.").arg(sensTable.size()));

        // 5. Compute total degradation per pixel
        //    D_x = sum_r( sum_y( w_r * r_y * i_rxy * S_xr ) )
        //    where i_rxy = distance decay function
        const double k = 0.5;   // half-saturation constant
        const double z = 2.5;   // scaling parameter

        std::vector<double> degradation(total, 0.0);
        const auto& lcData = landCover->data(0);

        for (int t = 0; t < numThreats; ++t) {
            reportProgress(0.15 + 0.55 * static_cast<double>(t) / numThreats,
                           QString("Processing threat %1/%2...").arg(t + 1).arg(numThreats));

            const auto& tData = threatRasters[t]->data(0);
            int maxDist = threatInfos[t].maxDist;
            double weight = threatInfos[t].weight;
            bool expDecay = threatInfos[t].exponential;

            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    int64_t idx = static_cast<int64_t>(r) * cols + c;

                    if (hasND && lcData[idx] == noData)
                        continue;

                    // Look up sensitivity of this pixel's class to this threat
                    int classId = static_cast<int>(std::round(lcData[idx]));
                    auto sensIt = sensTable.find(classId);
                    if (sensIt == sensTable.end() || sensIt.value().habitat <= 0.0)
                        continue;

                    double sens = (t < static_cast<int>(sensIt.value().sensitivity.size()))
                                      ? sensIt.value().sensitivity[t]
                                      : 0.0;
                    if (sens <= 0.0)
                        continue;

                    // Sum threat influence from surrounding pixels
                    int rStart = std::max(0, r - maxDist);
                    int rEnd = std::min(rows - 1, r + maxDist);
                    int cStart = std::max(0, c - maxDist);
                    int cEnd = std::min(cols - 1, c + maxDist);

                    double threatSum = 0.0;
                    for (int wr = rStart; wr <= rEnd; ++wr) {
                        for (int wc = cStart; wc <= cEnd; ++wc) {
                            int64_t wIdx = static_cast<int64_t>(wr) * cols + wc;
                            double tVal = tData[wIdx];
                            if (tVal <= 0.0)
                                continue;

                            double dist = std::sqrt(static_cast<double>((wr - r) * (wr - r) +
                                                                         (wc - c) * (wc - c)));
                            if (dist > maxDist)
                                continue;

                            double decay;
                            if (expDecay) {
                                decay = std::exp(-2.99 * dist / maxDist);
                            } else {
                                decay = 1.0 - (dist / maxDist);
                            }

                            threatSum += tVal * decay;
                        }
                    }

                    degradation[idx] += weight * sens * threatSum;
                }
            }
        }

        reportProgress(0.75, "Computing habitat quality...");

        // 6. Compute habitat quality: Q = H * (1 - D^z / (D^z + k^z))
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(landCover->geoTransform());
        output.setProjection(landCover->projection());
        output.setNoDataValue(noData);

        auto& outData = output.data(0);
        int64_t validCount = 0;
        double sumQ = 0.0;

        double kz = std::pow(k, z);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && lcData[i] == noData) {
                outData[i] = noData;
                continue;
            }

            int classId = static_cast<int>(std::round(lcData[i]));
            auto sensIt = sensTable.find(classId);
            double habitat = (sensIt != sensTable.end()) ? sensIt.value().habitat : 0.0;

            double D = degradation[i];
            double Dz = std::pow(D, z);
            double quality = habitat * (1.0 - Dz / (Dz + kz));

            outData[i] = std::max(0.0, std::min(1.0, quality));
            sumQ += outData[i];
            ++validCount;

            if (i % 500000 == 0)
                reportProgress(0.75 + 0.15 * static_cast<double>(i) / total);
        }

        reportProgress(0.92, "Writing output...");

        // 7. Write output
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        double meanQ = (validCount > 0) ? sumQ / validCount : 0.0;
        reportProgress(1.0,
            QString("Done. Mean habitat quality: %1, valid pixels: %2")
                .arg(meanQ, 0, 'f', 4)
                .arg(validCount));

        return true;
    }
};

REGISTER_MODULE(HabitatModule)

} // namespace aplaceholder
