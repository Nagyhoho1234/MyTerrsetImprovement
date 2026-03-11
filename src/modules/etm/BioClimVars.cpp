#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class BioClimVarsModule : public Module {
public:
    QString name() const override { return "BIOCLIMVARS"; }
    QString description() const override {
        return "Compute bioclimatic variables (BIO1-BIO19) from monthly temperature "
               "and precipitation data. These variables represent annual trends, "
               "seasonality, and extreme conditions commonly used in species "
               "distribution modeling.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("monthly_temp", "Monthly temperature rasters (comma-sep, 12)",
                "Comma-separated list of 12 monthly mean temperature raster file paths (Jan-Dec)"),
            ParameterDef::file("monthly_precip", "Monthly precipitation rasters (comma-sep, 12)",
                "Comma-separated list of 12 monthly precipitation raster file paths (Jan-Dec)"),
            ParameterDef::output("output_prefix", "Output filename prefix",
                "Output rasters will be named prefix_bio01 through prefix_bio19"),
        };
    }

    bool execute() override {
        QString tempParam = parameter("monthly_temp").toString();
        QString precipParam = parameter("monthly_precip").toString();
        QString prefix = parameter("output_prefix").toString();

        QStringList tempPaths = tempParam.split(",", Qt::SkipEmptyParts);
        QStringList precipPaths = precipParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : tempPaths) p = p.trimmed();
        for (auto& p : precipPaths) p = p.trimmed();

        if (tempPaths.size() != 12 || precipPaths.size() != 12) {
            setError("Exactly 12 monthly temperature and 12 monthly precipitation rasters required.");
            return false;
        }

        // Read all monthly rasters
        reportProgress(0.0, "Reading monthly rasters...");
        std::vector<std::unique_ptr<Raster>> tempRst(12), precipRst(12);
        int cols = 0, rows = 0;

        for (int m = 0; m < 12; ++m) {
            tempRst[m] = GdalIO::read(tempPaths[m]);
            precipRst[m] = GdalIO::read(precipPaths[m]);
            if (!tempRst[m] || !precipRst[m]) {
                setError("Failed to read month " + QString::number(m + 1) + " raster.");
                return false;
            }
            if (m == 0) {
                cols = tempRst[0]->cols();
                rows = tempRst[0]->rows();
            } else {
                if (tempRst[m]->cols() != cols || tempRst[m]->rows() != rows ||
                    precipRst[m]->cols() != cols || precipRst[m]->rows() != rows) {
                    setError("All monthly rasters must have the same dimensions.");
                    return false;
                }
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = tempRst[0]->hasNoData();
        double noData = tempRst[0]->noDataValue();
        double outNoData = -9999.0;

        // Collect data pointers
        std::vector<const std::vector<double>*> tData(12), pData(12);
        for (int m = 0; m < 12; ++m) {
            tData[m] = &tempRst[m]->data(0);
            pData[m] = &precipRst[m]->data(0);
        }

        // Create 19 output rasters
        static const char* bioNames[] = {
            "bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07",
            "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14",
            "bio15", "bio16", "bio17", "bio18", "bio19"
        };

        std::vector<Raster> outputs;
        outputs.reserve(19);
        for (int v = 0; v < 19; ++v) {
            outputs.emplace_back(cols, rows, 1, DataType::Float32);
            outputs[v].setGeoTransform(tempRst[0]->geoTransform());
            outputs[v].setProjection(tempRst[0]->projection());
            outputs[v].setNoDataValue(outNoData);
        }

        reportProgress(0.1, "Computing bioclimatic variables...");

        for (int64_t i = 0; i < total; ++i) {
            // Check nodata
            bool valid = true;
            if (hasND) {
                for (int m = 0; m < 12 && valid; ++m) {
                    if ((*tData[m])[i] == noData || (*pData[m])[i] == noData)
                        valid = false;
                }
            }

            if (!valid) {
                for (int v = 0; v < 19; ++v)
                    outputs[v].data(0)[i] = outNoData;
                continue;
            }

            double t[12], p[12];
            for (int m = 0; m < 12; ++m) {
                t[m] = (*tData[m])[i];
                p[m] = (*pData[m])[i];
            }

            // Basic monthly stats
            double tMin = t[0], tMax = t[0], tSum = 0.0;
            double pMin = p[0], pMax = p[0], pSum = 0.0;
            int tMinMonth = 0, tMaxMonth = 0, pMinMonth = 0, pMaxMonth = 0;

            for (int m = 0; m < 12; ++m) {
                tSum += t[m];
                pSum += p[m];
                if (t[m] < tMin) { tMin = t[m]; tMinMonth = m; }
                if (t[m] > tMax) { tMax = t[m]; tMaxMonth = m; }
                if (p[m] < pMin) { pMin = p[m]; pMinMonth = m; }
                if (p[m] > pMax) { pMax = p[m]; pMaxMonth = m; }
            }

            double tMean = tSum / 12.0;
            double pMean = pSum / 12.0;

            // Temperature range per month (diurnal range proxy not available, use monthly range)
            // BIO2: Mean Diurnal Range - approximate as mean of monthly ranges
            // Without daily data, use the difference between consecutive months as proxy
            double meanRange = tMax - tMin;  // simplified

            // BIO4: Temperature seasonality (std dev * 100)
            double tSumSq = 0.0;
            for (int m = 0; m < 12; ++m)
                tSumSq += (t[m] - tMean) * (t[m] - tMean);
            double tStdDev = std::sqrt(tSumSq / 12.0);

            // Quarter sums (rolling 3-month)
            double qTemp[12], qPrecip[12];
            for (int m = 0; m < 12; ++m) {
                qTemp[m] = t[m] + t[(m + 1) % 12] + t[(m + 2) % 12];
                qPrecip[m] = p[m] + p[(m + 1) % 12] + p[(m + 2) % 12];
            }

            double qTempMin = qTemp[0], qTempMax = qTemp[0];
            double qPrecipMin = qPrecip[0], qPrecipMax = qPrecip[0];
            int qTempMinIdx = 0, qTempMaxIdx = 0;
            int qPrecipMinIdx = 0, qPrecipMaxIdx = 0;

            for (int m = 1; m < 12; ++m) {
                if (qTemp[m] < qTempMin) { qTempMin = qTemp[m]; qTempMinIdx = m; }
                if (qTemp[m] > qTempMax) { qTempMax = qTemp[m]; qTempMaxIdx = m; }
                if (qPrecip[m] < qPrecipMin) { qPrecipMin = qPrecip[m]; qPrecipMinIdx = m; }
                if (qPrecip[m] > qPrecipMax) { qPrecipMax = qPrecip[m]; qPrecipMaxIdx = m; }
            }

            // BIO15: Precipitation seasonality (CV)
            double pSumSq = 0.0;
            for (int m = 0; m < 12; ++m)
                pSumSq += (p[m] - pMean) * (p[m] - pMean);
            double pStdDev = std::sqrt(pSumSq / 12.0);
            double precipCV = (pMean > 0.0) ? (pStdDev / pMean) * 100.0 : 0.0;

            // Assign BIO variables
            outputs[0].data(0)[i] = tMean;                          // BIO1: Annual Mean Temp
            outputs[1].data(0)[i] = meanRange;                      // BIO2: Mean Diurnal Range
            outputs[2].data(0)[i] = (tMax - tMin > 0) ?             // BIO3: Isothermality
                (meanRange / (tMax - tMin)) * 100.0 : 0.0;
            outputs[3].data(0)[i] = tStdDev * 100.0;                // BIO4: Temp Seasonality
            outputs[4].data(0)[i] = tMax;                            // BIO5: Max Temp Warmest Month
            outputs[5].data(0)[i] = tMin;                            // BIO6: Min Temp Coldest Month
            outputs[6].data(0)[i] = tMax - tMin;                     // BIO7: Temp Annual Range
            outputs[7].data(0)[i] = qTemp[qPrecipMaxIdx] / 3.0;     // BIO8: Mean Temp Wettest Quarter
            outputs[8].data(0)[i] = qTemp[qPrecipMinIdx] / 3.0;     // BIO9: Mean Temp Driest Quarter
            outputs[9].data(0)[i] = qTempMax / 3.0;                  // BIO10: Mean Temp Warmest Quarter
            outputs[10].data(0)[i] = qTempMin / 3.0;                 // BIO11: Mean Temp Coldest Quarter
            outputs[11].data(0)[i] = pSum;                            // BIO12: Annual Precipitation
            outputs[12].data(0)[i] = pMax;                            // BIO13: Precip Wettest Month
            outputs[13].data(0)[i] = pMin;                            // BIO14: Precip Driest Month
            outputs[14].data(0)[i] = precipCV;                        // BIO15: Precip Seasonality
            outputs[15].data(0)[i] = qPrecipMax;                     // BIO16: Precip Wettest Quarter
            outputs[16].data(0)[i] = qPrecipMin;                     // BIO17: Precip Driest Quarter
            outputs[17].data(0)[i] = qPrecip[qTempMaxIdx];           // BIO18: Precip Warmest Quarter
            outputs[18].data(0)[i] = qPrecip[qTempMinIdx];           // BIO19: Precip Coldest Quarter

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.7 * static_cast<double>(i) / total);
        }

        // Write outputs
        reportProgress(0.8, "Writing output rasters...");
        for (int v = 0; v < 19; ++v) {
            QString outPath = prefix + "_" + bioNames[v];
            if (!GdalIO::write(outputs[v], outPath)) {
                setError("Failed to write: " + outPath);
                return false;
            }
        }

        reportProgress(1.0, "Bioclimatic variables computation complete.");
        return true;
    }
};

REGISTER_MODULE(BioClimVarsModule)

} // namespace aplaceholder
