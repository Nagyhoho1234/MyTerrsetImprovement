#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class TCorModule : public Module {
public:
    QString name() const override { return "TCOR"; }
    QString description() const override {
        return "Temporal correlation. Computes per-pixel Pearson correlation coefficient "
               "and p-value between two time series of raster images.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("series1", "Time series 1 rasters (comma-separated)",
                "Comma-separated list of raster file paths for the first time series"),
            ParameterDef::file("series2", "Time series 2 rasters (comma-separated)",
                "Comma-separated list of raster file paths for the second time series"),
            ParameterDef::output("output_r", "Output correlation (r) raster"),
            ParameterDef::output("output_p", "Output p-value raster"),
        };
    }

    bool execute() override {
        QString s1Param = parameter("series1").toString();
        QString s2Param = parameter("series2").toString();
        QString outRPath = parameter("output_r").toString();
        QString outPPath = parameter("output_p").toString();

        QStringList s1Paths = s1Param.split(",", Qt::SkipEmptyParts);
        QStringList s2Paths = s2Param.split(",", Qt::SkipEmptyParts);
        for (auto& p : s1Paths) p = p.trimmed();
        for (auto& p : s2Paths) p = p.trimmed();

        if (s1Paths.size() != s2Paths.size()) {
            setError("Both time series must have the same number of time steps.");
            return false;
        }

        int numSteps = s1Paths.size();
        if (numSteps < 3) {
            setError("At least 3 time steps are required for correlation analysis.");
            return false;
        }

        // Read all rasters
        reportProgress(0.0, "Reading time series...");
        std::vector<std::unique_ptr<Raster>> rst1(numSteps), rst2(numSteps);
        int cols = 0, rows = 0;

        for (int t = 0; t < numSteps; ++t) {
            rst1[t] = GdalIO::read(s1Paths[t]);
            rst2[t] = GdalIO::read(s2Paths[t]);
            if (!rst1[t] || !rst2[t]) {
                setError("Failed to read time step " + QString::number(t + 1));
                return false;
            }
            if (t == 0) {
                cols = rst1[0]->cols();
                rows = rst1[0]->rows();
            } else {
                if (rst1[t]->cols() != cols || rst1[t]->rows() != rows ||
                    rst2[t]->cols() != cols || rst2[t]->rows() != rows) {
                    setError("All rasters must have the same dimensions.");
                    return false;
                }
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = rst1[0]->hasNoData();
        double noData = rst1[0]->noDataValue();
        double outNoData = -9999.0;

        std::vector<const std::vector<double>*> d1(numSteps), d2(numSteps);
        for (int t = 0; t < numSteps; ++t) {
            d1[t] = &rst1[t]->data(0);
            d2[t] = &rst2[t]->data(0);
        }

        Raster outR(cols, rows, 1, DataType::Float32);
        Raster outP(cols, rows, 1, DataType::Float32);
        outR.setGeoTransform(rst1[0]->geoTransform());
        outR.setProjection(rst1[0]->projection());
        outR.setNoDataValue(outNoData);
        outP.setGeoTransform(rst1[0]->geoTransform());
        outP.setProjection(rst1[0]->projection());
        outP.setNoDataValue(outNoData);

        auto& rData = outR.data(0);
        auto& pData = outP.data(0);

        reportProgress(0.2, "Computing correlations...");

        for (int64_t i = 0; i < total; ++i) {
            double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
            int n = 0;

            for (int t = 0; t < numSteps; ++t) {
                double x = (*d1[t])[i];
                double y = (*d2[t])[i];
                if (hasND && (x == noData || y == noData)) continue;
                sumX += x;
                sumY += y;
                sumXY += x * y;
                sumX2 += x * x;
                sumY2 += y * y;
                ++n;
            }

            if (n < 3) {
                rData[i] = outNoData;
                pData[i] = outNoData;
                continue;
            }

            double denom = std::sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));
            if (denom < 1e-15) {
                rData[i] = 0.0;
                pData[i] = 1.0;
                continue;
            }

            double r = (n * sumXY - sumX * sumY) / denom;
            // Clamp r to [-1, 1]
            r = std::max(-1.0, std::min(1.0, r));
            rData[i] = r;

            // Approximate p-value using t-distribution
            // t = r * sqrt(n-2) / sqrt(1 - r*r)
            double r2 = r * r;
            if (r2 >= 1.0) {
                pData[i] = 0.0;
            } else {
                double tStat = std::abs(r) * std::sqrt(static_cast<double>(n - 2) / (1.0 - r2));
                int df = n - 2;
                // Approximate two-tailed p-value using the beta incomplete function approximation
                // For large df: p ~ 2 * exp(-0.5 * t^2) (normal approximation)
                // Better approximation for moderate df:
                double x = df / (df + tStat * tStat);
                // Simple approximation: p = 2 * (1 - Phi(|t|)) for df > 30
                // For smaller df, use a cruder approximation
                if (df > 30) {
                    // Normal approximation
                    double z = tStat;
                    double p = std::erfc(z / std::sqrt(2.0));
                    pData[i] = p;
                } else {
                    // Rough approximation via regularized incomplete beta
                    double p = std::pow(x, df * 0.5) * df / (df + tStat * tStat);
                    pData[i] = std::min(1.0, std::max(0.0, p));
                }
            }

            if (i % 1000000 == 0)
                reportProgress(0.2 + 0.7 * static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");
        if (!GdalIO::write(outR, outRPath)) {
            setError("Failed to write correlation output.");
            return false;
        }
        if (!GdalIO::write(outP, outPPath)) {
            setError("Failed to write p-value output.");
            return false;
        }

        reportProgress(1.0, "Temporal correlation complete.");
        return true;
    }
};

REGISTER_MODULE(TCorModule)

} // namespace aplaceholder
