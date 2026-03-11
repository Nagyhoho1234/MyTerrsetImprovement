#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>

namespace aplaceholder {

class CorrelateModule : public Module {
public:
    QString name() const override { return "CORRELATE"; }
    QString description() const override {
        return "Pearson correlation between two raster images. "
               "Computes r, R-squared, and p-value. Optionally produces "
               "a residual raster from the linear fit.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First input raster image"),
            ParameterDef::file("input2", "Second input raster image"),
            ParameterDef::output("output_report", "Output report text file"),
            ParameterDef::output("output_residual", "Output residual raster (optional)",
                "Residual raster from linear fit of input2 on input1"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input1").toString());
        auto r2 = GdalIO::read(parameter("input2").toString());
        if (!r1 || !r2) {
            setError("Failed to read input rasters");
            return false;
        }

        if (r1->cols() != r2->cols() || r1->rows() != r2->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& d1 = r1->data(0);
        const auto& d2 = r2->data(0);
        double noData1 = r1->noDataValue();
        double noData2 = r2->noDataValue();
        bool hasND1 = r1->hasNoData();
        bool hasND2 = r2->hasNoData();

        reportProgress(0.0, "Computing correlation statistics...");

        // Accumulate sums for Pearson correlation
        double sumX = 0.0, sumY = 0.0, sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
        int64_t N = 0;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND1 && d1[i] == noData1) continue;
            if (hasND2 && d2[i] == noData2) continue;
            double x = d1[i], y = d2[i];
            sumX += x;
            sumY += y;
            sumXX += x * x;
            sumYY += y * y;
            sumXY += x * y;
            ++N;
        }

        if (N < 3) {
            setError("Insufficient valid pixel pairs for correlation (need at least 3)");
            return false;
        }

        double Nd = static_cast<double>(N);
        double meanX = sumX / Nd;
        double meanY = sumY / Nd;
        double varX = sumXX / Nd - meanX * meanX;
        double varY = sumYY / Nd - meanY * meanY;
        double covXY = sumXY / Nd - meanX * meanY;

        double r = 0.0;
        if (varX > 0.0 && varY > 0.0) {
            r = covXY / std::sqrt(varX * varY);
        }
        double rSquared = r * r;

        // t-test for significance of r
        double tStat = 0.0;
        double pValue = 1.0;
        if (std::abs(r) < 1.0 && N > 2) {
            tStat = r * std::sqrt((Nd - 2.0) / (1.0 - r * r));
            // Approximate two-tailed p-value using normal approximation for large N
            double absT = std::abs(tStat);
            double df = Nd - 2.0;
            // Use incomplete beta approximation: p = 1 - betai(df/2, 0.5, df/(df+t^2))
            // Simplified: for large df, t ~ N(0,1)
            // More accurate approximation using the formula:
            double x = df / (df + tStat * tStat);
            // Regularized incomplete beta function approximation
            // For large df, use normal approximation
            if (df > 30) {
                // Normal approximation
                pValue = 2.0 * std::erfc(absT / std::sqrt(2.0)) / 2.0;
                pValue = std::erfc(absT / std::sqrt(2.0));
            } else {
                // Simple approximation for smaller df
                double a = df / 2.0;
                pValue = std::pow(x, a) / (a * std::exp(std::lgamma(a) + std::lgamma(0.5) - std::lgamma(a + 0.5)));
                // Clamp
                if (pValue > 1.0) pValue = 1.0;
                if (pValue < 0.0) pValue = 0.0;
            }
        } else if (std::abs(r) >= 1.0) {
            pValue = 0.0;
        }

        // Linear regression coefficients for residual: y = a + b*x
        double b = (varX > 0.0) ? covXY / varX : 0.0;
        double a = meanY - b * meanX;

        reportProgress(0.5, "Writing report...");

        // Write report
        QString reportPath = parameter("output_report").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        outFile << "Pearson Correlation Analysis\n";
        outFile << "============================\n\n";
        outFile << "Input 1: " << parameter("input1").toString().toStdString() << "\n";
        outFile << "Input 2: " << parameter("input2").toString().toStdString() << "\n\n";
        outFile << "Number of valid pixel pairs (N): " << N << "\n\n";
        outFile << "Mean of Input 1: " << meanX << "\n";
        outFile << "Mean of Input 2: " << meanY << "\n\n";
        outFile << "Pearson r: " << r << "\n";
        outFile << "R-squared: " << rSquared << "\n";
        outFile << "t-statistic: " << tStat << "\n";
        outFile << "p-value (approx): " << pValue << "\n\n";
        outFile << "Linear fit: Y = " << a << " + " << b << " * X\n";
        outFile.close();

        // Optional residual raster
        QString residualPath = parameter("output_residual").toString();
        if (!residualPath.isEmpty()) {
            reportProgress(0.7, "Computing residual raster...");

            Raster residual(cols, rows, 1, DataType::Float64);
            residual.setGeoTransform(r1->geoTransform());
            residual.setProjection(r1->projection());
            double outNoData = -9999.0;
            residual.setNoDataValue(outNoData);
            auto& resData = residual.data(0);

            for (int64_t i = 0; i < total; ++i) {
                if ((hasND1 && d1[i] == noData1) || (hasND2 && d2[i] == noData2)) {
                    resData[i] = outNoData;
                } else {
                    double predicted = a + b * d1[i];
                    resData[i] = d2[i] - predicted;
                }
            }

            if (!GdalIO::write(residual, residualPath)) {
                setError("Failed to write residual raster");
                return false;
            }
        }

        reportProgress(1.0,
            QString("Pearson r = %1, R² = %2, p = %3")
                .arg(r, 0, 'f', 6)
                .arg(rSquared, 0, 'f', 6)
                .arg(pValue, 0, 'g', 4));

        return true;
    }
};

REGISTER_MODULE(CorrelateModule)

} // namespace aplaceholder
