#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>

namespace aplaceholder {

class RegressModule : public Module {
public:
    QString name() const override { return "REGRESS"; }
    QString description() const override {
        return "Simple linear regression between two raster images. "
               "Computes Y = a + bX with slope, intercept, R-squared, "
               "standard error, and p-value. Outputs predicted and residual rasters.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dependent", "Dependent variable raster (Y)"),
            ParameterDef::file("independent", "Independent variable raster (X)"),
            ParameterDef::output("output_report", "Output report text file"),
            ParameterDef::output("output_predicted", "Output predicted values raster"),
            ParameterDef::output("output_residual", "Output residual raster"),
        };
    }

    bool execute() override {
        auto rY = GdalIO::read(parameter("dependent").toString());
        auto rX = GdalIO::read(parameter("independent").toString());
        if (!rY || !rX) {
            setError("Failed to read input rasters");
            return false;
        }

        if (rY->cols() != rX->cols() || rY->rows() != rX->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = rY->cols(), rows = rY->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& dY = rY->data(0);
        const auto& dX = rX->data(0);
        double noDataY = rY->noDataValue();
        double noDataX = rX->noDataValue();
        bool hasNDY = rY->hasNoData();
        bool hasNDX = rX->hasNoData();

        reportProgress(0.0, "Computing regression statistics...");

        double sumX = 0.0, sumY = 0.0, sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
        int64_t N = 0;

        for (int64_t i = 0; i < total; ++i) {
            if (hasNDY && dY[i] == noDataY) continue;
            if (hasNDX && dX[i] == noDataX) continue;
            double x = dX[i], y = dY[i];
            sumX += x;
            sumY += y;
            sumXX += x * x;
            sumYY += y * y;
            sumXY += x * y;
            ++N;
        }

        if (N < 3) {
            setError("Insufficient valid pixel pairs for regression (need at least 3)");
            return false;
        }

        double Nd = static_cast<double>(N);
        double meanX = sumX / Nd;
        double meanY = sumY / Nd;
        double Sxx = sumXX - Nd * meanX * meanX;
        double Syy = sumYY - Nd * meanY * meanY;
        double Sxy = sumXY - Nd * meanX * meanY;

        if (std::abs(Sxx) < 1e-15) {
            setError("Independent variable has zero variance");
            return false;
        }

        double slope = Sxy / Sxx;
        double intercept = meanY - slope * meanX;

        // Compute R-squared and standard error
        double SSres = 0.0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasNDY && dY[i] == noDataY) continue;
            if (hasNDX && dX[i] == noDataX) continue;
            double predicted = intercept + slope * dX[i];
            double residual = dY[i] - predicted;
            SSres += residual * residual;
        }

        double SStot = Syy;
        double rSquared = (SStot > 0.0) ? 1.0 - SSres / SStot : 0.0;
        double stdError = std::sqrt(SSres / (Nd - 2.0));
        double seSlope = stdError / std::sqrt(Sxx);

        // t-statistic for slope
        double tStat = (seSlope > 0.0) ? slope / seSlope : 0.0;

        // p-value approximation
        double pValue = 1.0;
        if (std::abs(tStat) > 0.0 && N > 2) {
            double absT = std::abs(tStat);
            double df = Nd - 2.0;
            if (df > 30) {
                pValue = std::erfc(absT / std::sqrt(2.0));
            } else {
                double x = df / (df + tStat * tStat);
                double a2 = df / 2.0;
                pValue = std::pow(x, a2) / (a2 * std::exp(std::lgamma(a2) + std::lgamma(0.5) - std::lgamma(a2 + 0.5)));
                if (pValue > 1.0) pValue = 1.0;
                if (pValue < 0.0) pValue = 0.0;
            }
        }

        reportProgress(0.4, "Writing report...");

        // Write report
        QString reportPath = parameter("output_report").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        outFile << "Simple Linear Regression Analysis\n";
        outFile << "==================================\n\n";
        outFile << "Dependent (Y): " << parameter("dependent").toString().toStdString() << "\n";
        outFile << "Independent (X): " << parameter("independent").toString().toStdString() << "\n\n";
        outFile << "Number of valid pixel pairs (N): " << N << "\n\n";
        outFile << "Regression equation: Y = " << intercept << " + " << slope << " * X\n\n";
        outFile << "Slope (b): " << slope << "\n";
        outFile << "Intercept (a): " << intercept << "\n";
        outFile << "R-squared: " << rSquared << "\n";
        outFile << "Standard error of estimate: " << stdError << "\n";
        outFile << "Standard error of slope: " << seSlope << "\n";
        outFile << "t-statistic (slope): " << tStat << "\n";
        outFile << "p-value (approx): " << pValue << "\n";
        outFile.close();

        reportProgress(0.6, "Writing predicted raster...");

        // Output predicted raster
        double outNoData = -9999.0;

        Raster predicted(cols, rows, 1, DataType::Float64);
        predicted.setGeoTransform(rY->geoTransform());
        predicted.setProjection(rY->projection());
        predicted.setNoDataValue(outNoData);
        auto& predData = predicted.data(0);

        Raster residual(cols, rows, 1, DataType::Float64);
        residual.setGeoTransform(rY->geoTransform());
        residual.setProjection(rY->projection());
        residual.setNoDataValue(outNoData);
        auto& resData = residual.data(0);

        for (int64_t i = 0; i < total; ++i) {
            if ((hasNDY && dY[i] == noDataY) || (hasNDX && dX[i] == noDataX)) {
                predData[i] = outNoData;
                resData[i] = outNoData;
            } else {
                double pred = intercept + slope * dX[i];
                predData[i] = pred;
                resData[i] = dY[i] - pred;
            }
        }

        if (!GdalIO::write(predicted, parameter("output_predicted").toString())) {
            setError("Failed to write predicted raster");
            return false;
        }

        reportProgress(0.8, "Writing residual raster...");

        if (!GdalIO::write(residual, parameter("output_residual").toString())) {
            setError("Failed to write residual raster");
            return false;
        }

        reportProgress(1.0,
            QString("Y = %1 + %2*X, R² = %3, SE = %4")
                .arg(intercept, 0, 'f', 6)
                .arg(slope, 0, 'f', 6)
                .arg(rSquared, 0, 'f', 6)
                .arg(stdError, 0, 'f', 6));

        return true;
    }
};

REGISTER_MODULE(RegressModule)

} // namespace aplaceholder
