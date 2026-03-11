#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>

namespace aplaceholder {

class ValidateModule : public Module {
public:
    QString name() const override { return "VALIDATE"; }
    QString description() const override {
        return "Validation statistics comparing predicted vs observed rasters. "
               "Computes RMSE, MAE, R-squared, bias, and Nash-Sutcliffe efficiency.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("predicted", "Predicted values raster"),
            ParameterDef::file("observed", "Observed values raster"),
            ParameterDef::output("output_report", "Output validation report text file"),
        };
    }

    bool execute() override {
        auto rP = GdalIO::read(parameter("predicted").toString());
        auto rO = GdalIO::read(parameter("observed").toString());
        if (!rP || !rO) {
            setError("Failed to read input rasters");
            return false;
        }

        if (rP->cols() != rO->cols() || rP->rows() != rO->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = rP->cols(), rows = rP->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& dP = rP->data(0);
        const auto& dO = rO->data(0);
        double noDataP = rP->noDataValue();
        double noDataO = rO->noDataValue();
        bool hasNDP = rP->hasNoData();
        bool hasNDO = rO->hasNoData();

        reportProgress(0.0, "Computing validation statistics...");

        // First pass: accumulate sums
        double sumP = 0.0, sumO = 0.0;
        double sumSqErr = 0.0, sumAbsErr = 0.0;
        int64_t N = 0;

        for (int64_t i = 0; i < total; ++i) {
            if (hasNDP && dP[i] == noDataP) continue;
            if (hasNDO && dO[i] == noDataO) continue;
            double p = dP[i], o = dO[i];
            double err = p - o;
            sumP += p;
            sumO += o;
            sumSqErr += err * err;
            sumAbsErr += std::abs(err);
            ++N;
        }

        if (N < 1) {
            setError("No valid pixel pairs found");
            return false;
        }

        double Nd = static_cast<double>(N);
        double meanP = sumP / Nd;
        double meanO = sumO / Nd;
        double RMSE = std::sqrt(sumSqErr / Nd);
        double MAE = sumAbsErr / Nd;
        double bias = meanP - meanO;

        reportProgress(0.5, "Computing R² and Nash-Sutcliffe...");

        // Second pass for R² and Nash-Sutcliffe
        double sumSqDevO = 0.0;
        double sumSqDevP = 0.0;
        double sumCrossDev = 0.0;

        for (int64_t i = 0; i < total; ++i) {
            if (hasNDP && dP[i] == noDataP) continue;
            if (hasNDO && dO[i] == noDataO) continue;
            double devO = dO[i] - meanO;
            double devP = dP[i] - meanP;
            sumSqDevO += devO * devO;
            sumSqDevP += devP * devP;
            sumCrossDev += devO * devP;
        }

        // Pearson R²
        double rSquared = 0.0;
        if (sumSqDevO > 0.0 && sumSqDevP > 0.0) {
            double r = sumCrossDev / std::sqrt(sumSqDevO * sumSqDevP);
            rSquared = r * r;
        }

        // Nash-Sutcliffe Efficiency: 1 - (sum(O-P)^2 / sum(O-meanO)^2)
        double NSE = (sumSqDevO > 0.0) ? 1.0 - sumSqErr / sumSqDevO : -999.0;

        reportProgress(0.8, "Writing report...");

        // Write report
        QString reportPath = parameter("output_report").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        outFile << "Validation Statistics Report\n";
        outFile << "============================\n\n";
        outFile << "Predicted: " << parameter("predicted").toString().toStdString() << "\n";
        outFile << "Observed: " << parameter("observed").toString().toStdString() << "\n\n";
        outFile << "Number of valid pixel pairs (N): " << N << "\n\n";
        outFile << "Mean predicted: " << meanP << "\n";
        outFile << "Mean observed: " << meanO << "\n\n";
        outFile << "RMSE (Root Mean Squared Error): " << RMSE << "\n";
        outFile << "MAE (Mean Absolute Error): " << MAE << "\n";
        outFile << "R-squared: " << rSquared << "\n";
        outFile << "Bias (mean predicted - mean observed): " << bias << "\n";
        outFile << "Nash-Sutcliffe Efficiency (NSE): " << NSE << "\n\n";

        outFile << "Interpretation:\n";
        outFile << "  NSE = 1: Perfect prediction\n";
        outFile << "  NSE = 0: Predictions are as accurate as the mean of observations\n";
        outFile << "  NSE < 0: Mean of observations is a better predictor\n";
        outFile.close();

        reportProgress(1.0,
            QString("RMSE=%1, MAE=%2, R²=%3, Bias=%4, NSE=%5")
                .arg(RMSE, 0, 'f', 6)
                .arg(MAE, 0, 'f', 6)
                .arg(rSquared, 0, 'f', 6)
                .arg(bias, 0, 'f', 6)
                .arg(NSE, 0, 'f', 6));

        return true;
    }
};

REGISTER_MODULE(ValidateModule)

} // namespace aplaceholder
