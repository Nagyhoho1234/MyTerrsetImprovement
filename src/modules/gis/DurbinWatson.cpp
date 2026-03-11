#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>

namespace aplaceholder {

class DurbinWatsonModule : public Module {
public:
    QString name() const override { return "DURBINWATSON"; }
    QString description() const override {
        return "Durbin-Watson test for autocorrelation in regression residuals. "
               "Tests whether successive residual values (in row-major scan order) "
               "are correlated, indicating spatial or serial autocorrelation.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("residual_raster", "Residual raster from regression"),
            ParameterDef::output("output_report", "Output Durbin-Watson report text file"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("residual_raster").toString());
        if (!raster) {
            setError("Failed to read residual raster");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& data = raster->data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        reportProgress(0.0, "Computing Durbin-Watson statistic...");

        // Collect valid residuals in scan order
        std::vector<double> residuals;
        residuals.reserve(total);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            residuals.push_back(data[i]);
        }

        int64_t N = static_cast<int64_t>(residuals.size());
        if (N < 3) {
            setError("Insufficient valid residual pixels (need at least 3)");
            return false;
        }

        reportProgress(0.3, "Computing DW statistic...");

        // DW = sum((e_t - e_{t-1})^2) / sum(e_t^2)
        double sumSqDiff = 0.0;
        double sumSq = residuals[0] * residuals[0];
        double sumResiduals = residuals[0];
        double sumSqResiduals = residuals[0] * residuals[0];

        for (int64_t t = 1; t < N; ++t) {
            double diff = residuals[t] - residuals[t - 1];
            sumSqDiff += diff * diff;
            sumSq += residuals[t] * residuals[t];
            sumResiduals += residuals[t];
            sumSqResiduals += residuals[t] * residuals[t];

            if (t % 1000000 == 0)
                reportProgress(0.3 + 0.5 * static_cast<double>(t) / N);
        }

        double DW = (sumSq > 0.0) ? sumSqDiff / sumSq : 0.0;

        // Estimated first-order autocorrelation: rho ≈ 1 - DW/2
        double rho = 1.0 - DW / 2.0;

        // Mean of residuals
        double meanResidual = sumResiduals / N;

        reportProgress(0.9, "Writing report...");

        // Write report
        QString reportPath = parameter("output_report").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        outFile << "Durbin-Watson Autocorrelation Test\n";
        outFile << "===================================\n\n";
        outFile << "Residual raster: " << parameter("residual_raster").toString().toStdString() << "\n\n";
        outFile << "Number of valid residuals (N): " << N << "\n";
        outFile << "Mean residual: " << meanResidual << "\n";
        outFile << "Sum of squared residuals: " << sumSq << "\n\n";
        outFile << "Durbin-Watson statistic (DW): " << DW << "\n";
        outFile << "Estimated autocorrelation (rho): " << rho << "\n\n";

        outFile << "Interpretation:\n";
        outFile << "  DW ≈ 2.0: No autocorrelation\n";
        outFile << "  DW < 2.0: Positive autocorrelation (DW → 0)\n";
        outFile << "  DW > 2.0: Negative autocorrelation (DW → 4)\n\n";

        if (DW < 1.5)
            outFile << "Result: Evidence of positive autocorrelation in residuals\n";
        else if (DW > 2.5)
            outFile << "Result: Evidence of negative autocorrelation in residuals\n";
        else
            outFile << "Result: No strong evidence of autocorrelation\n";

        outFile.close();

        reportProgress(1.0,
            QString("DW = %1, rho = %2")
                .arg(DW, 0, 'f', 6)
                .arg(rho, 0, 'f', 6));

        return true;
    }
};

REGISTER_MODULE(DurbinWatsonModule)

} // namespace aplaceholder
