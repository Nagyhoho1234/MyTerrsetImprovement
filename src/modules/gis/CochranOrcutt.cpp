#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <vector>

namespace aplaceholder {

class CochranOrcuttModule : public Module {
public:
    QString name() const override { return "COCHRANORCUTT"; }
    QString description() const override {
        return "Cochrane-Orcutt iterative procedure for correcting autocorrelated "
               "regression residuals. Estimates the autocorrelation coefficient and "
               "transforms the data to produce corrected regression coefficients.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dependent", "Dependent variable raster (Y)"),
            ParameterDef::file("independent", "Independent variable raster (X)"),
            ParameterDef::output("output_coefficients", "Output coefficients report text file"),
            ParameterDef::output("output_residual", "Output corrected residual raster"),
            ParameterDef::integer("max_iterations", "Maximum iterations",
                10, 1, 100, "Maximum number of Cochrane-Orcutt iterations"),
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
        int maxIter = parameter("max_iterations").toInt();
        if (maxIter <= 0) maxIter = 10;

        reportProgress(0.0, "Collecting valid pixel pairs in scan order...");

        // Collect valid pairs in scan order for serial correlation
        std::vector<double> x, y;
        for (int64_t i = 0; i < total; ++i) {
            if (hasNDY && dY[i] == noDataY) continue;
            if (hasNDX && dX[i] == noDataX) continue;
            x.push_back(dX[i]);
            y.push_back(dY[i]);
        }

        int64_t N = static_cast<int64_t>(x.size());
        if (N < 4) {
            setError("Insufficient valid pixel pairs (need at least 4)");
            return false;
        }

        // Helper: OLS regression on vectors
        auto ols = [](const std::vector<double>& xv, const std::vector<double>& yv,
                      double& slope, double& intercept) {
            int64_t n = static_cast<int64_t>(xv.size());
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (int64_t i = 0; i < n; ++i) {
                sx += xv[i]; sy += yv[i];
                sxx += xv[i] * xv[i]; sxy += xv[i] * yv[i];
            }
            double nd = static_cast<double>(n);
            double denom = sxx - sx * sx / nd;
            if (std::abs(denom) < 1e-15) {
                slope = 0.0;
                intercept = sy / nd;
            } else {
                slope = (sxy - sx * sy / nd) / denom;
                intercept = (sy - slope * sx) / nd;
            }
        };

        reportProgress(0.1, "Running Cochrane-Orcutt iterations...");

        // Initial OLS
        double slope, intercept;
        ols(x, y, slope, intercept);

        // Compute residuals
        std::vector<double> residuals(N);
        for (int64_t i = 0; i < N; ++i) {
            residuals[i] = y[i] - intercept - slope * x[i];
        }

        double rho = 0.0;
        double prevRho = 0.0;
        double convergenceTol = 1e-6;

        struct IterResult { int iter; double rho; double slope; double intercept; };
        std::vector<IterResult> iterResults;

        for (int iter = 0; iter < maxIter; ++iter) {
            // Step 1: Estimate rho from residuals
            double sumNum = 0.0, sumDen = 0.0;
            for (int64_t t = 1; t < N; ++t) {
                sumNum += residuals[t] * residuals[t - 1];
                sumDen += residuals[t - 1] * residuals[t - 1];
            }
            rho = (sumDen > 0.0) ? sumNum / sumDen : 0.0;

            // Step 2: Transform data using rho
            std::vector<double> yTrans(N - 1), xTrans(N - 1);
            for (int64_t t = 1; t < N; ++t) {
                yTrans[t - 1] = y[t] - rho * y[t - 1];
                xTrans[t - 1] = x[t] - rho * x[t - 1];
            }

            // Step 3: OLS on transformed data
            ols(xTrans, yTrans, slope, intercept);

            // The original intercept: a_original = intercept / (1 - rho)
            double originalIntercept = (std::abs(1.0 - rho) > 1e-15) ?
                intercept / (1.0 - rho) : intercept;

            iterResults.push_back({iter + 1, rho, slope, originalIntercept});

            // Recompute residuals with original variables
            for (int64_t i = 0; i < N; ++i) {
                residuals[i] = y[i] - originalIntercept - slope * x[i];
            }

            reportProgress(0.1 + 0.7 * static_cast<double>(iter + 1) / maxIter,
                QString("Iteration %1: rho = %2").arg(iter + 1).arg(rho, 0, 'f', 6));

            // Check convergence
            if (iter > 0 && std::abs(rho - prevRho) < convergenceTol) {
                break;
            }
            prevRho = rho;
        }

        reportProgress(0.85, "Writing report and residual raster...");

        // Final coefficients
        double finalSlope = iterResults.back().slope;
        double finalIntercept = iterResults.back().intercept;
        double finalRho = iterResults.back().rho;

        // Write report
        QString reportPath = parameter("output_coefficients").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        outFile << "Cochrane-Orcutt Iterative Estimation\n";
        outFile << "=====================================\n\n";
        outFile << "Dependent (Y): " << parameter("dependent").toString().toStdString() << "\n";
        outFile << "Independent (X): " << parameter("independent").toString().toStdString() << "\n\n";
        outFile << "Number of valid pixel pairs (N): " << N << "\n";
        outFile << "Max iterations: " << maxIter << "\n";
        outFile << "Convergence tolerance: " << convergenceTol << "\n\n";

        outFile << "Iteration History:\n";
        outFile << "Iter\tRho\tSlope\tIntercept\n";
        for (const auto& ir : iterResults) {
            outFile << ir.iter << "\t" << ir.rho << "\t" << ir.slope << "\t" << ir.intercept << "\n";
        }

        outFile << "\nFinal Corrected Regression:\n";
        outFile << "  Y = " << finalIntercept << " + " << finalSlope << " * X\n";
        outFile << "  Autocorrelation coefficient (rho): " << finalRho << "\n";
        outFile << "  Converged after " << iterResults.size() << " iterations\n";
        outFile.close();

        // Write corrected residual raster
        double outNoData = -9999.0;
        Raster outResidual(cols, rows, 1, DataType::Float64);
        outResidual.setGeoTransform(rY->geoTransform());
        outResidual.setProjection(rY->projection());
        outResidual.setNoDataValue(outNoData);
        auto& resData = outResidual.data(0);

        for (int64_t i = 0; i < total; ++i) {
            if ((hasNDY && dY[i] == noDataY) || (hasNDX && dX[i] == noDataX)) {
                resData[i] = outNoData;
            } else {
                resData[i] = dY[i] - finalIntercept - finalSlope * dX[i];
            }
        }

        if (!GdalIO::write(outResidual, parameter("output_residual").toString())) {
            setError("Failed to write corrected residual raster");
            return false;
        }

        reportProgress(1.0,
            QString("Cochrane-Orcutt: rho=%1, Y=%2+%3*X, %4 iterations")
                .arg(finalRho, 0, 'f', 6)
                .arg(finalIntercept, 0, 'f', 6)
                .arg(finalSlope, 0, 'f', 6)
                .arg(iterResults.size()));

        return true;
    }
};

REGISTER_MODULE(CochranOrcuttModule)

} // namespace aplaceholder
