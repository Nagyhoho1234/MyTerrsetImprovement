#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>

namespace aplaceholder {

class AutocorrModule : public Module {
public:
    QString name() const override { return "AUTOCORR"; }
    QString description() const override {
        return "Moran's I spatial autocorrelation analysis. "
               "Measures the degree of spatial clustering or dispersion "
               "in a raster image using queen contiguity at a given lag distance.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::real("lag_distance", "Lag distance", 1.0, 1.0, 1000.0,
                "Distance lag in cell units for neighbor definition (1.0 = immediate neighbors)"),
            ParameterDef::output("output", "Output statistics text file",
                "Text file with Moran's I results"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& data = raster->data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();
        double lagDist = parameter("lag_distance").toDouble();
        int lag = static_cast<int>(std::round(lagDist));
        if (lag < 1) lag = 1;

        reportProgress(0.0, "Computing mean...");

        // Compute mean of valid cells
        double sum = 0.0;
        int64_t N = 0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            sum += data[i];
            ++N;
        }

        if (N < 2) {
            setError("Insufficient valid cells for autocorrelation analysis");
            return false;
        }

        double mean = sum / N;

        reportProgress(0.1, "Computing Moran's I...");

        // Compute Moran's I using queen contiguity within lag distance
        // I = (N/W) * sum_i(sum_j(w_ij * (x_i - mean)(x_j - mean))) / sum_i((x_i - mean)^2)
        //
        // For queen contiguity at lag distance d, neighbors are cells within
        // Chebyshev distance <= lag (i.e., max(|dr|, |dc|) <= lag)

        double numerator = 0.0;   // sum_i(sum_j(w_ij * zi * zj))
        double denominator = 0.0; // sum_i(zi^2)
        double W = 0.0;           // total weight sum

        // Precompute deviations from mean
        std::vector<double> z(total);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData)
                z[i] = noData;
            else
                z[i] = data[i] - mean;
        }

        // Compute denominator: sum of squared deviations
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            denominator += z[i] * z[i];
        }

        if (std::abs(denominator) < 1e-15) {
            setError("All values are identical; cannot compute autocorrelation");
            return false;
        }

        // Compute numerator and W using queen contiguity within lag
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx_i = static_cast<int64_t>(r) * cols + c;
                if (hasND && data[idx_i] == noData) continue;

                double zi = z[idx_i];

                // Iterate over neighbors within Chebyshev distance <= lag
                for (int dr = -lag; dr <= lag; ++dr) {
                    for (int dc = -lag; dc <= lag; ++dc) {
                        if (dr == 0 && dc == 0) continue; // skip self

                        int nr = r + dr;
                        int nc = c + dc;
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;

                        int64_t idx_j = static_cast<int64_t>(nr) * cols + nc;
                        if (hasND && data[idx_j] == noData) continue;

                        numerator += zi * z[idx_j];
                        W += 1.0;
                    }
                }
            }

            if (r % 100 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(r) / rows);
        }

        // Moran's I = (N / W) * (numerator / denominator)
        double moranI = (static_cast<double>(N) / W) * (numerator / denominator);

        // Expected value under null hypothesis of no spatial autocorrelation
        double expectedI = -1.0 / (N - 1);

        // Variance of I under normality assumption (randomization)
        // Var(I) = [N^2 * S1 - N * S2 + 3 * W^2] / [(W^2)(N^2 - 1)] - E(I)^2
        // For regular lattice with queen contiguity, we use a simpler approximation:
        // Under randomization: Var(I) approx 1/(N-1) for large N with regular weights
        // More precise: compute S1 and S2

        // S1 = 0.5 * sum_i sum_j (w_ij + w_ji)^2
        // For symmetric binary weights: w_ij + w_ji = 2 when connected, so S1 = 2 * W
        double S1 = 2.0 * W;

        // S2 = sum_i (sum_j w_ij + sum_j w_ji)^2
        // For symmetric binary weights, each connected cell i has degree d_i neighbors
        // so sum_j w_ij = d_i, and S2 = sum_i (2 * d_i)^2 = 4 * sum_i d_i^2
        // We need to compute sum of squared degrees
        double sumDegSq = 0.0;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx_i = static_cast<int64_t>(r) * cols + c;
                if (hasND && data[idx_i] == noData) continue;

                int degree = 0;
                for (int dr = -lag; dr <= lag; ++dr) {
                    for (int dc = -lag; dc <= lag; ++dc) {
                        if (dr == 0 && dc == 0) continue;
                        int nr = r + dr, nc = c + dc;
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                        int64_t idx_j = static_cast<int64_t>(nr) * cols + nc;
                        if (hasND && data[idx_j] == noData) continue;
                        ++degree;
                    }
                }
                sumDegSq += static_cast<double>(degree) * degree;
            }
        }
        double S2 = 4.0 * sumDegSq;

        double Nd = static_cast<double>(N);
        double W2 = W * W;
        double varianceI = (Nd * Nd * S1 - Nd * S2 + 3.0 * W2) /
                           (W2 * (Nd * Nd - 1.0)) - (expectedI * expectedI);

        double stddevI = std::sqrt(std::max(varianceI, 0.0));
        double zScore = (stddevI > 1e-15) ? (moranI - expectedI) / stddevI : 0.0;

        reportProgress(0.95, "Writing results...");

        // Write results to text file
        QString outputPath = parameter("output").toString();
        std::ofstream outFile(outputPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output file: " + outputPath);
            return false;
        }

        outFile << "Moran's I Spatial Autocorrelation Analysis\n";
        outFile << "==========================================\n\n";
        outFile << "Input: " << parameter("input").toString().toStdString() << "\n";
        outFile << "Lag distance: " << lag << " cell(s)\n\n";
        outFile << "Number of valid cells (N): " << N << "\n";
        outFile << "Total weight sum (W): " << W << "\n";
        outFile << "Mean value: " << mean << "\n\n";
        outFile << "Moran's I: " << moranI << "\n";
        outFile << "Expected I [E(I)]: " << expectedI << "\n";
        outFile << "Variance of I: " << varianceI << "\n";
        outFile << "Std. deviation of I: " << stddevI << "\n";
        outFile << "Z-score: " << zScore << "\n\n";

        if (zScore > 2.576)
            outFile << "Result: Significant positive spatial autocorrelation (p < 0.01)\n";
        else if (zScore > 1.96)
            outFile << "Result: Significant positive spatial autocorrelation (p < 0.05)\n";
        else if (zScore < -2.576)
            outFile << "Result: Significant negative spatial autocorrelation (p < 0.01)\n";
        else if (zScore < -1.96)
            outFile << "Result: Significant negative spatial autocorrelation (p < 0.05)\n";
        else
            outFile << "Result: No significant spatial autocorrelation at p = 0.05\n";

        outFile << "\nInterpretation:\n";
        outFile << "  I > E(I): Positive autocorrelation (clustering of similar values)\n";
        outFile << "  I < E(I): Negative autocorrelation (dispersion / checkerboard pattern)\n";
        outFile << "  I ~ E(I): Random spatial pattern\n";
        outFile.close();

        // Report summary via progress callback
        reportProgress(1.0,
            QString("Moran's I = %1, E(I) = %2, z-score = %3")
                .arg(moranI, 0, 'f', 6)
                .arg(expectedI, 0, 'f', 6)
                .arg(zScore, 0, 'f', 4));

        return true;
    }
};

REGISTER_MODULE(AutocorrModule)

} // namespace aplaceholder
