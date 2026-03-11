#include "Module.h"
#include "ModuleRegistry.h"
#include <QFile>
#include <QTextStream>
#include <cmath>
#include <vector>

namespace aplaceholder {

class WeightModule : public Module {
public:
    QString name() const override { return "WEIGHT"; }
    QString description() const override {
        return "Calculates criterion weights for Multi-Criteria Evaluation using "
               "Saaty's Analytical Hierarchy Process (AHP) pairwise comparison method. "
               "Computes the principal eigenvector of a pairwise comparison matrix.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Pairwise comparison matrix file",
                "Text file containing N x N pairwise comparison matrix"),
            ParameterDef::output("output", "Output weights file",
                "Text file to write computed weights and consistency ratio"),
            ParameterDef::integer("num_criteria", "Number of criteria",
                3, 2, 15, "Number of criteria (rows/columns in the matrix)"),
        };
    }

    bool execute() override {
        int n = parameter("num_criteria").toInt();
        QString inputPath = parameter("input").toString();
        QString outputPath = parameter("output").toString();

        // Read pairwise comparison matrix
        std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));
        if (!readMatrix(inputPath, matrix, n)) {
            return false;
        }

        reportProgress(0.1, "Computing eigenvector...");

        // Compute principal eigenvector using the power method
        std::vector<double> weights(n, 1.0 / n);
        const int maxIterations = 1000;
        const double tolerance = 1e-10;

        for (int iter = 0; iter < maxIterations; ++iter) {
            // Multiply matrix by current vector
            std::vector<double> newWeights(n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    newWeights[i] += matrix[i][j] * weights[j];
                }
            }

            // Normalize to sum = 1
            double sum = 0.0;
            for (int i = 0; i < n; ++i) sum += newWeights[i];
            if (sum == 0.0) {
                setError("Degenerate matrix: eigenvector sum is zero");
                return false;
            }
            for (int i = 0; i < n; ++i) newWeights[i] /= sum;

            // Check convergence
            double maxDiff = 0.0;
            for (int i = 0; i < n; ++i) {
                maxDiff = std::max(maxDiff, std::abs(newWeights[i] - weights[i]));
            }

            weights = newWeights;
            if (maxDiff < tolerance) break;

            if (iter % 100 == 0)
                reportProgress(0.1 + 0.6 * (static_cast<double>(iter) / maxIterations));
        }

        reportProgress(0.7, "Computing consistency ratio...");

        // Compute principal eigenvalue (lambda_max)
        // lambda_max = sum of (A*w)_i / w_i for each i, divided by n
        double lambdaMax = 0.0;
        for (int i = 0; i < n; ++i) {
            double aw_i = 0.0;
            for (int j = 0; j < n; ++j) {
                aw_i += matrix[i][j] * weights[j];
            }
            if (weights[i] > 0.0) {
                lambdaMax += aw_i / weights[i];
            }
        }
        lambdaMax /= n;

        // Consistency Index
        double ci = (n > 1) ? (lambdaMax - n) / (n - 1) : 0.0;

        // Random Index table (Saaty, for matrices of size 1-15)
        static const double ri[] = {
            0.00, 0.00, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41,
            1.45, 1.49, 1.51, 1.48, 1.56, 1.57, 1.59
        };
        double riVal = (n >= 1 && n <= 15) ? ri[n - 1] : 1.59;

        // Consistency Ratio
        double cr = (riVal > 0.0) ? ci / riVal : 0.0;

        reportProgress(0.9, "Writing output...");

        // Write output file
        if (!writeOutput(outputPath, weights, n, lambdaMax, ci, cr)) {
            return false;
        }

        reportProgress(1.0, "Complete");
        return true;
    }

private:
    bool readMatrix(const QString& path, std::vector<std::vector<double>>& matrix, int n) {
        QFile file(path);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open matrix file: " + path);
            return false;
        }

        QTextStream in(&file);
        int row = 0;
        while (!in.atEnd() && row < n) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty()) continue;

            // Support both comma and whitespace delimiters
            QStringList parts;
            if (line.contains(',')) {
                parts = line.split(",", Qt::SkipEmptyParts);
            } else {
                // Split on whitespace
                parts = line.simplified().split(' ', Qt::SkipEmptyParts);
            }

            if (parts.size() < n) {
                setError("Row " + QString::number(row + 1) +
                         " has fewer than " + QString::number(n) + " values");
                return false;
            }

            for (int col = 0; col < n; ++col) {
                bool ok;
                // Handle fraction notation like "1/3"
                QString val = parts[col].trimmed();
                if (val.contains('/')) {
                    QStringList frac = val.split('/');
                    if (frac.size() == 2) {
                        double num = frac[0].toDouble(&ok);
                        if (!ok) {
                            setError("Invalid numerator in fraction at row " +
                                     QString::number(row + 1));
                            return false;
                        }
                        double den = frac[1].toDouble(&ok);
                        if (!ok || den == 0.0) {
                            setError("Invalid denominator in fraction at row " +
                                     QString::number(row + 1));
                            return false;
                        }
                        matrix[row][col] = num / den;
                    }
                } else {
                    matrix[row][col] = val.toDouble(&ok);
                    if (!ok) {
                        setError("Invalid value at row " + QString::number(row + 1) +
                                 ", col " + QString::number(col + 1));
                        return false;
                    }
                }
            }
            ++row;
        }

        if (row < n) {
            setError("Matrix file has fewer than " + QString::number(n) + " rows");
            return false;
        }

        return true;
    }

    bool writeOutput(const QString& path, const std::vector<double>& weights,
                     int n, double lambdaMax, double ci, double cr) {
        QFile file(path);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file for writing: " + path);
            return false;
        }

        QTextStream out(&file);
        out << "AHP Weight Calculation Results\n";
        out << "==============================\n\n";
        out << "Number of criteria: " << n << "\n\n";

        out << "Weights:\n";
        for (int i = 0; i < n; ++i) {
            out << "  Criterion " << (i + 1) << ": "
                << QString::number(weights[i], 'f', 6) << "\n";
        }

        out << "\nPrincipal eigenvalue (lambda_max): "
            << QString::number(lambdaMax, 'f', 6) << "\n";
        out << "Consistency Index (CI): "
            << QString::number(ci, 'f', 6) << "\n";
        out << "Consistency Ratio (CR): "
            << QString::number(cr, 'f', 6) << "\n\n";

        if (cr > 0.10) {
            out << "WARNING: CR exceeds 0.10. The pairwise comparison matrix\n"
                << "is inconsistent and should be re-evaluated.\n";
        } else {
            out << "CR is acceptable (< 0.10). Weights are consistent.\n";
        }

        return true;
    }
};

REGISTER_MODULE(WeightModule)

} // namespace aplaceholder
