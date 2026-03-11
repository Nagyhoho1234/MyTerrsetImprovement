#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <cmath>

namespace aplaceholder {

class MarkovChainModule : public Module {
public:
    QString name() const override { return "MARKOV_CHAIN"; }
    QString description() const override {
        return "Computes a Markov chain transition probability matrix from two "
               "land cover images of different dates. Projects future transition "
               "probabilities to a specified prediction date.";
    }
    QString category() const override { return "Land Change Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("earlier_image", "Earlier land cover image"),
            ParameterDef::file("later_image", "Later land cover image"),
            ParameterDef::output("output_matrix", "Output transition probability matrix"),
            ParameterDef::integer("prediction_date", "Prediction date (years from later image)", 10, 1, 500,
                "Number of time periods to project forward"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("earlier_image").toString());
        auto r2 = GdalIO::read(parameter("later_image").toString());
        if (!r1 || !r2) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        if (cols != r2->cols() || rows != r2->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        const auto& d1 = r1->data(0);
        const auto& d2 = r2->data(0);
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int predPeriods = parameter("prediction_date").toInt();

        // Collect all unique classes
        std::set<int> classSet;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (d1[i] == noData || d2[i] == noData))
                continue;
            classSet.insert(static_cast<int>(d1[i]));
            classSet.insert(static_cast<int>(d2[i]));
        }

        std::vector<int> classes(classSet.begin(), classSet.end());
        int n = static_cast<int>(classes.size());
        if (n == 0) {
            setError("No valid pixels found");
            return false;
        }

        // Map class values to indices
        std::map<int, int> classIndex;
        for (int i = 0; i < n; ++i)
            classIndex[classes[i]] = i;

        // Count transitions
        reportProgress(0.0, "Counting transitions...");
        std::vector<std::vector<int64_t>> counts(n, std::vector<int64_t>(n, 0));
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (d1[i] == noData || d2[i] == noData))
                continue;
            int from = classIndex[static_cast<int>(d1[i])];
            int to = classIndex[static_cast<int>(d2[i])];
            counts[from][to]++;

            if (i % 1000000 == 0)
                reportProgress(0.3 * static_cast<double>(i) / total);
        }

        // Normalize rows to get single-step transition probability matrix
        reportProgress(0.3, "Computing transition probabilities...");
        std::vector<std::vector<double>> P(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            int64_t rowSum = 0;
            for (int j = 0; j < n; ++j)
                rowSum += counts[i][j];
            if (rowSum > 0) {
                for (int j = 0; j < n; ++j)
                    P[i][j] = static_cast<double>(counts[i][j]) / rowSum;
            } else {
                // If class had no pixels at t1, set identity (persists)
                P[i][i] = 1.0;
            }
        }

        // Matrix exponentiation: P^t using repeated squaring
        reportProgress(0.5, "Projecting to future date...");
        std::vector<std::vector<double>> result = matrixPow(P, predPeriods, n);

        // Write output CSV
        reportProgress(0.8, "Writing transition matrix...");
        QString outputPath = parameter("output_matrix").toString();
        std::ofstream csv(outputPath.toStdString());
        if (!csv.is_open()) {
            setError("Failed to open output file for writing");
            return false;
        }

        // Header row
        csv << "From\\To";
        for (int j = 0; j < n; ++j)
            csv << "," << classes[j];
        csv << "\n";

        // Data rows
        for (int i = 0; i < n; ++i) {
            csv << classes[i];
            for (int j = 0; j < n; ++j)
                csv << "," << result[i][j];
            csv << "\n";
        }

        csv.close();

        reportProgress(1.0, "Markov chain analysis complete.");
        return true;
    }

private:
    // Matrix multiplication for square matrices
    std::vector<std::vector<double>> matrixMul(
        const std::vector<std::vector<double>>& A,
        const std::vector<std::vector<double>>& B,
        int n) const
    {
        std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k) {
                if (A[i][k] == 0.0) continue;
                for (int j = 0; j < n; ++j)
                    C[i][j] += A[i][k] * B[k][j];
            }
        return C;
    }

    // Matrix exponentiation by repeated squaring: M^t
    std::vector<std::vector<double>> matrixPow(
        const std::vector<std::vector<double>>& M,
        int t, int n) const
    {
        // Identity matrix
        std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i)
            result[i][i] = 1.0;

        std::vector<std::vector<double>> base = M;
        while (t > 0) {
            if (t & 1)
                result = matrixMul(result, base, n);
            base = matrixMul(base, base, n);
            t >>= 1;
        }
        return result;
    }
};

REGISTER_MODULE(MarkovChainModule)

} // namespace aplaceholder
