#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <vector>

namespace aplaceholder {

class KendallModule : public Module {
public:
    QString name() const override { return "KENDALL"; }
    QString description() const override {
        return "Kendall rank correlation (tau-b) between two raster images. "
               "A non-parametric measure of association based on concordant "
               "and discordant pairs, with correction for ties.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First input raster image"),
            ParameterDef::file("input2", "Second input raster image"),
            ParameterDef::output("output_report", "Output report text file"),
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

        reportProgress(0.0, "Collecting valid pixel pairs...");

        // Collect valid pairs
        std::vector<double> x, y;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND1 && d1[i] == noData1) continue;
            if (hasND2 && d2[i] == noData2) continue;
            x.push_back(d1[i]);
            y.push_back(d2[i]);
        }

        int64_t N = static_cast<int64_t>(x.size());
        if (N < 2) {
            setError("Insufficient valid pixel pairs for Kendall correlation (need at least 2)");
            return false;
        }

        reportProgress(0.1, "Computing Kendall tau-b...");

        // Compute concordant, discordant pairs and ties
        // For large N, use O(N*N) — may be slow for very large rasters
        // Limit to first 50000 pairs if too many
        int64_t maxPairs = 50000;
        bool sampled = false;
        if (N > maxPairs) {
            // Use systematic sampling
            double step = static_cast<double>(N) / maxPairs;
            std::vector<double> sx, sy;
            sx.reserve(maxPairs);
            sy.reserve(maxPairs);
            for (int64_t i = 0; i < maxPairs; ++i) {
                int64_t idx = static_cast<int64_t>(i * step);
                sx.push_back(x[idx]);
                sy.push_back(y[idx]);
            }
            x = std::move(sx);
            y = std::move(sy);
            N = maxPairs;
            sampled = true;
        }

        int64_t concordant = 0, discordant = 0;
        int64_t tiedX = 0, tiedY = 0, tiedXY = 0;

        for (int64_t i = 0; i < N - 1; ++i) {
            for (int64_t j = i + 1; j < N; ++j) {
                double dx = x[j] - x[i];
                double dy = y[j] - y[i];

                if (dx == 0.0 && dy == 0.0) {
                    ++tiedXY;
                } else if (dx == 0.0) {
                    ++tiedX;
                } else if (dy == 0.0) {
                    ++tiedY;
                } else if ((dx > 0 && dy > 0) || (dx < 0 && dy < 0)) {
                    ++concordant;
                } else {
                    ++discordant;
                }
            }

            if (i % 1000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / (N - 1));
        }

        // Kendall tau-b = (C - D) / sqrt((C+D+T_x)(C+D+T_y))
        // where T_x = pairs tied only on x, T_y = pairs tied only on y
        double n0 = concordant + discordant + tiedX + tiedY + tiedXY;
        double n1 = concordant + discordant + tiedX;  // pairs not tied on Y (or tied on both excluded)
        double n2 = concordant + discordant + tiedY;  // pairs not tied on X

        double tauB = 0.0;
        double denom = std::sqrt(static_cast<double>(n1) * n2);
        if (denom > 0.0) {
            tauB = static_cast<double>(concordant - discordant) / denom;
        }

        // Z-score for significance (for large N, tau ~ N(0, 2(2n+5)/(9n(n-1))))
        double varTau = (2.0 * (2.0 * N + 5.0)) / (9.0 * N * (N - 1.0));
        double zScore = (varTau > 0.0) ? tauB / std::sqrt(varTau) : 0.0;
        double pValue = std::erfc(std::abs(zScore) / std::sqrt(2.0));

        reportProgress(0.95, "Writing report...");

        // Write report
        QString reportPath = parameter("output_report").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        outFile << "Kendall Rank Correlation (tau-b) Analysis\n";
        outFile << "==========================================\n\n";
        outFile << "Input 1: " << parameter("input1").toString().toStdString() << "\n";
        outFile << "Input 2: " << parameter("input2").toString().toStdString() << "\n\n";
        outFile << "Number of pixel pairs (N): " << N << "\n";
        if (sampled) {
            outFile << "Note: Systematic sampling was applied (original N was larger)\n";
        }
        outFile << "\nConcordant pairs: " << concordant << "\n";
        outFile << "Discordant pairs: " << discordant << "\n";
        outFile << "Pairs tied on X only: " << tiedX << "\n";
        outFile << "Pairs tied on Y only: " << tiedY << "\n";
        outFile << "Pairs tied on both: " << tiedXY << "\n\n";
        outFile << "Kendall tau-b: " << tauB << "\n";
        outFile << "Z-score: " << zScore << "\n";
        outFile << "p-value (two-tailed, approx): " << pValue << "\n\n";

        if (pValue < 0.01)
            outFile << "Result: Significant correlation at p < 0.01\n";
        else if (pValue < 0.05)
            outFile << "Result: Significant correlation at p < 0.05\n";
        else
            outFile << "Result: No significant correlation at p = 0.05\n";

        outFile.close();

        reportProgress(1.0,
            QString("Kendall tau-b = %1, z = %2, p = %3")
                .arg(tauB, 0, 'f', 6)
                .arg(zScore, 0, 'f', 4)
                .arg(pValue, 0, 'g', 4));

        return true;
    }
};

REGISTER_MODULE(KendallModule)

} // namespace aplaceholder
