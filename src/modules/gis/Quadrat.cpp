#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <fstream>

namespace aplaceholder {

class QuadratModule : public Module {
public:
    QString name() const override { return "QUADRAT"; }
    QString description() const override {
        return "Quadrat analysis for point pattern. Divides the study area into a grid, "
               "counts points per cell, and computes variance/mean ratio to test "
               "spatial randomness (VMR = 1 random, >1 clustered, <1 dispersed).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input point raster (non-zero = point)"),
            ParameterDef::output("output_file", "Output CSV file"),
            ParameterDef::integer("quadrat_size", "Quadrat size (pixels)", 10, 2, 10000,
                "Size of each quadrat cell in pixels"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) { setError("Failed to read input raster"); return false; }

        int cols = raster->cols(), rows = raster->rows();
        const auto& data = raster->data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();
        int qSize = parameter("quadrat_size").toInt();

        if (qSize < 2) {
            setError("Quadrat size must be at least 2");
            return false;
        }

        int qCols = (cols + qSize - 1) / qSize;
        int qRows = (rows + qSize - 1) / qSize;
        int numQuadrats = qCols * qRows;

        std::vector<int> counts(numQuadrats, 0);

        // Count points in each quadrat
        reportProgress(0.0, "Counting points per quadrat...");
        int64_t totalPoints = 0;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                double val = data[idx];

                if (hasND && val == noData) continue;
                if (val == 0) continue;

                int qr = r / qSize;
                int qc = c / qSize;
                counts[qr * qCols + qc]++;
                totalPoints++;
            }
            if (r % 100 == 0)
                reportProgress(0.7 * static_cast<double>(r) / rows);
        }

        // Compute statistics
        reportProgress(0.7, "Computing statistics...");
        double sum = 0.0;
        for (int i = 0; i < numQuadrats; ++i)
            sum += counts[i];

        double mean = sum / numQuadrats;

        double sumSqDev = 0.0;
        for (int i = 0; i < numQuadrats; ++i) {
            double dev = counts[i] - mean;
            sumSqDev += dev * dev;
        }
        double variance = (numQuadrats > 1) ? sumSqDev / (numQuadrats - 1) : 0.0;
        double vmr = (mean > 0) ? variance / mean : 0.0;

        // Chi-square statistic: (n-1) * VMR, df = n-1
        double chiSquare = (numQuadrats - 1) * vmr;

        // Write output CSV
        reportProgress(0.9, "Writing output...");
        QString outPath = parameter("output_file").toString();
        std::ofstream ofs(outPath.toStdString());
        if (!ofs.is_open()) {
            setError("Failed to open output file: " + outPath);
            return false;
        }

        ofs << "Quadrat Analysis Results\n";
        ofs << "Parameter,Value\n";
        ofs << "Raster columns," << cols << "\n";
        ofs << "Raster rows," << rows << "\n";
        ofs << "Quadrat size (pixels)," << qSize << "\n";
        ofs << "Number of quadrats," << numQuadrats << "\n";
        ofs << "Quadrat grid (cols x rows)," << qCols << " x " << qRows << "\n";
        ofs << "Total points," << totalPoints << "\n";
        ofs << "Mean points per quadrat," << mean << "\n";
        ofs << "Variance," << variance << "\n";
        ofs << "Variance/Mean ratio (VMR)," << vmr << "\n";
        ofs << "Chi-square statistic," << chiSquare << "\n";
        ofs << "Degrees of freedom," << (numQuadrats - 1) << "\n";
        ofs << "\nInterpretation:\n";
        ofs << "VMR = 1: Complete spatial randomness (Poisson)\n";
        ofs << "VMR > 1: Clustered pattern\n";
        ofs << "VMR < 1: Dispersed/regular pattern\n";
        ofs << "\nQuadrat Counts\n";
        ofs << "QuadratRow,QuadratCol,Count\n";
        for (int qr = 0; qr < qRows; ++qr) {
            for (int qc = 0; qc < qCols; ++qc) {
                ofs << qr << "," << qc << "," << counts[qr * qCols + qc] << "\n";
            }
        }

        ofs.close();
        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(QuadratModule)

} // namespace aplaceholder
