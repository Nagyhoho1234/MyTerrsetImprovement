#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>

namespace aplaceholder {

class ScatterModule : public Module {
public:
    QString name() const override { return "SCATTER"; }
    QString description() const override {
        return "Generate scatter plot data from two raster images. "
               "Samples pairs of values and outputs as CSV for plotting.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First input raster (X axis)"),
            ParameterDef::file("input2", "Second input raster (Y axis)"),
            ParameterDef::output("output_file", "Output CSV file"),
            ParameterDef::integer("max_points", "Maximum number of sample points",
                10000, 100, 1000000,
                "Maximum number of pixel pairs to sample"),
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
        int maxPoints = parameter("max_points").toInt();
        if (maxPoints <= 0) maxPoints = 10000;

        reportProgress(0.0, "Collecting valid pixel indices...");

        // Collect valid pixel indices
        std::vector<int64_t> validIndices;
        validIndices.reserve(std::min(total, static_cast<int64_t>(maxPoints * 2)));
        for (int64_t i = 0; i < total; ++i) {
            if (hasND1 && d1[i] == noData1) continue;
            if (hasND2 && d2[i] == noData2) continue;
            validIndices.push_back(i);
        }

        if (validIndices.empty()) {
            setError("No valid pixel pairs found");
            return false;
        }

        reportProgress(0.3, "Sampling points...");

        // Sample if needed
        int64_t sampleSize = static_cast<int64_t>(validIndices.size());
        if (sampleSize > maxPoints) {
            std::mt19937_64 rng(42); // fixed seed for reproducibility
            std::shuffle(validIndices.begin(), validIndices.end(), rng);
            validIndices.resize(maxPoints);
            sampleSize = maxPoints;
        }

        reportProgress(0.5, "Writing CSV...");

        // Write CSV output
        QString outputPath = parameter("output_file").toString();
        std::ofstream outFile(outputPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output file: " + outputPath);
            return false;
        }

        outFile << "X,Y\n";
        for (int64_t idx = 0; idx < sampleSize; ++idx) {
            int64_t i = validIndices[idx];
            outFile << d1[i] << "," << d2[i] << "\n";

            if (idx % 10000 == 0)
                reportProgress(0.5 + 0.5 * static_cast<double>(idx) / sampleSize);
        }

        outFile.close();

        reportProgress(1.0,
            QString("Scatter data exported: %1 points from %2 valid pixels")
                .arg(sampleSize)
                .arg(static_cast<int64_t>(validIndices.size())));

        return true;
    }
};

REGISTER_MODULE(ScatterModule)

} // namespace aplaceholder
