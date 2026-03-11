#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class HyperSigModule : public Module {
public:
    QString name() const override { return "HYPERSIG"; }
    QString description() const override {
        return "Hyperspectral signature development from training areas. Creates full spectral "
               "profile signatures (mean, min, max, standard deviation per band) from training "
               "regions defined in a raster image.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::raster("training_raster", "Training site raster",
                "Raster image with integer class identifiers (0 = background)"),
            ParameterDef::output("output_sig_file", "Output signature file (CSV)",
                "CSV output with full spectral profiles per class"),
        };
    }

    bool execute() override {
        // --------------------------------------------------------------------
        // 1. Read input bands
        // --------------------------------------------------------------------
        QString bandsParam = parameter("bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths) p = p.trimmed();

        if (bandPaths.isEmpty()) {
            setError("No input bands specified");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> bandRasters;
        for (const auto& path : bandPaths) {
            auto r = GdalIO::read(path);
            if (!r) {
                setError("Failed to read band image: " + path);
                return false;
            }
            bandRasters.push_back(std::move(r));
        }

        int numBands = static_cast<int>(bandRasters.size());
        int cols = bandRasters[0]->cols();
        int rows = bandRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = bandRasters[0]->hasNoData();
        double noData = bandRasters[0]->noDataValue();

        for (int b = 1; b < numBands; ++b) {
            if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                setError("All input bands must have the same dimensions");
                return false;
            }
        }

        // --------------------------------------------------------------------
        // 2. Read training raster
        // --------------------------------------------------------------------
        QString trainPath = parameter("training_raster").toString();
        auto trainRaster = GdalIO::read(trainPath);
        if (!trainRaster) {
            setError("Failed to read training raster: " + trainPath);
            return false;
        }

        if (trainRaster->cols() != cols || trainRaster->rows() != rows) {
            setError("Training raster dimensions must match input band dimensions");
            return false;
        }

        const auto& trainData = trainRaster->data(0);

        reportProgress(0.1, "Collecting training statistics...");

        // --------------------------------------------------------------------
        // 3. Collect band data pointers
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // --------------------------------------------------------------------
        // 4. Identify unique classes and accumulate statistics
        // --------------------------------------------------------------------
        struct ClassStats {
            int classId = 0;
            int64_t count = 0;
            std::vector<double> sum;       // per band
            std::vector<double> sumSq;     // per band
            std::vector<double> minVal;    // per band
            std::vector<double> maxVal;    // per band
        };

        // Map from class ID to stats
        std::vector<ClassStats> classStats;
        auto findOrCreate = [&](int classId) -> ClassStats& {
            for (auto& cs : classStats) {
                if (cs.classId == classId) return cs;
            }
            classStats.push_back({});
            auto& cs = classStats.back();
            cs.classId = classId;
            cs.count = 0;
            cs.sum.assign(numBands, 0.0);
            cs.sumSq.assign(numBands, 0.0);
            cs.minVal.assign(numBands, std::numeric_limits<double>::max());
            cs.maxVal.assign(numBands, -std::numeric_limits<double>::max());
            return cs;
        };

        for (int64_t idx = 0; idx < total; ++idx) {
            int classId = static_cast<int>(trainData[idx]);
            if (classId <= 0) continue; // 0 = background, skip
            if (hasND && (*bands[0])[idx] == noData) continue;

            auto& cs = findOrCreate(classId);
            cs.count++;

            for (int b = 0; b < numBands; ++b) {
                double val = (*bands[b])[idx];
                cs.sum[b] += val;
                cs.sumSq[b] += val * val;
                if (val < cs.minVal[b]) cs.minVal[b] = val;
                if (val > cs.maxVal[b]) cs.maxVal[b] = val;
            }

            if (idx % 1000000 == 0)
                reportProgress(0.1 + 0.6 * static_cast<double>(idx) / total);
        }

        // Sort classes by ID
        std::sort(classStats.begin(), classStats.end(),
                  [](const ClassStats& a, const ClassStats& b) {
                      return a.classId < b.classId;
                  });

        if (classStats.empty()) {
            setError("No training pixels found (all class IDs are 0 or background)");
            return false;
        }

        reportProgress(0.8, "Writing signature file...");

        // --------------------------------------------------------------------
        // 5. Write output signature file (CSV)
        // --------------------------------------------------------------------
        QString outputPath = parameter("output_sig_file").toString();
        std::ofstream outFile(outputPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output file for writing: " + outputPath);
            return false;
        }

        outFile << "# Hyperspectral Signature File\n";
        outFile << "# Bands: " << numBands << "\n";
        outFile << "# Classes: " << classStats.size() << "\n";
        outFile << "#\n";

        for (const auto& cs : classStats) {
            outFile << "# Class " << cs.classId << " (n=" << cs.count << ")\n";

            // Mean
            outFile << "CLASS," << cs.classId << ",MEAN";
            for (int b = 0; b < numBands; ++b) {
                double mean = cs.sum[b] / cs.count;
                outFile << "," << mean;
            }
            outFile << "\n";

            // Min
            outFile << "CLASS," << cs.classId << ",MIN";
            for (int b = 0; b < numBands; ++b)
                outFile << "," << cs.minVal[b];
            outFile << "\n";

            // Max
            outFile << "CLASS," << cs.classId << ",MAX";
            for (int b = 0; b < numBands; ++b)
                outFile << "," << cs.maxVal[b];
            outFile << "\n";

            // StdDev
            outFile << "CLASS," << cs.classId << ",STDDEV";
            for (int b = 0; b < numBands; ++b) {
                double mean = cs.sum[b] / cs.count;
                double variance = (cs.sumSq[b] / cs.count) - (mean * mean);
                if (variance < 0.0) variance = 0.0;
                outFile << "," << std::sqrt(variance);
            }
            outFile << "\n";
        }

        outFile.close();

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(HyperSigModule)

} // namespace aplaceholder
