#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

namespace aplaceholder {

class MakeSigModule : public Module {
public:
    QString name() const override { return "MAKESIG"; }
    QString description() const override {
        return "Create Training Signatures. Computes per-class mean vectors and "
               "covariance matrices from multi-band rasters and a training site "
               "raster, producing signature files compatible with MAXLIKE.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated paths)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::file("training", "Training sites raster",
                "Integer raster with class IDs (0 = background)"),
            ParameterDef::output("output", "Output signature file",
                "CSV signature file compatible with MAXLIKE"),
        };
    }

    bool execute() override {
        // Parse comma-separated band file paths
        QString bandsParam = parameter("bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths)
            p = p.trimmed();

        if (bandPaths.isEmpty()) {
            setError("No band images specified");
            return false;
        }

        int numBands = bandPaths.size();

        // Read training sites raster
        auto training = GdalIO::read(parameter("training").toString());
        if (!training) {
            setError("Failed to read training sites raster");
            return false;
        }

        int cols = training->cols(), rows = training->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& trainData = training->data(0);
        bool hasNDTrain = training->hasNoData();
        double ndTrain = training->noDataValue();

        // Read all band rasters
        reportProgress(0.0, "Reading band images...");
        std::vector<std::unique_ptr<Raster>> bandRasters(numBands);
        for (int b = 0; b < numBands; ++b) {
            bandRasters[b] = GdalIO::read(bandPaths[b]);
            if (!bandRasters[b]) {
                setError("Failed to read band image: " + bandPaths[b]);
                return false;
            }
            if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                setError("Band image dimensions do not match training raster: " + bandPaths[b]);
                return false;
            }
        }

        // Collect band data pointers
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // Identify classes from training raster (non-zero, non-nodata)
        reportProgress(0.1, "Identifying classes and collecting pixel data...");

        struct ClassData {
            int classId;
            int64_t count = 0;
            std::vector<double> sum;       // length numBands
            std::vector<double> min;       // length numBands
            std::vector<double> max;       // length numBands
            // For covariance: we accumulate sum of (x_i * x_j) for each band pair
            std::vector<double> crossSum;  // numBands * numBands, row-major
        };

        std::map<int, ClassData> classMap;

        for (int64_t i = 0; i < total; ++i) {
            if (hasNDTrain && trainData[i] == ndTrain) continue;
            int classId = static_cast<int>(trainData[i]);
            if (classId == 0) continue; // background

            auto& cd = classMap[classId];
            if (cd.count == 0) {
                cd.classId = classId;
                cd.sum.resize(numBands, 0.0);
                cd.min.resize(numBands, std::numeric_limits<double>::max());
                cd.max.resize(numBands, -std::numeric_limits<double>::max());
                cd.crossSum.resize(numBands * numBands, 0.0);
            }

            for (int b = 0; b < numBands; ++b) {
                double val = (*bands[b])[i];
                cd.sum[b] += val;
                if (val < cd.min[b]) cd.min[b] = val;
                if (val > cd.max[b]) cd.max[b] = val;
            }

            // Accumulate cross products for covariance
            for (int b1 = 0; b1 < numBands; ++b1) {
                double v1 = (*bands[b1])[i];
                for (int b2 = b1; b2 < numBands; ++b2) {
                    double v2 = (*bands[b2])[i];
                    cd.crossSum[b1 * numBands + b2] += v1 * v2;
                    if (b1 != b2)
                        cd.crossSum[b2 * numBands + b1] += v1 * v2;
                }
            }

            cd.count++;

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.6 * static_cast<double>(i) / total);
        }

        if (classMap.empty()) {
            setError("No training pixels found (training raster has no non-zero values)");
            return false;
        }

        reportProgress(0.75, "Computing signatures...");

        // Compute mean and covariance for each class
        // Covariance: cov(b1,b2) = (1/N) * sum(x_b1 * x_b2) - mean_b1 * mean_b2
        // Using population covariance (consistent with typical remote sensing usage)

        struct ClassSignature {
            int classId;
            std::vector<double> mean;
            std::vector<double> covMatrix; // numBands x numBands, row-major
        };

        std::vector<ClassSignature> signatures;

        for (auto& [id, cd] : classMap) {
            if (cd.count < numBands) {
                reportProgress(-1.0, QString("Warning: Class %1 has only %2 pixels "
                    "(fewer than %3 bands). Signature may be unreliable.")
                    .arg(id).arg(cd.count).arg(numBands));
            }

            ClassSignature sig;
            sig.classId = id;
            sig.mean.resize(numBands);
            sig.covMatrix.resize(numBands * numBands);

            double N = static_cast<double>(cd.count);
            for (int b = 0; b < numBands; ++b)
                sig.mean[b] = cd.sum[b] / N;

            for (int b1 = 0; b1 < numBands; ++b1) {
                for (int b2 = 0; b2 < numBands; ++b2) {
                    sig.covMatrix[b1 * numBands + b2] =
                        cd.crossSum[b1 * numBands + b2] / N
                        - sig.mean[b1] * sig.mean[b2];
                }
            }

            signatures.push_back(std::move(sig));
        }

        reportProgress(0.9, "Writing signature file...");

        // Write signature file in CSV format compatible with MaxLike:
        // Each line: class_id, mean_b1, mean_b2, ..., mean_bN, cov_11, cov_12, ..., cov_NN
        // Lines starting with # are comments
        QString outPath = parameter("output").toString();
        std::ofstream outFile(outPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output signature file: " + outPath);
            return false;
        }

        outFile << "# Signature file generated by MAKESIG\n";
        outFile << "# Number of bands: " << numBands << "\n";
        outFile << "# Band files:";
        for (int b = 0; b < numBands; ++b)
            outFile << " " << bandPaths[b].toStdString();
        outFile << "\n";
        outFile << "# Number of classes: " << signatures.size() << "\n";
        outFile << "# Format: class_id, mean_b1..mean_bN, cov_11..cov_NN\n";

        outFile << std::fixed << std::setprecision(8);

        for (const auto& sig : signatures) {
            outFile << sig.classId;
            for (int b = 0; b < numBands; ++b)
                outFile << ", " << sig.mean[b];
            for (int i = 0; i < numBands * numBands; ++i)
                outFile << ", " << sig.covMatrix[i];
            outFile << "\n";
        }

        outFile.close();

        QString statsMsg = QString("Signature file created with %1 classes across %2 bands.")
                           .arg(signatures.size()).arg(numBands);
        reportProgress(1.0, statsMsg);
        return true;
    }
};

REGISTER_MODULE(MakeSigModule)

} // namespace aplaceholder
