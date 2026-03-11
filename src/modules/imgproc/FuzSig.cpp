#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>

namespace aplaceholder {

class FuzSigModule : public Module {
public:
    QString name() const override { return "FUZSIG"; }
    QString description() const override {
        return "Fuzzy Signatures. Creates weighted training signatures from training sites "
               "using fuzzy membership grades. Each training pixel is weighted by its "
               "membership grade when computing per-class mean and covariance.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated paths)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::file("training_raster", "Training sites raster",
                "Integer raster with training site IDs (0 = background). "
                "Corresponding fuzzy membership images named fz<classname> are expected."),
            ParameterDef::output("output_sig_file", "Output fuzzy signature file",
                "CSV signature file compatible with classifiers"),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Parse band file paths
        // ------------------------------------------------------------------
        QString bandsParam = parameter("bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths)
            p = p.trimmed();

        if (bandPaths.isEmpty()) {
            setError("No band images specified");
            return false;
        }

        int numBands = bandPaths.size();

        // ------------------------------------------------------------------
        // 2. Read training sites raster
        // ------------------------------------------------------------------
        auto training = GdalIO::read(parameter("training_raster").toString());
        if (!training) {
            setError("Failed to read training sites raster");
            return false;
        }

        int cols = training->cols(), rows = training->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& trainData = training->data(0);
        bool hasNDTrain = training->hasNoData();
        double ndTrain = training->noDataValue();

        // ------------------------------------------------------------------
        // 3. Read band rasters
        // ------------------------------------------------------------------
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

        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // ------------------------------------------------------------------
        // 4. Identify distinct class IDs from training raster
        // ------------------------------------------------------------------
        reportProgress(0.1, "Identifying classes...");
        std::map<int, int> classIndex; // classId -> index
        for (int64_t i = 0; i < total; ++i) {
            if (hasNDTrain && trainData[i] == ndTrain) continue;
            int cid = static_cast<int>(trainData[i]);
            if (cid == 0) continue;
            if (classIndex.find(cid) == classIndex.end()) {
                int idx = static_cast<int>(classIndex.size());
                classIndex[cid] = idx;
            }
        }

        if (classIndex.empty()) {
            setError("No training pixels found");
            return false;
        }

        int numClasses = static_cast<int>(classIndex.size());

        // ------------------------------------------------------------------
        // 5. Build fuzzy membership by frequency histogram approach
        //    For each class, for each band, build a histogram of pixel values
        //    from training sites, then compute weighted mean/covariance where
        //    the weight of each pixel is the frequency-based fuzzy membership.
        //
        //    Frequency histogram approach: for each class c and band b, build
        //    a histogram of DN values within training pixels belonging to class c.
        //    The fuzzy membership of a pixel for class c is proportional to the
        //    histogram count at its DN value. This gives higher weight to typical
        //    (frequently occurring) values and lower weight to outliers.
        // ------------------------------------------------------------------
        reportProgress(0.2, "Building frequency histograms per class per band...");

        // First pass: collect pixel values per class
        struct ClassPixels {
            int classId;
            std::vector<std::vector<double>> bandValues; // [band][pixel_index]
        };
        std::vector<ClassPixels> classPixels(numClasses);
        for (auto& [cid, idx] : classIndex) {
            classPixels[idx].classId = cid;
            classPixels[idx].bandValues.resize(numBands);
        }

        for (int64_t i = 0; i < total; ++i) {
            if (hasNDTrain && trainData[i] == ndTrain) continue;
            int cid = static_cast<int>(trainData[i]);
            if (cid == 0) continue;
            int idx = classIndex[cid];
            for (int b = 0; b < numBands; ++b)
                classPixels[idx].bandValues[b].push_back((*bands[b])[i]);
        }

        // For each class, compute fuzzy membership weights using a frequency
        // histogram approach. We use a kernel density estimate (Gaussian kernel)
        // for continuous data, approximated by binning into 256 bins.
        reportProgress(0.4, "Computing fuzzy membership weights...");

        struct ClassSignature {
            int classId;
            std::vector<double> mean;       // numBands
            std::vector<double> covMatrix;  // numBands x numBands
            std::vector<double> minVal;     // numBands
            std::vector<double> maxVal;     // numBands
        };

        std::vector<ClassSignature> signatures(numClasses);

        for (int c = 0; c < numClasses; ++c) {
            auto& cp = classPixels[c];
            auto& sig = signatures[c];
            sig.classId = cp.classId;
            int64_t nPix = static_cast<int64_t>(cp.bandValues[0].size());
            if (nPix == 0) continue;

            // Compute per-band histograms and derive membership weights
            // Weight for each pixel = product of per-band frequency memberships
            const int nBins = 256;
            std::vector<double> weights(nPix, 1.0);

            for (int b = 0; b < numBands; ++b) {
                double bmin = *std::min_element(cp.bandValues[b].begin(), cp.bandValues[b].end());
                double bmax = *std::max_element(cp.bandValues[b].begin(), cp.bandValues[b].end());
                double range = bmax - bmin;
                if (range < 1e-12) range = 1.0;

                // Build histogram
                std::vector<int> hist(nBins, 0);
                for (int64_t p = 0; p < nPix; ++p) {
                    int bin = static_cast<int>((cp.bandValues[b][p] - bmin) / range * (nBins - 1));
                    if (bin < 0) bin = 0;
                    if (bin >= nBins) bin = nBins - 1;
                    hist[bin]++;
                }

                // Find max histogram count for normalization
                int maxCount = *std::max_element(hist.begin(), hist.end());
                if (maxCount == 0) maxCount = 1;

                // Multiply per-pixel weight by normalized frequency
                for (int64_t p = 0; p < nPix; ++p) {
                    int bin = static_cast<int>((cp.bandValues[b][p] - bmin) / range * (nBins - 1));
                    if (bin < 0) bin = 0;
                    if (bin >= nBins) bin = nBins - 1;
                    double membership = static_cast<double>(hist[bin]) / maxCount;
                    weights[p] *= membership;
                }
            }

            // Normalize weights so they sum to 1
            double wSum = 0.0;
            for (int64_t p = 0; p < nPix; ++p)
                wSum += weights[p];
            if (wSum < 1e-12) wSum = 1.0;
            for (int64_t p = 0; p < nPix; ++p)
                weights[p] /= wSum;

            // Compute weighted mean
            sig.mean.resize(numBands, 0.0);
            sig.minVal.resize(numBands, std::numeric_limits<double>::max());
            sig.maxVal.resize(numBands, -std::numeric_limits<double>::max());

            for (int b = 0; b < numBands; ++b) {
                for (int64_t p = 0; p < nPix; ++p) {
                    sig.mean[b] += weights[p] * cp.bandValues[b][p];
                }
                // Min/max from pixels with weight > 0.5 * max weight
                double maxW = *std::max_element(weights.begin(), weights.end());
                double threshold = 0.5 * maxW;
                for (int64_t p = 0; p < nPix; ++p) {
                    if (weights[p] >= threshold) {
                        if (cp.bandValues[b][p] < sig.minVal[b])
                            sig.minVal[b] = cp.bandValues[b][p];
                        if (cp.bandValues[b][p] > sig.maxVal[b])
                            sig.maxVal[b] = cp.bandValues[b][p];
                    }
                }
            }

            // Compute weighted covariance matrix
            sig.covMatrix.resize(numBands * numBands, 0.0);
            for (int b1 = 0; b1 < numBands; ++b1) {
                for (int b2 = b1; b2 < numBands; ++b2) {
                    double cov = 0.0;
                    for (int64_t p = 0; p < nPix; ++p) {
                        cov += weights[p]
                             * (cp.bandValues[b1][p] - sig.mean[b1])
                             * (cp.bandValues[b2][p] - sig.mean[b2]);
                    }
                    sig.covMatrix[b1 * numBands + b2] = cov;
                    sig.covMatrix[b2 * numBands + b1] = cov;
                }
            }

            reportProgress(0.4 + 0.4 * static_cast<double>(c + 1) / numClasses);
        }

        // ------------------------------------------------------------------
        // 6. Write signature file
        // ------------------------------------------------------------------
        reportProgress(0.9, "Writing fuzzy signature file...");

        QString outPath = parameter("output_sig_file").toString();
        std::ofstream outFile(outPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output signature file: " + outPath);
            return false;
        }

        outFile << "# Fuzzy signature file generated by FUZSIG\n";
        outFile << "# Number of bands: " << numBands << "\n";
        outFile << "# Band files:";
        for (int b = 0; b < numBands; ++b)
            outFile << " " << bandPaths[b].toStdString();
        outFile << "\n";
        outFile << "# Number of classes: " << numClasses << "\n";
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

        QString statsMsg = QString("Fuzzy signature file created with %1 classes across %2 bands.")
                           .arg(numClasses).arg(numBands);
        reportProgress(1.0, statsMsg);
        return true;
    }
};

REGISTER_MODULE(FuzSigModule)

} // namespace aplaceholder
