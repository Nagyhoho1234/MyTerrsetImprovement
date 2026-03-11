#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <random>

namespace aplaceholder {

class HyperAutoSigModule : public Module {
public:
    QString name() const override { return "HYPERAUTOSIG"; }
    QString description() const override {
        return "Automatic endmember extraction from hyperspectral data. Uses a combination "
               "of N-FINDR and Pixel Purity Index (PPI) algorithms to identify spectrally "
               "pure endmember pixels.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::output("output_endmember_file", "Output endmember file (CSV)",
                "CSV output: each row is an extracted endmember, columns are band values"),
            ParameterDef::integer("num_endmembers", "Number of endmembers", 5, 2, 256,
                "Number of endmembers to extract"),
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

        int numEM = parameter("num_endmembers").toInt();
        if (numEM > numBands) {
            setError(QString("Number of endmembers (%1) cannot exceed number of bands (%2)")
                     .arg(numEM).arg(numBands));
            return false;
        }

        // --------------------------------------------------------------------
        // 2. Collect band data and build valid pixel list
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        std::vector<int64_t> validPixels;
        validPixels.reserve(total);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) continue;
            validPixels.push_back(i);
        }

        int64_t numValid = static_cast<int64_t>(validPixels.size());
        if (numValid < numEM) {
            setError("Not enough valid pixels to extract the requested number of endmembers");
            return false;
        }

        reportProgress(0.05, "Running Pixel Purity Index (PPI)...");

        // --------------------------------------------------------------------
        // 3. PPI: Project all pixels onto random unit vectors and count
        //    how often each pixel appears as an extreme (min or max)
        // --------------------------------------------------------------------
        std::vector<int> ppiCount(numValid, 0);
        std::mt19937 rng(42); // Fixed seed for reproducibility
        std::normal_distribution<double> normDist(0.0, 1.0);

        int numProjections = std::min(static_cast<int64_t>(10000), numValid * 2);
        std::vector<double> projVec(numBands);

        for (int p = 0; p < numProjections; ++p) {
            // Generate random unit vector
            double mag = 0.0;
            for (int b = 0; b < numBands; ++b) {
                projVec[b] = normDist(rng);
                mag += projVec[b] * projVec[b];
            }
            mag = std::sqrt(mag);
            if (mag < 1e-12) continue;
            for (int b = 0; b < numBands; ++b)
                projVec[b] /= mag;

            // Project all valid pixels
            double minProj = std::numeric_limits<double>::max();
            double maxProj = -std::numeric_limits<double>::max();
            int64_t minIdx = 0, maxIdx = 0;

            for (int64_t v = 0; v < numValid; ++v) {
                int64_t pixIdx = validPixels[v];
                double proj = 0.0;
                for (int b = 0; b < numBands; ++b)
                    proj += projVec[b] * (*bands[b])[pixIdx];

                if (proj < minProj) { minProj = proj; minIdx = v; }
                if (proj > maxProj) { maxProj = proj; maxIdx = v; }
            }

            ppiCount[minIdx]++;
            ppiCount[maxIdx]++;

            if (p % 1000 == 0)
                reportProgress(0.05 + 0.30 * static_cast<double>(p) / numProjections);
        }

        reportProgress(0.35, "Selecting candidate endmembers...");

        // --------------------------------------------------------------------
        // 4. Select top PPI candidates as initial endmembers for N-FINDR
        // --------------------------------------------------------------------
        // Sort by PPI count descending, take top candidates
        int numCandidates = std::min(static_cast<int64_t>(numEM * 20), numValid);
        std::vector<int64_t> candidateOrder(numValid);
        for (int64_t i = 0; i < numValid; ++i)
            candidateOrder[i] = i;

        std::partial_sort(candidateOrder.begin(),
                          candidateOrder.begin() + numCandidates,
                          candidateOrder.end(),
                          [&](int64_t a, int64_t b) {
                              return ppiCount[a] > ppiCount[b];
                          });

        // Initialize endmember set with top PPI pixels
        std::vector<int64_t> emIndices(numEM);
        for (int e = 0; e < numEM; ++e)
            emIndices[e] = validPixels[candidateOrder[e]];

        reportProgress(0.4, "Refining endmembers with N-FINDR...");

        // --------------------------------------------------------------------
        // 5. N-FINDR: Iteratively replace endmembers to maximize simplex volume
        // --------------------------------------------------------------------
        // Helper: compute volume of simplex formed by numEM points in numBands-space
        // Using the determinant of the augmented matrix approach
        // Volume proportional to |det([e1..eP; 1..1])| where P = numEM
        auto computeVolume = [&](const std::vector<int64_t>& emIdx) -> double {
            int p = static_cast<int>(emIdx.size());
            // Build p x p matrix (use first p-1 bands + row of 1s)
            // For general case, build the matrix M where M[i][j] = band j value of endmember i
            // and add a row of 1s, then compute determinant
            int dim = p;
            std::vector<double> mat(dim * dim, 0.0);
            for (int i = 0; i < p; ++i) {
                for (int j = 0; j < p - 1; ++j) {
                    int bandIdx = j < numBands ? j : 0;
                    mat[i * dim + j] = (*bands[bandIdx])[emIdx[i]];
                }
                mat[i * dim + (p - 1)] = 1.0;
            }

            // Compute determinant via LU decomposition (partial pivoting)
            std::vector<double> lu(mat);
            double det = 1.0;
            for (int i = 0; i < dim; ++i) {
                // Find pivot
                int pivotRow = i;
                double pivotVal = std::abs(lu[i * dim + i]);
                for (int k = i + 1; k < dim; ++k) {
                    double val = std::abs(lu[k * dim + i]);
                    if (val > pivotVal) {
                        pivotVal = val;
                        pivotRow = k;
                    }
                }

                if (pivotVal < 1e-30) return 0.0;

                if (pivotRow != i) {
                    for (int j = 0; j < dim; ++j)
                        std::swap(lu[i * dim + j], lu[pivotRow * dim + j]);
                    det = -det;
                }

                det *= lu[i * dim + i];
                double invPivot = 1.0 / lu[i * dim + i];
                for (int k = i + 1; k < dim; ++k) {
                    lu[k * dim + i] *= invPivot;
                    for (int j = i + 1; j < dim; ++j)
                        lu[k * dim + j] -= lu[k * dim + i] * lu[i * dim + j];
                }
            }

            return std::abs(det);
        };

        // N-FINDR iterations
        const int maxNfindrIter = 50;
        for (int iter = 0; iter < maxNfindrIter; ++iter) {
            bool replaced = false;

            for (int c = 0; c < numCandidates; ++c) {
                int64_t candPix = validPixels[candidateOrder[c]];

                for (int e = 0; e < numEM; ++e) {
                    int64_t origPix = emIndices[e];
                    if (candPix == origPix) continue;

                    double origVol = computeVolume(emIndices);

                    emIndices[e] = candPix;
                    double newVol = computeVolume(emIndices);

                    if (newVol > origVol) {
                        replaced = true; // keep replacement
                    } else {
                        emIndices[e] = origPix; // revert
                    }
                }
            }

            if (!replaced) break;

            reportProgress(0.4 + 0.4 * static_cast<double>(iter + 1) / maxNfindrIter);
        }

        reportProgress(0.85, "Writing endmember file...");

        // --------------------------------------------------------------------
        // 6. Write endmembers to CSV
        // --------------------------------------------------------------------
        QString outputPath = parameter("output_endmember_file").toString();
        std::ofstream outFile(outputPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output file for writing: " + outputPath);
            return false;
        }

        outFile << "# Automatically extracted endmembers (" << numEM << " endmembers, "
                << numBands << " bands)\n";

        for (int e = 0; e < numEM; ++e) {
            int64_t pixIdx = emIndices[e];
            for (int b = 0; b < numBands; ++b) {
                if (b > 0) outFile << ",";
                outFile << (*bands[b])[pixIdx];
            }
            outFile << "\n";
        }

        outFile.close();

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(HyperAutoSigModule)

} // namespace aplaceholder
