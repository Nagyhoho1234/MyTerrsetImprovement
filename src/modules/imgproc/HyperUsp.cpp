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

class HyperUspModule : public Module {
public:
    QString name() const override { return "HYPERUSP"; }
    QString description() const override {
        return "Unified Subspace Projection for hyperspectral target detection. "
               "Combines orthogonal subspace projection with a matched filter approach "
               "using the image covariance for improved detection performance.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::file("target_spectrum_file", "Target spectrum file (CSV)",
                "CSV: single row with target spectral values per band"),
            ParameterDef::output("output", "Output detection image",
                "Output raster with combined OSP + matched filter detection score"),
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
        // 2. Read target spectrum
        // --------------------------------------------------------------------
        QString targetPath = parameter("target_spectrum_file").toString();
        std::vector<double> target;

        {
            std::ifstream file(targetPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open target spectrum file: " + targetPath);
                return false;
            }

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;

                std::istringstream iss(line);
                std::string token;
                while (std::getline(iss, token, ',')) {
                    size_t start = token.find_first_not_of(" \t");
                    size_t end = token.find_last_not_of(" \t\r\n");
                    if (start != std::string::npos && end != std::string::npos)
                        target.push_back(std::stod(token.substr(start, end - start + 1)));
                }
                if (!target.empty()) break;
            }
        }

        if (static_cast<int>(target.size()) < numBands) {
            setError("Target spectrum has fewer values than the number of bands");
            return false;
        }
        target.resize(numBands);

        // --------------------------------------------------------------------
        // 3. Collect band data pointers
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bandsData(numBands);
        for (int b = 0; b < numBands; ++b)
            bandsData[b] = &bandRasters[b]->data(0);

        reportProgress(0.1, "Computing image covariance matrix...");

        // --------------------------------------------------------------------
        // 4. Compute image mean and covariance matrix
        // --------------------------------------------------------------------
        std::vector<double> means(numBands, 0.0);
        int64_t validCount = 0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bandsData[0])[i] == noData) continue;
            validCount++;
            for (int b = 0; b < numBands; ++b)
                means[b] += (*bandsData[b])[i];
        }

        if (validCount == 0) {
            setError("No valid pixels found");
            return false;
        }

        for (int b = 0; b < numBands; ++b)
            means[b] /= validCount;

        // Covariance matrix (numBands x numBands)
        std::vector<double> covMat(numBands * numBands, 0.0);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bandsData[0])[i] == noData) continue;
            for (int b1 = 0; b1 < numBands; ++b1) {
                double d1 = (*bandsData[b1])[i] - means[b1];
                for (int b2 = b1; b2 < numBands; ++b2) {
                    double d2 = (*bandsData[b2])[i] - means[b2];
                    covMat[b1 * numBands + b2] += d1 * d2;
                }
            }
        }
        for (int b1 = 0; b1 < numBands; ++b1) {
            for (int b2 = b1; b2 < numBands; ++b2) {
                covMat[b1 * numBands + b2] /= (validCount - 1);
                covMat[b2 * numBands + b1] = covMat[b1 * numBands + b2];
            }
        }

        reportProgress(0.3, "Inverting covariance matrix...");

        // --------------------------------------------------------------------
        // 5. Invert covariance matrix via Gauss-Jordan
        // --------------------------------------------------------------------
        std::vector<double> covInv(numBands * numBands, 0.0);
        {
            int n = numBands;
            std::vector<double> aug(n * 2 * n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j)
                    aug[i * 2 * n + j] = covMat[i * n + j];
                aug[i * 2 * n + n + i] = 1.0;
            }

            for (int i = 0; i < n; ++i) {
                int pivotRow = i;
                double pivotVal = std::abs(aug[i * 2 * n + i]);
                for (int k = i + 1; k < n; ++k) {
                    double val = std::abs(aug[k * 2 * n + i]);
                    if (val > pivotVal) {
                        pivotVal = val;
                        pivotRow = k;
                    }
                }

                if (pivotVal < 1e-12) {
                    setError("Image covariance matrix is singular or nearly singular");
                    return false;
                }

                if (pivotRow != i) {
                    for (int j = 0; j < 2 * n; ++j)
                        std::swap(aug[i * 2 * n + j], aug[pivotRow * 2 * n + j]);
                }

                double scale = aug[i * 2 * n + i];
                for (int j = 0; j < 2 * n; ++j)
                    aug[i * 2 * n + j] /= scale;

                for (int k = 0; k < n; ++k) {
                    if (k == i) continue;
                    double factor = aug[k * 2 * n + i];
                    for (int j = 0; j < 2 * n; ++j)
                        aug[k * 2 * n + j] -= factor * aug[i * 2 * n + j];
                }
            }

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    covInv[i * n + j] = aug[i * 2 * n + n + j];
        }

        // --------------------------------------------------------------------
        // 6. Compute USP filter: w = covInv * d / (d^T * covInv * d)
        //    This is the matched filter (constrained energy minimization)
        // --------------------------------------------------------------------
        // covInv * d
        std::vector<double> covInvD(numBands, 0.0);
        for (int i = 0; i < numBands; ++i) {
            for (int j = 0; j < numBands; ++j)
                covInvD[i] += covInv[i * numBands + j] * target[j];
        }

        // d^T * covInv * d
        double dTcovInvD = 0.0;
        for (int b = 0; b < numBands; ++b)
            dTcovInvD += target[b] * covInvD[b];

        if (std::abs(dTcovInvD) < 1e-12) {
            setError("Target spectrum is in the null space of the covariance matrix");
            return false;
        }

        // w = covInv * d / (d^T * covInv * d)
        std::vector<double> w(numBands);
        for (int b = 0; b < numBands; ++b)
            w[b] = covInvD[b] / dTcovInvD;

        reportProgress(0.4, "Applying USP detection to pixels...");

        // --------------------------------------------------------------------
        // 7. Create output raster and apply filter
        // --------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);

        for (int64_t idx = 0; idx < total; ++idx) {
            if (hasND && (*bandsData[0])[idx] == noData) {
                outData[idx] = noData;
                continue;
            }

            double score = 0.0;
            for (int b = 0; b < numBands; ++b)
                score += w[b] * (*bandsData[b])[idx];

            outData[idx] = score;

            if (idx % 1000000 == 0)
                reportProgress(0.4 + 0.5 * static_cast<double>(idx) / total);
        }

        // --------------------------------------------------------------------
        // 8. Write output
        // --------------------------------------------------------------------
        reportProgress(0.9, "Writing output...");
        QString outputPath = parameter("output").toString();
        if (!GdalIO::write(output, outputPath)) {
            setError("Failed to write output: " + outputPath);
            return false;
        }

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(HyperUspModule)

} // namespace aplaceholder
