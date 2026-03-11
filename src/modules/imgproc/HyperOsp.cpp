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

class HyperOspModule : public Module {
public:
    QString name() const override { return "HYPEROSP"; }
    QString description() const override {
        return "Orthogonal Subspace Projection for hyperspectral target detection. "
               "Projects out known background endmember signatures to isolate the target signal.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::file("target_spectrum_file", "Target spectrum file (CSV)",
                "CSV: single row with target spectral values per band"),
            ParameterDef::file("background_endmembers_file", "Background endmembers file (CSV)",
                "CSV: each row is a background endmember, columns are band values"),
            ParameterDef::output("output", "Output detection image",
                "Output raster expressing degree of support for target presence"),
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
        // 2. Read target spectrum (CSV, single row)
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
                if (!target.empty()) break; // use first valid row
            }
        }

        if (static_cast<int>(target.size()) < numBands) {
            setError("Target spectrum has fewer values than the number of bands");
            return false;
        }
        target.resize(numBands);

        // --------------------------------------------------------------------
        // 3. Read background endmembers (CSV)
        // --------------------------------------------------------------------
        QString bgPath = parameter("background_endmembers_file").toString();
        std::vector<std::vector<double>> bgEndmembers;

        {
            std::ifstream file(bgPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open background endmembers file: " + bgPath);
                return false;
            }

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;

                std::istringstream iss(line);
                std::vector<double> values;
                std::string token;
                while (std::getline(iss, token, ',')) {
                    size_t start = token.find_first_not_of(" \t");
                    size_t end = token.find_last_not_of(" \t\r\n");
                    if (start != std::string::npos && end != std::string::npos)
                        values.push_back(std::stod(token.substr(start, end - start + 1)));
                }

                if (static_cast<int>(values.size()) >= numBands) {
                    values.resize(numBands);
                    bgEndmembers.push_back(std::move(values));
                }
            }
        }

        int numBG = static_cast<int>(bgEndmembers.size());
        if (numBG == 0) {
            setError("No valid background endmembers found");
            return false;
        }

        reportProgress(0.05, "Computing orthogonal subspace projection...");

        // --------------------------------------------------------------------
        // 4. Construct background matrix U (numBands x numBG) and compute
        //    projector P_perp = I - U * (U^T U)^{-1} * U^T
        // --------------------------------------------------------------------
        // U^T U (numBG x numBG)
        std::vector<double> UtU(numBG * numBG, 0.0);
        for (int i = 0; i < numBG; ++i) {
            for (int j = 0; j < numBG; ++j) {
                double sum = 0.0;
                for (int b = 0; b < numBands; ++b)
                    sum += bgEndmembers[i][b] * bgEndmembers[j][b];
                UtU[i * numBG + j] = sum;
            }
        }

        // Invert UtU via Gauss-Jordan
        std::vector<double> UtUinv(numBG * numBG, 0.0);
        {
            std::vector<double> aug(numBG * 2 * numBG, 0.0);
            for (int i = 0; i < numBG; ++i) {
                for (int j = 0; j < numBG; ++j)
                    aug[i * 2 * numBG + j] = UtU[i * numBG + j];
                aug[i * 2 * numBG + numBG + i] = 1.0;
            }

            for (int i = 0; i < numBG; ++i) {
                int pivotRow = i;
                double pivotVal = std::abs(aug[i * 2 * numBG + i]);
                for (int k = i + 1; k < numBG; ++k) {
                    double val = std::abs(aug[k * 2 * numBG + i]);
                    if (val > pivotVal) {
                        pivotVal = val;
                        pivotRow = k;
                    }
                }

                if (pivotVal < 1e-12) {
                    setError("Background endmember matrix is singular");
                    return false;
                }

                if (pivotRow != i) {
                    for (int j = 0; j < 2 * numBG; ++j)
                        std::swap(aug[i * 2 * numBG + j], aug[pivotRow * 2 * numBG + j]);
                }

                double scale = aug[i * 2 * numBG + i];
                for (int j = 0; j < 2 * numBG; ++j)
                    aug[i * 2 * numBG + j] /= scale;

                for (int k = 0; k < numBG; ++k) {
                    if (k == i) continue;
                    double factor = aug[k * 2 * numBG + i];
                    for (int j = 0; j < 2 * numBG; ++j)
                        aug[k * 2 * numBG + j] -= factor * aug[i * 2 * numBG + j];
                }
            }

            for (int i = 0; i < numBG; ++i)
                for (int j = 0; j < numBG; ++j)
                    UtUinv[i * numBG + j] = aug[i * 2 * numBG + numBG + j];
        }

        // Compute P_perp * d (the OSP operator applied to target)
        // P_perp = I - U * (U^T U)^{-1} * U^T
        // We need: q = P_perp * d  (numBands vector)
        // First compute: (U^T U)^{-1} * U^T * d  (numBG vector)
        std::vector<double> Utd(numBG, 0.0);
        for (int i = 0; i < numBG; ++i) {
            for (int b = 0; b < numBands; ++b)
                Utd[i] += bgEndmembers[i][b] * target[b];
        }

        std::vector<double> invUtd(numBG, 0.0);
        for (int i = 0; i < numBG; ++i) {
            for (int j = 0; j < numBG; ++j)
                invUtd[i] += UtUinv[i * numBG + j] * Utd[j];
        }

        // U * invUtd (numBands vector)
        std::vector<double> UinvUtd(numBands, 0.0);
        for (int b = 0; b < numBands; ++b) {
            for (int i = 0; i < numBG; ++i)
                UinvUtd[b] += bgEndmembers[i][b] * invUtd[i];
        }

        // q = d - U * (U^T U)^{-1} * U^T * d
        std::vector<double> q(numBands);
        for (int b = 0; b < numBands; ++b)
            q[b] = target[b] - UinvUtd[b];

        reportProgress(0.1, "Applying OSP to pixels...");

        // --------------------------------------------------------------------
        // 5. Collect band data pointers
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bandsData(numBands);
        for (int b = 0; b < numBands; ++b)
            bandsData[b] = &bandRasters[b]->data(0);

        // --------------------------------------------------------------------
        // 6. Create output raster
        // --------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);

        // --------------------------------------------------------------------
        // 7. For each pixel, compute OSP score: q^T * pixel
        // --------------------------------------------------------------------
        for (int64_t idx = 0; idx < total; ++idx) {
            if (hasND && (*bandsData[0])[idx] == noData) {
                outData[idx] = noData;
                continue;
            }

            double score = 0.0;
            for (int b = 0; b < numBands; ++b)
                score += q[b] * (*bandsData[b])[idx];

            outData[idx] = score;

            if (idx % 1000000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(idx) / total);
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

REGISTER_MODULE(HyperOspModule)

} // namespace aplaceholder
