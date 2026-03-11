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

class UnmixModule : public Module {
public:
    QString name() const override { return "UNMIX"; }
    QString description() const override {
        return "Linear spectral unmixing. Solves for endmember fractions per pixel using "
               "constrained least squares (sum-to-one, non-negativity).";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::file("endmember_file", "Endmember file (CSV)",
                "CSV: each row is an endmember, columns are band values"),
            ParameterDef::file("output_prefix", "Output prefix",
                "Prefix for fraction raster outputs (one per endmember + residual)"),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Read input bands
        // ------------------------------------------------------------------
        QString bandsParam = parameter("input_bands").toString();
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

        // ------------------------------------------------------------------
        // 2. Read endmember file
        // ------------------------------------------------------------------
        // CSV format: each row is an endmember, columns are spectral values per band
        QString emPath = parameter("endmember_file").toString();

        // endmembers[e][b] = spectral value of endmember e in band b
        std::vector<std::vector<double>> endmembers;

        {
            std::ifstream file(emPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open endmember file: " + emPath);
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
                    values.resize(numBands); // trim to numBands
                    endmembers.push_back(std::move(values));
                }
            }
        }

        int numEM = static_cast<int>(endmembers.size());
        if (numEM == 0) {
            setError("No valid endmembers found in endmember file");
            return false;
        }

        if (numEM > numBands) {
            setError(QString("Number of endmembers (%1) exceeds number of bands (%2). "
                             "Standard unmixing requires endmembers <= bands.")
                     .arg(numEM).arg(numBands));
            return false;
        }

        reportProgress(0.05, "Precomputing unmixing matrix...");

        // ------------------------------------------------------------------
        // 3. Precompute (E^T * E)^{-1} * E^T for least-squares solution
        // ------------------------------------------------------------------
        // E is numBands x numEM matrix (endmembers as columns)
        // pixel = E * f + residual
        // Unconstrained solution: f = (E^T E)^{-1} E^T pixel

        // Build E^T * E (numEM x numEM)
        std::vector<double> EtE(numEM * numEM, 0.0);
        for (int i = 0; i < numEM; ++i) {
            for (int j = 0; j < numEM; ++j) {
                double sum = 0.0;
                for (int b = 0; b < numBands; ++b)
                    sum += endmembers[i][b] * endmembers[j][b];
                EtE[i * numEM + j] = sum;
            }
        }

        // Invert EtE using Gauss-Jordan elimination
        std::vector<double> EtEinv(numEM * numEM, 0.0);
        {
            // Augmented matrix [EtE | I]
            std::vector<double> aug(numEM * 2 * numEM, 0.0);
            for (int i = 0; i < numEM; ++i) {
                for (int j = 0; j < numEM; ++j)
                    aug[i * 2 * numEM + j] = EtE[i * numEM + j];
                aug[i * 2 * numEM + numEM + i] = 1.0;
            }

            for (int i = 0; i < numEM; ++i) {
                // Find pivot
                int pivotRow = i;
                double pivotVal = std::abs(aug[i * 2 * numEM + i]);
                for (int k = i + 1; k < numEM; ++k) {
                    double val = std::abs(aug[k * 2 * numEM + i]);
                    if (val > pivotVal) {
                        pivotVal = val;
                        pivotRow = k;
                    }
                }

                if (pivotVal < 1e-12) {
                    setError("Endmember matrix is singular or nearly singular");
                    return false;
                }

                // Swap rows
                if (pivotRow != i) {
                    for (int j = 0; j < 2 * numEM; ++j)
                        std::swap(aug[i * 2 * numEM + j], aug[pivotRow * 2 * numEM + j]);
                }

                // Scale pivot row
                double scale = aug[i * 2 * numEM + i];
                for (int j = 0; j < 2 * numEM; ++j)
                    aug[i * 2 * numEM + j] /= scale;

                // Eliminate column
                for (int k = 0; k < numEM; ++k) {
                    if (k == i) continue;
                    double factor = aug[k * 2 * numEM + i];
                    for (int j = 0; j < 2 * numEM; ++j)
                        aug[k * 2 * numEM + j] -= factor * aug[i * 2 * numEM + j];
                }
            }

            // Extract inverse
            for (int i = 0; i < numEM; ++i)
                for (int j = 0; j < numEM; ++j)
                    EtEinv[i * numEM + j] = aug[i * 2 * numEM + numEM + j];
        }

        // Compute pseudoinverse: (E^T E)^{-1} * E^T  (numEM x numBands)
        std::vector<double> pinv(numEM * numBands, 0.0);
        for (int i = 0; i < numEM; ++i) {
            for (int b = 0; b < numBands; ++b) {
                double sum = 0.0;
                for (int k = 0; k < numEM; ++k)
                    sum += EtEinv[i * numEM + k] * endmembers[k][b];
                pinv[i * numBands + b] = sum;
            }
        }

        reportProgress(0.1, "Unmixing pixels...");

        // ------------------------------------------------------------------
        // 4. Collect band data pointers
        // ------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // ------------------------------------------------------------------
        // 5. Create output rasters (one per endmember + residual)
        // ------------------------------------------------------------------
        std::vector<Raster> fractionRasters;
        for (int e = 0; e <= numEM; ++e) { // +1 for residual
            fractionRasters.emplace_back(cols, rows, 1, DataType::Float64);
            fractionRasters.back().setGeoTransform(bandRasters[0]->geoTransform());
            fractionRasters.back().setProjection(bandRasters[0]->projection());
            fractionRasters.back().setNoDataValue(noData);
        }

        // Get data pointers
        std::vector<std::vector<double>*> fracData(numEM + 1);
        for (int e = 0; e <= numEM; ++e)
            fracData[e] = &fractionRasters[e].data(0);

        // ------------------------------------------------------------------
        // 6. For each pixel, solve for fractions with constraints
        // ------------------------------------------------------------------
        std::vector<double> pixelVals(numBands);
        std::vector<double> fractions(numEM);

        for (int64_t idx = 0; idx < total; ++idx) {
            if (hasND && (*bands[0])[idx] == noData) {
                for (int e = 0; e <= numEM; ++e)
                    (*fracData[e])[idx] = noData;
                continue;
            }

            // Get pixel spectral values
            for (int b = 0; b < numBands; ++b)
                pixelVals[b] = (*bands[b])[idx];

            // Unconstrained least-squares solution: f = pinv * pixel
            for (int e = 0; e < numEM; ++e) {
                double sum = 0.0;
                for (int b = 0; b < numBands; ++b)
                    sum += pinv[e * numBands + b] * pixelVals[b];
                fractions[e] = sum;
            }

            // Apply non-negativity constraint
            for (int e = 0; e < numEM; ++e) {
                if (fractions[e] < 0.0) fractions[e] = 0.0;
            }

            // Normalize to sum-to-one
            double sumF = 0.0;
            for (int e = 0; e < numEM; ++e)
                sumF += fractions[e];

            if (sumF > 1e-12) {
                for (int e = 0; e < numEM; ++e)
                    fractions[e] /= sumF;
            }

            // Store fractions
            for (int e = 0; e < numEM; ++e)
                (*fracData[e])[idx] = fractions[e];

            // Compute residual: sum of squared differences between observed and modeled
            double residual = 0.0;
            for (int b = 0; b < numBands; ++b) {
                double modeled = 0.0;
                for (int e = 0; e < numEM; ++e)
                    modeled += fractions[e] * endmembers[e][b];
                double diff = pixelVals[b] - modeled;
                residual += diff * diff;
            }
            (*fracData[numEM])[idx] = std::sqrt(residual);

            if (idx % 1000000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(idx) / total);
        }

        // ------------------------------------------------------------------
        // 7. Write output rasters
        // ------------------------------------------------------------------
        reportProgress(0.9, "Writing output rasters...");
        QString prefix = parameter("output_prefix").toString();

        for (int e = 0; e < numEM; ++e) {
            QString path = prefix + QString("_em%1.tif").arg(e + 1);
            if (!GdalIO::write(fractionRasters[e], path)) {
                setError("Failed to write fraction raster: " + path);
                return false;
            }
        }

        // Write residual
        QString resPath = prefix + "_residual.tif";
        if (!GdalIO::write(fractionRasters[numEM], resPath)) {
            setError("Failed to write residual raster: " + resPath);
            return false;
        }

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(UnmixModule)

} // namespace aplaceholder
