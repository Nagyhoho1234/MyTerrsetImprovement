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

class MahalClassModule : public Module {
public:
    QString name() const override { return "MAHALCLASS"; }
    QString description() const override {
        return "Mahalanobis Distance classifier. Assigns each pixel to the class "
               "with the smallest Mahalanobis distance, accounting for covariance "
               "structure of each class. Outputs typicality values (0-1) per class.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated paths)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::file("signature_file", "Signature file (CSV)"),
            ParameterDef::output("output", "Output classified image"),
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
        // 2. Read band rasters
        // ------------------------------------------------------------------
        reportProgress(0.0, "Reading band images...");
        std::vector<std::unique_ptr<Raster>> bandRasters(numBands);
        int cols = 0, rows = 0;

        for (int b = 0; b < numBands; ++b) {
            bandRasters[b] = GdalIO::read(bandPaths[b]);
            if (!bandRasters[b]) {
                setError("Failed to read band image: " + bandPaths[b]);
                return false;
            }
            if (b == 0) {
                cols = bandRasters[b]->cols();
                rows = bandRasters[b]->rows();
            } else if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                setError("Band dimensions do not match: " + bandPaths[b]);
                return false;
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = bandRasters[0]->hasNoData();
        double noData = bandRasters[0]->noDataValue();

        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // ------------------------------------------------------------------
        // 3. Read signature file (same format as MaxLike/MinDist)
        // ------------------------------------------------------------------
        QString sigPath = parameter("signature_file").toString();
        int n = numBands;

        struct ClassSig {
            int classId;
            std::vector<double> mean;
            std::vector<double> covMatrix;
            std::vector<double> covInverse;
        };

        std::vector<ClassSig> classes;

        {
            std::ifstream file(sigPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open signature file: " + sigPath);
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

                int expectedMin = 1 + n + n * n;
                if (static_cast<int>(values.size()) < expectedMin) continue;

                ClassSig sig;
                sig.classId = static_cast<int>(values[0]);
                sig.mean.resize(n);
                for (int b = 0; b < n; ++b)
                    sig.mean[b] = values[1 + b];
                sig.covMatrix.resize(n * n);
                for (int i = 0; i < n * n; ++i)
                    sig.covMatrix[i] = values[1 + n + i];

                classes.push_back(std::move(sig));
            }
        }

        if (classes.empty()) {
            setError("No valid class signatures found");
            return false;
        }

        // ------------------------------------------------------------------
        // 4. Precompute inverse covariance matrices using Gauss-Jordan
        // ------------------------------------------------------------------
        reportProgress(0.1, "Computing inverse covariance matrices...");

        for (auto& cls : classes) {
            std::vector<double> A(n * 2 * n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j)
                    A[i * 2 * n + j] = cls.covMatrix[i * n + j];
                A[i * 2 * n + n + i] = 1.0;
            }

            for (int i = 0; i < n; ++i) {
                int maxRow = i;
                for (int k = i + 1; k < n; ++k)
                    if (std::fabs(A[k * 2 * n + i]) > std::fabs(A[maxRow * 2 * n + i]))
                        maxRow = k;
                if (maxRow != i)
                    for (int j = 0; j < 2 * n; ++j)
                        std::swap(A[i * 2 * n + j], A[maxRow * 2 * n + j]);
                double pivot = A[i * 2 * n + i];
                if (std::fabs(pivot) < 1e-15) {
                    A[i * 2 * n + i] += 1e-6;
                    pivot = A[i * 2 * n + i];
                }
                for (int j = 0; j < 2 * n; ++j)
                    A[i * 2 * n + j] /= pivot;
                for (int k = 0; k < n; ++k) {
                    if (k == i) continue;
                    double factor = A[k * 2 * n + i];
                    for (int j = 0; j < 2 * n; ++j)
                        A[k * 2 * n + j] -= factor * A[i * 2 * n + j];
                }
            }

            cls.covInverse.resize(n * n);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    cls.covInverse[i * n + j] = A[i * 2 * n + n + j];
        }

        // ------------------------------------------------------------------
        // 5. Classify each pixel using Mahalanobis distance
        // ------------------------------------------------------------------
        reportProgress(0.2, "Classifying pixels...");

        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);
        std::vector<double> diff(n);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                out[i] = 0;
                continue;
            }

            double bestMahal = std::numeric_limits<double>::max();
            int bestClass = 0;

            for (const auto& cls : classes) {
                // Compute (x - mu)
                for (int b = 0; b < n; ++b)
                    diff[b] = (*bands[b])[i] - cls.mean[b];

                // Mahalanobis distance squared: (x-mu)^T * Sigma^{-1} * (x-mu)
                double mahal = 0.0;
                for (int b1 = 0; b1 < n; ++b1) {
                    double rowSum = 0.0;
                    for (int b2 = 0; b2 < n; ++b2)
                        rowSum += cls.covInverse[b1 * n + b2] * diff[b2];
                    mahal += diff[b1] * rowSum;
                }

                if (mahal < bestMahal) {
                    bestMahal = mahal;
                    bestClass = cls.classId;
                }
            }

            out[i] = bestClass;

            if (i % 1000000 == 0)
                reportProgress(0.2 + 0.7 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(MahalClassModule)

} // namespace aplaceholder
