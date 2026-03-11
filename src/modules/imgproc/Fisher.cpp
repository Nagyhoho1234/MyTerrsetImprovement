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

class FisherModule : public Module {
public:
    QString name() const override { return "FISHER"; }
    QString description() const override {
        return "Linear Discriminant Analysis (Fisher) classifier. Projects multi-band data "
               "onto Fisher's discriminant axes that maximize between-class variance "
               "relative to within-class variance, then classifies to nearest class "
               "centroid in projected space.";
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
        // 3. Read signature file
        // ------------------------------------------------------------------
        QString sigPath = parameter("signature_file").toString();

        struct ClassSig {
            int classId;
            std::vector<double> mean;
            std::vector<double> covMatrix;
        };

        std::vector<ClassSig> classes;
        int n = numBands;

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

        int numClasses = static_cast<int>(classes.size());

        reportProgress(0.1, "Computing Fisher discriminant...");

        // ------------------------------------------------------------------
        // 4. Compute within-class scatter matrix (Sw) and between-class
        //    scatter matrix (Sb)
        // ------------------------------------------------------------------
        // Sw = sum of class covariance matrices (pooled)
        std::vector<double> Sw(n * n, 0.0);
        for (const auto& cls : classes)
            for (int i = 0; i < n * n; ++i)
                Sw[i] += cls.covMatrix[i];

        // Grand mean
        std::vector<double> grandMean(n, 0.0);
        for (const auto& cls : classes)
            for (int b = 0; b < n; ++b)
                grandMean[b] += cls.mean[b];
        for (int b = 0; b < n; ++b)
            grandMean[b] /= numClasses;

        // Sb = sum_c (mu_c - mu_grand)(mu_c - mu_grand)^T
        std::vector<double> Sb(n * n, 0.0);
        for (const auto& cls : classes) {
            std::vector<double> diff(n);
            for (int b = 0; b < n; ++b)
                diff[b] = cls.mean[b] - grandMean[b];
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    Sb[i * n + j] += diff[i] * diff[j];
        }

        // ------------------------------------------------------------------
        // 5. Compute Sw^{-1} * Sb using Gauss-Jordan on Sw
        // ------------------------------------------------------------------
        // Invert Sw
        std::vector<double> SwInv(n * 2 * n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                SwInv[i * 2 * n + j] = Sw[i * n + j];
            SwInv[i * 2 * n + n + i] = 1.0;
        }
        for (int i = 0; i < n; ++i) {
            // Partial pivoting
            int maxRow = i;
            for (int k = i + 1; k < n; ++k)
                if (std::fabs(SwInv[k * 2 * n + i]) > std::fabs(SwInv[maxRow * 2 * n + i]))
                    maxRow = k;
            if (maxRow != i)
                for (int j = 0; j < 2 * n; ++j)
                    std::swap(SwInv[i * 2 * n + j], SwInv[maxRow * 2 * n + j]);
            double pivot = SwInv[i * 2 * n + i];
            if (std::fabs(pivot) < 1e-15) {
                // Regularize
                SwInv[i * 2 * n + i] += 1e-6;
                pivot = SwInv[i * 2 * n + i];
            }
            for (int j = 0; j < 2 * n; ++j)
                SwInv[i * 2 * n + j] /= pivot;
            for (int k = 0; k < n; ++k) {
                if (k == i) continue;
                double factor = SwInv[k * 2 * n + i];
                for (int j = 0; j < 2 * n; ++j)
                    SwInv[k * 2 * n + j] -= factor * SwInv[i * 2 * n + j];
            }
        }
        std::vector<double> SwInvMat(n * n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                SwInvMat[i * n + j] = SwInv[i * 2 * n + n + j];

        // Compute M = Sw^{-1} * Sb
        std::vector<double> M(n * n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k)
                    M[i * n + j] += SwInvMat[i * n + k] * Sb[k * n + j];

        // ------------------------------------------------------------------
        // 6. Extract eigenvectors of M using power iteration
        //    We need min(numClasses-1, numBands) discriminant axes.
        // ------------------------------------------------------------------
        int numAxes = std::min(numClasses - 1, n);
        std::vector<std::vector<double>> eigenVecs(numAxes, std::vector<double>(n, 0.0));

        // Power iteration with deflation
        std::vector<double> Mcur(M);
        for (int ax = 0; ax < numAxes; ++ax) {
            std::vector<double> v(n, 1.0);
            // Normalize initial
            double norm = 0.0;
            for (int i = 0; i < n; ++i) norm += v[i] * v[i];
            norm = std::sqrt(norm);
            for (int i = 0; i < n; ++i) v[i] /= norm;

            for (int iter = 0; iter < 200; ++iter) {
                // w = Mcur * v
                std::vector<double> w(n, 0.0);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        w[i] += Mcur[i * n + j] * v[j];

                norm = 0.0;
                for (int i = 0; i < n; ++i) norm += w[i] * w[i];
                norm = std::sqrt(norm);
                if (norm < 1e-15) break;
                for (int i = 0; i < n; ++i) v[i] = w[i] / norm;
            }

            eigenVecs[ax] = v;

            // Deflate: Mcur = Mcur - lambda * v * v^T
            // lambda = v^T * Mcur * v
            std::vector<double> Mv(n, 0.0);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    Mv[i] += Mcur[i * n + j] * v[j];
            double lambda = 0.0;
            for (int i = 0; i < n; ++i)
                lambda += v[i] * Mv[i];

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    Mcur[i * n + j] -= lambda * v[i] * v[j];
        }

        // ------------------------------------------------------------------
        // 7. Project class centroids onto discriminant axes
        // ------------------------------------------------------------------
        std::vector<std::vector<double>> projectedMeans(numClasses, std::vector<double>(numAxes, 0.0));
        for (int c = 0; c < numClasses; ++c)
            for (int ax = 0; ax < numAxes; ++ax)
                for (int b = 0; b < n; ++b)
                    projectedMeans[c][ax] += eigenVecs[ax][b] * classes[c].mean[b];

        // ------------------------------------------------------------------
        // 8. Classify pixels: project each pixel and find nearest class centroid
        // ------------------------------------------------------------------
        reportProgress(0.3, "Classifying pixels...");

        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                out[i] = 0;
                continue;
            }

            // Project pixel onto discriminant axes
            std::vector<double> proj(numAxes);
            for (int ax = 0; ax < numAxes; ++ax) {
                proj[ax] = 0.0;
                for (int b = 0; b < n; ++b)
                    proj[ax] += eigenVecs[ax][b] * (*bands[b])[i];
            }

            // Find nearest class centroid in projected space
            double bestDist = std::numeric_limits<double>::max();
            int bestClass = 0;

            for (int c = 0; c < numClasses; ++c) {
                double dist = 0.0;
                for (int ax = 0; ax < numAxes; ++ax) {
                    double d = proj[ax] - projectedMeans[c][ax];
                    dist += d * d;
                }
                if (dist < bestDist) {
                    bestDist = dist;
                    bestClass = classes[c].classId;
                }
            }

            out[i] = bestClass;

            if (i % 1000000 == 0)
                reportProgress(0.3 + 0.6 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(FisherModule)

} // namespace aplaceholder
