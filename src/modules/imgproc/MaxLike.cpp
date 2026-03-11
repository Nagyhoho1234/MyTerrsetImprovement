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

class MaxLikeModule : public Module {
public:
    QString name() const override { return "MAXLIKE"; }
    QString description() const override {
        return "Maximum Likelihood Classification. Assigns each pixel to the class "
               "with the highest probability based on training site signatures.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input multi-band image"),
            ParameterDef::file("signature_file", "Signature file (CSV)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::combo("prior_probabilities", "Prior probabilities",
                {"Equal", "From signature file"}, 0,
                "Method for assigning prior probabilities to classes"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input image");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int numBands = raster->bands();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        // Parse signature file
        // Expected CSV format:
        //   class_id, mean_b1, mean_b2, ..., mean_bN, cov_11, cov_12, ..., cov_NN
        // One row per class. Optionally a prior probability at the end if
        // "From signature file" is selected.
        QString sigPath = parameter("signature_file").toString();
        int priorMode = parameter("prior_probabilities").toInt();

        struct ClassSig {
            int classId;
            std::vector<double> mean;       // length = numBands
            std::vector<double> covMatrix;  // numBands x numBands, row-major
            double prior;
            // Precomputed for classification
            std::vector<double> covInverse; // numBands x numBands
            double logDetCov;
        };

        std::vector<ClassSig> classes;

        // Read signature file
        {
            std::ifstream file(sigPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open signature file: " + sigPath);
                return false;
            }

            std::string line;
            while (std::getline(file, line)) {
                // Skip empty lines and comments
                if (line.empty() || line[0] == '#') continue;

                std::istringstream iss(line);
                std::vector<double> values;
                std::string token;
                while (std::getline(iss, token, ',')) {
                    // Trim whitespace
                    size_t start = token.find_first_not_of(" \t");
                    size_t end = token.find_last_not_of(" \t\r\n");
                    if (start != std::string::npos && end != std::string::npos)
                        values.push_back(std::stod(token.substr(start, end - start + 1)));
                }

                // Expected: 1 (class_id) + numBands (means) + numBands*numBands (covariance)
                // Optionally + 1 (prior)
                int expectedMin = 1 + numBands + numBands * numBands;
                if (static_cast<int>(values.size()) < expectedMin) continue;

                ClassSig sig;
                sig.classId = static_cast<int>(values[0]);
                sig.mean.resize(numBands);
                for (int b = 0; b < numBands; ++b)
                    sig.mean[b] = values[1 + b];

                sig.covMatrix.resize(numBands * numBands);
                for (int i = 0; i < numBands * numBands; ++i)
                    sig.covMatrix[i] = values[1 + numBands + i];

                // Prior probability (if provided and requested)
                if (priorMode == 1 && static_cast<int>(values.size()) > expectedMin)
                    sig.prior = values[expectedMin];
                else
                    sig.prior = 1.0; // will be normalized later

                classes.push_back(std::move(sig));
            }
        }

        if (classes.empty()) {
            setError("No valid class signatures found in signature file");
            return false;
        }

        // Normalize priors
        if (priorMode == 0) {
            // Equal priors
            double p = 1.0 / classes.size();
            for (auto& cls : classes)
                cls.prior = p;
        } else {
            // Normalize priors from file
            double sumPrior = 0.0;
            for (const auto& cls : classes)
                sumPrior += cls.prior;
            if (sumPrior > 0) {
                for (auto& cls : classes)
                    cls.prior /= sumPrior;
            }
        }

        reportProgress(0.1, "Computing covariance inverses (Cholesky)...");

        // Precompute inverse covariance matrices and log-determinants
        // using Cholesky decomposition: Sigma = L * L^T
        int n = numBands;
        for (auto& cls : classes) {
            // Cholesky decomposition: compute lower triangular L
            std::vector<double> L(n * n, 0.0);
            bool choleskyOk = true;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j <= i; ++j) {
                    double sum = cls.covMatrix[i * n + j];
                    for (int k = 0; k < j; ++k)
                        sum -= L[i * n + k] * L[j * n + k];
                    if (i == j) {
                        if (sum <= 0) {
                            // Matrix not positive definite; add small regularization
                            sum = 1e-6;
                            choleskyOk = false;
                        }
                        L[i * n + j] = std::sqrt(sum);
                    } else {
                        L[i * n + j] = sum / L[j * n + j];
                    }
                }
            }

            // Log determinant = 2 * sum(log(diag(L)))
            cls.logDetCov = 0.0;
            for (int i = 0; i < n; ++i)
                cls.logDetCov += 2.0 * std::log(L[i * n + i]);

            // Inverse of L (lower triangular inverse)
            std::vector<double> Linv(n * n, 0.0);
            for (int i = 0; i < n; ++i) {
                Linv[i * n + i] = 1.0 / L[i * n + i];
                for (int j = i + 1; j < n; ++j) {
                    double sum = 0.0;
                    for (int k = i; k < j; ++k)
                        sum -= L[j * n + k] * Linv[k * n + i];
                    Linv[j * n + i] = sum / L[j * n + j];
                }
            }

            // Sigma^{-1} = (L^T)^{-1} * L^{-1} = Linv^T * Linv
            cls.covInverse.resize(n * n, 0.0);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    double sum = 0.0;
                    // covInv[i][j] = sum_k Linv[k][i] * Linv[k][j]
                    // (Linv^T * Linv), but Linv is lower triangular
                    // so Linv[k][i] is nonzero only for k >= i
                    // and Linv[k][j] is nonzero only for k >= j
                    int kStart = std::max(i, j);
                    for (int k = kStart; k < n; ++k)
                        sum += Linv[k * n + i] * Linv[k * n + j];
                    cls.covInverse[i * n + j] = sum;
                }
        }

        reportProgress(0.2, "Classifying pixels...");

        // Collect band data pointers
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &raster->data(b);

        // Create output classified image
        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);

        // For each pixel, compute log-posterior for each class and assign to max
        // log P(x|c) = -0.5 * n * log(2*pi) - 0.5 * log|Sigma_c|
        //              - 0.5 * (x - mu_c)^T * Sigma_c^{-1} * (x - mu_c)
        // log posterior = log P(x|c) + log P(c)
        double logTwoPiTerm = -0.5 * n * std::log(2.0 * M_PI);

        std::vector<double> diff(n);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                out[i] = 0; // unclassified
                continue;
            }

            double bestLogPost = -std::numeric_limits<double>::max();
            int bestClass = 0;

            for (const auto& cls : classes) {
                // Compute (x - mu)
                for (int b = 0; b < n; ++b)
                    diff[b] = (*bands[b])[i] - cls.mean[b];

                // Compute (x - mu)^T * Sigma^{-1} * (x - mu)
                double mahal = 0.0;
                for (int b1 = 0; b1 < n; ++b1) {
                    double rowSum = 0.0;
                    for (int b2 = 0; b2 < n; ++b2)
                        rowSum += cls.covInverse[b1 * n + b2] * diff[b2];
                    mahal += diff[b1] * rowSum;
                }

                double logPost = logTwoPiTerm
                                 - 0.5 * cls.logDetCov
                                 - 0.5 * mahal
                                 + std::log(cls.prior);

                if (logPost > bestLogPost) {
                    bestLogPost = logPost;
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

REGISTER_MODULE(MaxLikeModule)

} // namespace aplaceholder
