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

class BayClassModule : public Module {
public:
    QString name() const override { return "BAYCLASS"; }
    QString description() const override {
        return "Bayesian soft classifier with prior probabilities. Computes posterior "
               "probability for each class per pixel using Bayes' theorem with "
               "multivariate normal likelihoods. Outputs one probability raster per class "
               "plus a classification uncertainty image.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated paths)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::file("signature_file", "Signature file (CSV)"),
            ParameterDef::file("prior_probs", "Prior probabilities (comma-separated, optional)",
                "Comma-separated prior probabilities per class, in signature file order. "
                "Leave empty for equal priors."),
            ParameterDef::output("output", "Output prefix for probability images"),
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
        int n = numBands;

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
            std::vector<double> covInverse;
            double logDetCov;
            double prior;
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
                sig.prior = 1.0; // will be set later

                classes.push_back(std::move(sig));
            }
        }

        if (classes.empty()) {
            setError("No valid class signatures found");
            return false;
        }

        int numClasses = static_cast<int>(classes.size());

        // ------------------------------------------------------------------
        // 4. Parse prior probabilities
        // ------------------------------------------------------------------
        QString priorParam = parameter("prior_probs").toString().trimmed();
        if (!priorParam.isEmpty()) {
            QStringList priorStrs = priorParam.split(",", Qt::SkipEmptyParts);
            if (priorStrs.size() == numClasses) {
                double sumP = 0.0;
                for (int c = 0; c < numClasses; ++c) {
                    classes[c].prior = priorStrs[c].trimmed().toDouble();
                    sumP += classes[c].prior;
                }
                if (sumP > 0)
                    for (auto& cls : classes)
                        cls.prior /= sumP;
            } else {
                // Mismatch: use equal priors
                for (auto& cls : classes)
                    cls.prior = 1.0 / numClasses;
            }
        } else {
            for (auto& cls : classes)
                cls.prior = 1.0 / numClasses;
        }

        // ------------------------------------------------------------------
        // 5. Precompute inverse covariance and log-determinant (Cholesky)
        // ------------------------------------------------------------------
        reportProgress(0.1, "Computing covariance inverses...");

        for (auto& cls : classes) {
            std::vector<double> L(n * n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j <= i; ++j) {
                    double sum = cls.covMatrix[i * n + j];
                    for (int k = 0; k < j; ++k)
                        sum -= L[i * n + k] * L[j * n + k];
                    if (i == j) {
                        if (sum <= 0) sum = 1e-6;
                        L[i * n + j] = std::sqrt(sum);
                    } else {
                        L[i * n + j] = sum / L[j * n + j];
                    }
                }
            }

            cls.logDetCov = 0.0;
            for (int i = 0; i < n; ++i)
                cls.logDetCov += 2.0 * std::log(L[i * n + i]);

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

            cls.covInverse.resize(n * n, 0.0);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    double sum = 0.0;
                    int kStart = std::max(i, j);
                    for (int k = kStart; k < n; ++k)
                        sum += Linv[k * n + i] * Linv[k * n + j];
                    cls.covInverse[i * n + j] = sum;
                }
        }

        // ------------------------------------------------------------------
        // 6. Create output rasters: one probability image per class + uncertainty
        // ------------------------------------------------------------------
        reportProgress(0.2, "Classifying pixels...");

        QString prefix = parameter("output").toString();

        std::vector<Raster> probOutputs;
        probOutputs.reserve(numClasses);
        for (int c = 0; c < numClasses; ++c) {
            probOutputs.emplace_back(cols, rows, 1, DataType::Float32);
            probOutputs[c].setGeoTransform(bandRasters[0]->geoTransform());
            probOutputs[c].setProjection(bandRasters[0]->projection());
            probOutputs[c].setNoDataValue(-1.0);
        }

        Raster uncertOutput(cols, rows, 1, DataType::Float32);
        uncertOutput.setGeoTransform(bandRasters[0]->geoTransform());
        uncertOutput.setProjection(bandRasters[0]->projection());
        uncertOutput.setNoDataValue(-1.0);

        std::vector<std::vector<double>*> probData(numClasses);
        for (int c = 0; c < numClasses; ++c)
            probData[c] = &probOutputs[c].data(0);
        auto& uncertData = uncertOutput.data(0);

        double logTwoPiTerm = -0.5 * n * std::log(2.0 * M_PI);
        std::vector<double> diff(n);
        std::vector<double> logLikelihoods(numClasses);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                for (int c = 0; c < numClasses; ++c)
                    (*probData[c])[i] = -1.0;
                uncertData[i] = -1.0;
                continue;
            }

            // Compute log-likelihood for each class
            double maxLogPost = -std::numeric_limits<double>::max();
            for (int c = 0; c < numClasses; ++c) {
                const auto& cls = classes[c];
                for (int b = 0; b < n; ++b)
                    diff[b] = (*bands[b])[i] - cls.mean[b];

                double mahal = 0.0;
                for (int b1 = 0; b1 < n; ++b1) {
                    double rowSum = 0.0;
                    for (int b2 = 0; b2 < n; ++b2)
                        rowSum += cls.covInverse[b1 * n + b2] * diff[b2];
                    mahal += diff[b1] * rowSum;
                }

                logLikelihoods[c] = logTwoPiTerm - 0.5 * cls.logDetCov
                                    - 0.5 * mahal + std::log(cls.prior);
                if (logLikelihoods[c] > maxLogPost)
                    maxLogPost = logLikelihoods[c];
            }

            // Convert to posterior probabilities using log-sum-exp trick
            double sumExp = 0.0;
            for (int c = 0; c < numClasses; ++c)
                sumExp += std::exp(logLikelihoods[c] - maxLogPost);
            double logSumExp = maxLogPost + std::log(sumExp);

            double maxProb = 0.0;
            double sumProb = 0.0;
            for (int c = 0; c < numClasses; ++c) {
                double post = std::exp(logLikelihoods[c] - logSumExp);
                (*probData[c])[i] = post;
                if (post > maxProb) maxProb = post;
                sumProb += post;
            }

            // Classification uncertainty = 1 - (max - sum/n) / (1 - 1/n)
            double nC = static_cast<double>(numClasses);
            double uncert = 1.0;
            if (nC > 1.0)
                uncert = 1.0 - (maxProb - sumProb / nC) / (1.0 - 1.0 / nC);
            uncertData[i] = std::max(0.0, std::min(1.0, uncert));

            if (i % 1000000 == 0)
                reportProgress(0.2 + 0.7 * static_cast<double>(i) / total);
        }

        // ------------------------------------------------------------------
        // 7. Write outputs
        // ------------------------------------------------------------------
        reportProgress(0.95, "Writing output images...");

        for (int c = 0; c < numClasses; ++c) {
            QString outPath = QString("%1_class%2.tif").arg(prefix).arg(classes[c].classId);
            if (!GdalIO::write(probOutputs[c], outPath)) {
                setError("Failed to write probability image: " + outPath);
                return false;
            }
        }

        QString uncertPath = prefix + "_uncertainty.tif";
        if (!GdalIO::write(uncertOutput, uncertPath)) {
            setError("Failed to write uncertainty image: " + uncertPath);
            return false;
        }

        reportProgress(1.0, QString("Bayesian classification complete: %1 class probability images + uncertainty.")
                       .arg(numClasses));
        return true;
    }
};

REGISTER_MODULE(BayClassModule)

} // namespace aplaceholder
