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

class BelClassModule : public Module {
public:
    QString name() const override { return "BELCLASS"; }
    QString description() const override {
        return "Dempster-Shafer Belief classifier. Computes belief, plausibility, and "
               "belief interval per class for each pixel. Unlike Bayesian classification, "
               "explicitly recognizes ignorance (unknown classes). Outputs belief and "
               "plausibility images per class plus a Dempster-Shafer ignorance image.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated paths)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::file("signature_file", "Signature file (CSV)"),
            ParameterDef::output("output_prefix", "Output prefix for belief/plausibility images"),
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

        int numClasses = static_cast<int>(classes.size());

        // ------------------------------------------------------------------
        // 4. Precompute inverse covariance and log-determinant (Cholesky)
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
        // 5. Create output rasters
        // ------------------------------------------------------------------
        reportProgress(0.2, "Computing belief and plausibility...");

        QString prefix = parameter("output_prefix").toString();

        std::vector<Raster> beliefOutputs, plausOutputs;
        beliefOutputs.reserve(numClasses);
        plausOutputs.reserve(numClasses);
        for (int c = 0; c < numClasses; ++c) {
            beliefOutputs.emplace_back(cols, rows, 1, DataType::Float32);
            beliefOutputs[c].setGeoTransform(bandRasters[0]->geoTransform());
            beliefOutputs[c].setProjection(bandRasters[0]->projection());
            beliefOutputs[c].setNoDataValue(-1.0);

            plausOutputs.emplace_back(cols, rows, 1, DataType::Float32);
            plausOutputs[c].setGeoTransform(bandRasters[0]->geoTransform());
            plausOutputs[c].setProjection(bandRasters[0]->projection());
            plausOutputs[c].setNoDataValue(-1.0);
        }

        Raster ignoranceOutput(cols, rows, 1, DataType::Float32);
        ignoranceOutput.setGeoTransform(bandRasters[0]->geoTransform());
        ignoranceOutput.setProjection(bandRasters[0]->projection());
        ignoranceOutput.setNoDataValue(-1.0);

        std::vector<std::vector<double>*> beliefData(numClasses), plausData(numClasses);
        for (int c = 0; c < numClasses; ++c) {
            beliefData[c] = &beliefOutputs[c].data(0);
            plausData[c] = &plausOutputs[c].data(0);
        }
        auto& ignoranceData = ignoranceOutput.data(0);

        double logTwoPiTerm = -0.5 * n * std::log(2.0 * M_PI);
        std::vector<double> diff(n);
        std::vector<double> likelihoods(numClasses);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                for (int c = 0; c < numClasses; ++c) {
                    (*beliefData[c])[i] = -1.0;
                    (*plausData[c])[i] = -1.0;
                }
                ignoranceData[i] = -1.0;
                continue;
            }

            // Compute likelihoods p(e|h_c) for each class using MVN
            double maxLogLike = -std::numeric_limits<double>::max();
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

                likelihoods[c] = logTwoPiTerm - 0.5 * cls.logDetCov - 0.5 * mahal;
                if (likelihoods[c] > maxLogLike)
                    maxLogLike = likelihoods[c];
            }

            // Convert log-likelihoods to likelihoods (using log-sum-exp for stability)
            double sumLike = 0.0;
            for (int c = 0; c < numClasses; ++c) {
                likelihoods[c] = std::exp(likelihoods[c] - maxLogLike);
                sumLike += likelihoods[c];
            }

            // Normalize likelihoods to get basic probability assignments (BPA)
            // In D-S theory with singleton evidence, BPA for each class = likelihood
            // The remaining mass goes to the full frame (ignorance)
            // We model: m(h_c) = proportional to likelihood, m(Theta) = ignorance
            // Scale so that sum of singleton BPAs + ignorance = 1
            // Use a scaling factor: total evidence = sum of likelihoods
            // If total evidence is high, ignorance is low; if low, ignorance is high

            // Compute average likelihood relative to a uniform reference
            // Ignorance = proportion of mass not committed to any singleton
            double maxLike = 0.0;
            for (int c = 0; c < numClasses; ++c) {
                likelihoods[c] /= sumLike; // normalize to sum to 1
                if (likelihoods[c] > maxLike) maxLike = likelihoods[c];
            }

            // Dempster-Shafer: commitment factor based on how peaked the
            // distribution is. If all equal => maximum ignorance.
            // commitment = 1 - entropy_ratio
            double entropy = 0.0;
            double maxEntropy = std::log(static_cast<double>(numClasses));
            for (int c = 0; c < numClasses; ++c) {
                if (likelihoods[c] > 1e-15)
                    entropy -= likelihoods[c] * std::log(likelihoods[c]);
            }
            double commitment = (maxEntropy > 0) ? (1.0 - entropy / maxEntropy) : 0.0;
            commitment = std::max(0.0, std::min(1.0, commitment));

            // Belief(h_c) = commitment * likelihood(h_c)
            // Plausibility(h_c) = Belief(h_c) + ignorance = Belief(h_c) + (1 - commitment)
            // Belief interval = Plausibility - Belief = ignorance = (1 - commitment)
            double ignorance = 1.0 - commitment;

            for (int c = 0; c < numClasses; ++c) {
                double bel = commitment * likelihoods[c];
                double plaus = bel + ignorance;
                (*beliefData[c])[i] = bel;
                (*plausData[c])[i] = std::min(plaus, 1.0);
            }

            ignoranceData[i] = ignorance;

            if (i % 1000000 == 0)
                reportProgress(0.2 + 0.7 * static_cast<double>(i) / total);
        }

        // ------------------------------------------------------------------
        // 6. Write outputs
        // ------------------------------------------------------------------
        reportProgress(0.95, "Writing output images...");

        for (int c = 0; c < numClasses; ++c) {
            QString belPath = QString("%1_belief_class%2.tif").arg(prefix).arg(classes[c].classId);
            if (!GdalIO::write(beliefOutputs[c], belPath)) {
                setError("Failed to write belief image: " + belPath);
                return false;
            }

            QString plausPath = QString("%1_plaus_class%2.tif").arg(prefix).arg(classes[c].classId);
            if (!GdalIO::write(plausOutputs[c], plausPath)) {
                setError("Failed to write plausibility image: " + plausPath);
                return false;
            }
        }

        QString ignPath = prefix + "_ignorance.tif";
        if (!GdalIO::write(ignoranceOutput, ignPath)) {
            setError("Failed to write ignorance image: " + ignPath);
            return false;
        }

        reportProgress(1.0, QString("Dempster-Shafer classification complete: %1 classes, "
                       "belief + plausibility + ignorance images.").arg(numClasses));
        return true;
    }
};

REGISTER_MODULE(BelClassModule)

} // namespace aplaceholder
