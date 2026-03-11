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

class FuzClassModule : public Module {
public:
    QString name() const override { return "FUZCLASS"; }
    QString description() const override {
        return "Fuzzy Set classifier. Computes fuzzy membership values per class based "
               "on standardized distance from class mean. Membership is 1.0 at the mean "
               "and linearly decreases to 0.0 at a user-specified z-score distance. "
               "Outputs one membership raster per class plus an uncertainty image.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated paths)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::file("fuzzy_sig_file", "Signature file (CSV)",
                "CSV signature file with mean and covariance per class"),
            ParameterDef::output("output_prefix", "Output prefix for membership images"),
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
        QString sigPath = parameter("fuzzy_sig_file").toString();

        struct ClassSig {
            int classId;
            std::vector<double> mean;
            std::vector<double> stddev; // per-band standard deviations
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

                // Need at least class_id + means + covariance matrix
                int expectedMin = 1 + n + n * n;
                if (static_cast<int>(values.size()) < expectedMin) {
                    // Try simpler format: class_id + means + stddevs
                    int simpleMin = 1 + n + n;
                    if (static_cast<int>(values.size()) >= simpleMin) {
                        ClassSig sig;
                        sig.classId = static_cast<int>(values[0]);
                        sig.mean.resize(n);
                        sig.stddev.resize(n);
                        for (int b = 0; b < n; ++b)
                            sig.mean[b] = values[1 + b];
                        for (int b = 0; b < n; ++b) {
                            double sd = values[1 + n + b];
                            sig.stddev[b] = (sd > 1e-12) ? sd : 1e-12;
                        }
                        classes.push_back(std::move(sig));
                    }
                    continue;
                }

                ClassSig sig;
                sig.classId = static_cast<int>(values[0]);
                sig.mean.resize(n);
                for (int b = 0; b < n; ++b)
                    sig.mean[b] = values[1 + b];

                // Extract per-band standard deviations from covariance matrix diagonal
                sig.stddev.resize(n);
                for (int b = 0; b < n; ++b) {
                    double var = values[1 + n + b * n + b]; // diagonal element
                    sig.stddev[b] = (var > 1e-12) ? std::sqrt(var) : 1e-6;
                }

                classes.push_back(std::move(sig));
            }
        }

        if (classes.empty()) {
            setError("No valid class signatures found");
            return false;
        }

        int numClasses = static_cast<int>(classes.size());

        // z-score threshold where membership becomes zero
        // Default: 1.96 (2-sigma, 95% coverage)
        double zThreshold = 1.96;

        // ------------------------------------------------------------------
        // 4. Classify pixels with fuzzy membership
        // ------------------------------------------------------------------
        reportProgress(0.1, "Computing fuzzy memberships...");

        QString prefix = parameter("output_prefix").toString();

        std::vector<Raster> memberOutputs;
        memberOutputs.reserve(numClasses);
        for (int c = 0; c < numClasses; ++c) {
            memberOutputs.emplace_back(cols, rows, 1, DataType::Float32);
            memberOutputs[c].setGeoTransform(bandRasters[0]->geoTransform());
            memberOutputs[c].setProjection(bandRasters[0]->projection());
            memberOutputs[c].setNoDataValue(-1.0);
        }

        Raster uncertOutput(cols, rows, 1, DataType::Float32);
        uncertOutput.setGeoTransform(bandRasters[0]->geoTransform());
        uncertOutput.setProjection(bandRasters[0]->projection());
        uncertOutput.setNoDataValue(-1.0);

        std::vector<std::vector<double>*> memberData(numClasses);
        for (int c = 0; c < numClasses; ++c)
            memberData[c] = &memberOutputs[c].data(0);
        auto& uncertData = uncertOutput.data(0);

        std::vector<double> memberships(numClasses);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                for (int c = 0; c < numClasses; ++c)
                    (*memberData[c])[i] = -1.0;
                uncertData[i] = -1.0;
                continue;
            }

            // Compute standardized distance for each class
            double sumMembership = 0.0;
            double maxMembership = 0.0;

            for (int c = 0; c < numClasses; ++c) {
                const auto& cls = classes[c];

                // Standardized Euclidean distance
                double distSq = 0.0;
                for (int b = 0; b < n; ++b) {
                    double d = ((*bands[b])[i] - cls.mean[b]) / cls.stddev[b];
                    distSq += d * d;
                }
                double zDist = std::sqrt(distSq / n); // average z-score

                // Linear fuzzy membership: 1.0 at distance 0, 0.0 at zThreshold
                double membership = 1.0 - zDist / zThreshold;
                if (membership < 0.0) membership = 0.0;
                if (membership > 1.0) membership = 1.0;

                memberships[c] = membership;
                sumMembership += membership;
                if (membership > maxMembership)
                    maxMembership = membership;
            }

            // Normalize memberships so they sum to 1.0
            if (sumMembership > 0) {
                for (int c = 0; c < numClasses; ++c) {
                    memberships[c] /= sumMembership;
                    (*memberData[c])[i] = memberships[c];
                }
            } else {
                for (int c = 0; c < numClasses; ++c)
                    (*memberData[c])[i] = 0.0;
            }

            // Classification uncertainty
            double nC = static_cast<double>(numClasses);
            double normMax = (sumMembership > 0) ? maxMembership / sumMembership : 0.0;
            double uncert = 1.0;
            if (nC > 1.0 && sumMembership > 0)
                uncert = 1.0 - (normMax - 1.0 / nC) / (1.0 - 1.0 / nC);
            uncertData[i] = std::max(0.0, std::min(1.0, uncert));

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        // ------------------------------------------------------------------
        // 5. Write outputs
        // ------------------------------------------------------------------
        reportProgress(0.95, "Writing output images...");

        for (int c = 0; c < numClasses; ++c) {
            QString outPath = QString("%1_class%2.tif").arg(prefix).arg(classes[c].classId);
            if (!GdalIO::write(memberOutputs[c], outPath)) {
                setError("Failed to write membership image: " + outPath);
                return false;
            }
        }

        QString uncertPath = prefix + "_uncertainty.tif";
        if (!GdalIO::write(uncertOutput, uncertPath)) {
            setError("Failed to write uncertainty image: " + uncertPath);
            return false;
        }

        reportProgress(1.0, QString("Fuzzy classification complete: %1 membership images + uncertainty.")
                       .arg(numClasses));
        return true;
    }
};

REGISTER_MODULE(FuzClassModule)

} // namespace aplaceholder
