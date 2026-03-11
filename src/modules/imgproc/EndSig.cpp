#include "Module.h"
#include "ModuleRegistry.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class EndSigModule : public Module {
public:
    QString name() const override { return "ENDSIG"; }
    QString description() const override {
        return "Evaluate signature effectiveness. Computes spectral contrast metrics between "
               "endmembers including average spectral angle, minimum spectral angle, and "
               "spectral information divergence (SID).";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("endmember_file", "Endmember file (CSV)",
                "CSV: each row is an endmember, columns are band values"),
            ParameterDef::output("output_report", "Output report file",
                "Text report with spectral contrast metrics"),
        };
    }

    bool execute() override {
        // --------------------------------------------------------------------
        // 1. Read endmember file (CSV)
        // --------------------------------------------------------------------
        QString emPath = parameter("endmember_file").toString();
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
                    if (start != std::string::npos && end != std::string::npos) {
                        try {
                            values.push_back(std::stod(token.substr(start, end - start + 1)));
                        } catch (...) {}
                    }
                }

                if (!values.empty())
                    endmembers.push_back(std::move(values));
            }
        }

        int numEM = static_cast<int>(endmembers.size());
        if (numEM < 2) {
            setError("At least 2 endmembers required to evaluate spectral contrast");
            return false;
        }

        // Determine number of bands (minimum across all endmembers)
        int numBands = static_cast<int>(endmembers[0].size());
        for (int e = 1; e < numEM; ++e) {
            int sz = static_cast<int>(endmembers[e].size());
            if (sz < numBands) numBands = sz;
        }
        for (auto& em : endmembers) em.resize(numBands);

        if (numBands == 0) {
            setError("Endmembers have no spectral bands");
            return false;
        }

        reportProgress(0.1, "Computing spectral contrast metrics...");

        // --------------------------------------------------------------------
        // 2. Compute pairwise spectral angles
        // --------------------------------------------------------------------
        int numPairs = numEM * (numEM - 1) / 2;
        std::vector<double> angles;
        angles.reserve(numPairs);

        // Angle matrix (for report)
        std::vector<std::vector<double>> angleMatrix(numEM, std::vector<double>(numEM, 0.0));

        for (int i = 0; i < numEM; ++i) {
            for (int j = i + 1; j < numEM; ++j) {
                double dotProduct = 0.0;
                double magI = 0.0, magJ = 0.0;

                for (int b = 0; b < numBands; ++b) {
                    dotProduct += endmembers[i][b] * endmembers[j][b];
                    magI += endmembers[i][b] * endmembers[i][b];
                    magJ += endmembers[j][b] * endmembers[j][b];
                }

                magI = std::sqrt(magI);
                magJ = std::sqrt(magJ);

                double angle = 0.0;
                if (magI > 1e-12 && magJ > 1e-12) {
                    double cosAngle = dotProduct / (magI * magJ);
                    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                    angle = std::acos(cosAngle);
                }

                angleMatrix[i][j] = angle;
                angleMatrix[j][i] = angle;
                angles.push_back(angle);
            }
        }

        double avgAngle = 0.0;
        for (double a : angles) avgAngle += a;
        avgAngle /= angles.size();

        double minAngle = *std::min_element(angles.begin(), angles.end());
        double maxAngle = *std::max_element(angles.begin(), angles.end());

        reportProgress(0.5, "Computing spectral information divergence...");

        // --------------------------------------------------------------------
        // 3. Compute pairwise Spectral Information Divergence (SID)
        // --------------------------------------------------------------------
        // SID(x,y) = D(x||y) + D(y||x) where D is KL divergence
        // First convert each endmember to a probability distribution (normalize to sum=1)
        std::vector<std::vector<double>> emProb(numEM, std::vector<double>(numBands));
        for (int e = 0; e < numEM; ++e) {
            double sum = 0.0;
            for (int b = 0; b < numBands; ++b) {
                // Use absolute values for probability (spectra should be positive)
                emProb[e][b] = std::max(endmembers[e][b], 1e-10);
                sum += emProb[e][b];
            }
            if (sum > 1e-12) {
                for (int b = 0; b < numBands; ++b)
                    emProb[e][b] /= sum;
            }
        }

        std::vector<std::vector<double>> sidMatrix(numEM, std::vector<double>(numEM, 0.0));
        std::vector<double> sidValues;
        sidValues.reserve(numPairs);

        for (int i = 0; i < numEM; ++i) {
            for (int j = i + 1; j < numEM; ++j) {
                double dKL_ij = 0.0, dKL_ji = 0.0;

                for (int b = 0; b < numBands; ++b) {
                    double pi = emProb[i][b];
                    double pj = emProb[j][b];

                    if (pi > 1e-15 && pj > 1e-15) {
                        dKL_ij += pi * std::log(pi / pj);
                        dKL_ji += pj * std::log(pj / pi);
                    }
                }

                double sid = dKL_ij + dKL_ji;
                sidMatrix[i][j] = sid;
                sidMatrix[j][i] = sid;
                sidValues.push_back(sid);
            }
        }

        double avgSID = 0.0;
        for (double s : sidValues) avgSID += s;
        avgSID /= sidValues.size();

        double minSID = *std::min_element(sidValues.begin(), sidValues.end());

        reportProgress(0.8, "Writing report...");

        // --------------------------------------------------------------------
        // 4. Write output report
        // --------------------------------------------------------------------
        QString reportPath = parameter("output_report").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report for writing: " + reportPath);
            return false;
        }

        outFile << "===================================================\n";
        outFile << "  ENDSIG - Endmember Signature Effectiveness Report\n";
        outFile << "===================================================\n\n";

        outFile << "Number of endmembers: " << numEM << "\n";
        outFile << "Number of bands: " << numBands << "\n\n";

        outFile << "--- Spectral Angle Summary ---\n";
        outFile << "Average spectral angle: " << avgAngle << " radians ("
                << (avgAngle * 180.0 / M_PI) << " degrees)\n";
        outFile << "Minimum spectral angle: " << minAngle << " radians ("
                << (minAngle * 180.0 / M_PI) << " degrees)\n";
        outFile << "Maximum spectral angle: " << maxAngle << " radians ("
                << (maxAngle * 180.0 / M_PI) << " degrees)\n\n";

        outFile << "--- Spectral Angle Matrix (radians) ---\n";
        outFile << "EM";
        for (int j = 0; j < numEM; ++j)
            outFile << "\t" << (j + 1);
        outFile << "\n";
        for (int i = 0; i < numEM; ++i) {
            outFile << (i + 1);
            for (int j = 0; j < numEM; ++j)
                outFile << "\t" << angleMatrix[i][j];
            outFile << "\n";
        }
        outFile << "\n";

        outFile << "--- Spectral Information Divergence Summary ---\n";
        outFile << "Average SID: " << avgSID << "\n";
        outFile << "Minimum SID: " << minSID << "\n\n";

        outFile << "--- SID Matrix ---\n";
        outFile << "EM";
        for (int j = 0; j < numEM; ++j)
            outFile << "\t" << (j + 1);
        outFile << "\n";
        for (int i = 0; i < numEM; ++i) {
            outFile << (i + 1);
            for (int j = 0; j < numEM; ++j)
                outFile << "\t" << sidMatrix[i][j];
            outFile << "\n";
        }
        outFile << "\n";

        // Separability assessment
        outFile << "--- Separability Assessment ---\n";
        if (minAngle < 0.05)
            outFile << "WARNING: Some endmembers have very small spectral angles (<0.05 rad). "
                       "These may be difficult to distinguish.\n";
        else if (minAngle < 0.1)
            outFile << "CAUTION: Some endmembers have moderate spectral angles (0.05-0.1 rad). "
                       "Classification may be challenging.\n";
        else
            outFile << "GOOD: All endmember pairs have sufficient spectral separation "
                       "(min angle > 0.1 rad).\n";

        outFile.close();

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(EndSigModule)

} // namespace aplaceholder
