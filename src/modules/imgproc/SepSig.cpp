#include "Module.h"
#include "ModuleRegistry.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <limits>

namespace aplaceholder {

class SepSigModule : public Module {
public:
    QString name() const override { return "SEPSIG"; }
    QString description() const override {
        return "Signature Separability. Reads a signature file and computes pairwise "
               "Jeffries-Matusita separability between all classes. Outputs a "
               "separability matrix as CSV.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("sig_file", "Signature file (CSV)"),
            ParameterDef::output("output_file", "Output separability matrix (CSV)"),
        };
    }

private:
    struct ClassSig {
        int classId;
        std::vector<double> mean;
        std::vector<double> covMatrix;
    };

    bool readSigFile(const QString& path, std::vector<ClassSig>& classes, int& numBands) {
        std::ifstream file(path.toStdString());
        if (!file.is_open()) return false;

        std::string line;
        numBands = -1;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') {
                if (line.find("Number of bands:") != std::string::npos) {
                    std::string token = line.substr(line.find(':') + 1);
                    numBands = std::stoi(token);
                }
                continue;
            }

            std::istringstream iss(line);
            std::vector<double> values;
            std::string token;
            while (std::getline(iss, token, ',')) {
                size_t start = token.find_first_not_of(" \t");
                size_t end = token.find_last_not_of(" \t\r\n");
                if (start != std::string::npos && end != std::string::npos)
                    values.push_back(std::stod(token.substr(start, end - start + 1)));
            }

            if (values.size() < 2) continue;

            if (numBands < 0) {
                int vsize = static_cast<int>(values.size()) - 1;
                for (int n = 1; n < 100; ++n) {
                    if (n + n * n == vsize) { numBands = n; break; }
                }
                if (numBands < 0) return false;
            }

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

            classes.push_back(std::move(sig));
        }

        return !classes.empty() && numBands > 0;
    }

    double determinant(const std::vector<double>& mat, int n) {
        std::vector<double> A(mat);
        double det = 1.0;
        for (int i = 0; i < n; ++i) {
            int maxRow = i;
            for (int k = i + 1; k < n; ++k)
                if (std::fabs(A[k * n + i]) > std::fabs(A[maxRow * n + i]))
                    maxRow = k;
            if (maxRow != i) {
                for (int j = 0; j < n; ++j)
                    std::swap(A[i * n + j], A[maxRow * n + j]);
                det = -det;
            }
            if (std::fabs(A[i * n + i]) < 1e-15) return 0.0;
            det *= A[i * n + i];
            for (int k = i + 1; k < n; ++k) {
                double factor = A[k * n + i] / A[i * n + i];
                for (int j = i; j < n; ++j)
                    A[k * n + j] -= factor * A[i * n + j];
            }
        }
        return det;
    }

    bool invertMatrix(const std::vector<double>& mat, std::vector<double>& inv, int n) {
        std::vector<double> A(n * 2 * n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                A[i * 2 * n + j] = mat[i * n + j];
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
            if (std::fabs(pivot) < 1e-15) return false;
            for (int j = 0; j < 2 * n; ++j)
                A[i * 2 * n + j] /= pivot;
            for (int k = 0; k < n; ++k) {
                if (k == i) continue;
                double factor = A[k * 2 * n + i];
                for (int j = 0; j < 2 * n; ++j)
                    A[k * 2 * n + j] -= factor * A[i * 2 * n + j];
            }
        }
        inv.resize(n * n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                inv[i * n + j] = A[i * 2 * n + n + j];
        return true;
    }

    double jeffriesMatusita(const ClassSig& c1, const ClassSig& c2, int n) {
        std::vector<double> avgCov(n * n);
        for (int i = 0; i < n * n; ++i)
            avgCov[i] = (c1.covMatrix[i] + c2.covMatrix[i]) / 2.0;

        std::vector<double> avgCovInv;
        if (!invertMatrix(avgCov, avgCovInv, n))
            return 2.0;

        std::vector<double> diff(n);
        for (int i = 0; i < n; ++i)
            diff[i] = c1.mean[i] - c2.mean[i];

        double quad = 0.0;
        for (int i = 0; i < n; ++i) {
            double rowSum = 0.0;
            for (int j = 0; j < n; ++j)
                rowSum += avgCovInv[i * n + j] * diff[j];
            quad += diff[i] * rowSum;
        }

        double detAvg = determinant(avgCov, n);
        double det1 = determinant(c1.covMatrix, n);
        double det2 = determinant(c2.covMatrix, n);

        double logTerm = 0.0;
        if (det1 > 0 && det2 > 0 && detAvg > 0)
            logTerm = 0.5 * std::log(detAvg / std::sqrt(det1 * det2));

        double B = quad / 8.0 + logTerm;
        double JM = 2.0 * (1.0 - std::exp(-B));
        return std::min(JM, 2.0);
    }

public:
    bool execute() override {
        std::vector<ClassSig> classes;
        int numBands = 0;

        QString sigPath = parameter("sig_file").toString();
        if (!readSigFile(sigPath, classes, numBands)) {
            setError("Failed to read or parse signature file: " + sigPath);
            return false;
        }

        int numClasses = static_cast<int>(classes.size());
        if (numClasses < 2) {
            setError("At least two classes are required for separability analysis");
            return false;
        }

        reportProgress(0.1, "Computing pairwise J-M distances...");

        // Compute separability matrix
        std::vector<std::vector<double>> sepMatrix(numClasses, std::vector<double>(numClasses, 0.0));

        int totalPairs = numClasses * (numClasses - 1) / 2;
        int pairsDone = 0;

        for (int i = 0; i < numClasses; ++i) {
            for (int j = i + 1; j < numClasses; ++j) {
                double jm = jeffriesMatusita(classes[i], classes[j], numBands);
                sepMatrix[i][j] = jm;
                sepMatrix[j][i] = jm;
                pairsDone++;
                if (pairsDone % 5 == 0)
                    reportProgress(0.1 + 0.7 * static_cast<double>(pairsDone) / totalPairs);
            }
        }

        // Write output CSV
        reportProgress(0.9, "Writing separability matrix...");

        QString outPath = parameter("output_file").toString();
        std::ofstream out(outPath.toStdString());
        if (!out.is_open()) {
            setError("Failed to open output file: " + outPath);
            return false;
        }

        out << "# Signature Separability Matrix (Jeffries-Matusita Distance)\n";
        out << "# Source: " << sigPath.toStdString() << "\n";
        out << "# Number of bands: " << numBands << "\n";
        out << "# Number of classes: " << numClasses << "\n";
        out << "# J-M range: [0, 2]. Values > 1.9 = excellent separability.\n";
        out << "#\n";

        out << std::fixed << std::setprecision(4);

        // Header row
        out << "Class";
        for (int j = 0; j < numClasses; ++j)
            out << ", " << classes[j].classId;
        out << "\n";

        // Data rows
        for (int i = 0; i < numClasses; ++i) {
            out << classes[i].classId;
            for (int j = 0; j < numClasses; ++j)
                out << ", " << sepMatrix[i][j];
            out << "\n";
        }

        // Summary statistics
        double minJM = 2.0, maxJM = 0.0, avgJM = 0.0;
        int pairCount = 0;
        for (int i = 0; i < numClasses; ++i) {
            for (int j = i + 1; j < numClasses; ++j) {
                double v = sepMatrix[i][j];
                if (v < minJM) minJM = v;
                if (v > maxJM) maxJM = v;
                avgJM += v;
                pairCount++;
            }
        }
        if (pairCount > 0) avgJM /= pairCount;

        out << "\n# Summary:\n";
        out << "# Minimum J-M distance: " << minJM << "\n";
        out << "# Maximum J-M distance: " << maxJM << "\n";
        out << "# Average J-M distance: " << avgJM << "\n";

        out.close();

        reportProgress(1.0, QString("Separability matrix written: %1 classes, %2 pairs.")
                       .arg(numClasses).arg(pairCount));
        return true;
    }
};

REGISTER_MODULE(SepSigModule)

} // namespace aplaceholder
