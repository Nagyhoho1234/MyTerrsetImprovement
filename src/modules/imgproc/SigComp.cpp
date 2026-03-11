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

class SigCompModule : public Module {
public:
    QString name() const override { return "SIGCOMP"; }
    QString description() const override {
        return "Signature Comparison. Reads two signature files and computes "
               "Jeffries-Matusita distance and Transformed Divergence between "
               "all class pairs. Outputs a comparison report.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("sig_file1", "First signature file (CSV)"),
            ParameterDef::file("sig_file2", "Second signature file (CSV)"),
            ParameterDef::output("output_report", "Output comparison report"),
        };
    }

private:
    struct ClassSig {
        int classId;
        std::vector<double> mean;
        std::vector<double> covMatrix; // numBands x numBands, row-major
    };

    bool readSigFile(const QString& path, std::vector<ClassSig>& classes, int& numBands) {
        std::ifstream file(path.toStdString());
        if (!file.is_open()) return false;

        std::string line;
        numBands = -1;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') {
                // Try to parse number of bands from comment
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

            // If numBands not from comment, infer from line length:
            // values = 1 + N + N*N => solve for N
            if (numBands < 0) {
                int vsize = static_cast<int>(values.size()) - 1;
                // N + N*N = vsize => N*(1+N) = vsize
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

    // Compute determinant of n x n matrix using LU decomposition
    double determinant(const std::vector<double>& mat, int n) {
        std::vector<double> A(mat);
        double det = 1.0;
        for (int i = 0; i < n; ++i) {
            // Partial pivoting
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

    // Compute inverse of n x n matrix
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

    // Compute trace of n x n matrix
    double trace(const std::vector<double>& mat, int n) {
        double t = 0.0;
        for (int i = 0; i < n; ++i)
            t += mat[i * n + i];
        return t;
    }

    // Matrix multiply: C = A * B (all n x n)
    std::vector<double> matMul(const std::vector<double>& A, const std::vector<double>& B, int n) {
        std::vector<double> C(n * n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k)
                    C[i * n + j] += A[i * n + k] * B[k * n + j];
        return C;
    }

    // Jeffries-Matusita distance between two multivariate normal distributions
    // JM = 2(1 - e^{-B}) where B is the Bhattacharyya distance
    // B = (1/8)(mu1-mu2)^T * [(S1+S2)/2]^{-1} * (mu1-mu2)
    //   + (1/2) * ln( det((S1+S2)/2) / sqrt(det(S1)*det(S2)) )
    double jeffriesMatusita(const ClassSig& c1, const ClassSig& c2, int n) {
        // Average covariance
        std::vector<double> avgCov(n * n);
        for (int i = 0; i < n * n; ++i)
            avgCov[i] = (c1.covMatrix[i] + c2.covMatrix[i]) / 2.0;

        std::vector<double> avgCovInv;
        if (!invertMatrix(avgCov, avgCovInv, n))
            return 2.0; // maximum distance if singular

        // Mean difference
        std::vector<double> diff(n);
        for (int i = 0; i < n; ++i)
            diff[i] = c1.mean[i] - c2.mean[i];

        // Quadratic form: diff^T * avgCovInv * diff
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
        return std::min(JM, 2.0); // clamp to [0, 2]
    }

    // Transformed Divergence between two multivariate normal distributions
    // TD = 2000 * (1 - exp(-D/8)) where D is the divergence
    // D = 0.5 * tr((S1-S2)*(S2inv - S1inv)) + 0.5 * tr((S1inv+S2inv)*(mu1-mu2)*(mu1-mu2)^T)
    double transformedDivergence(const ClassSig& c1, const ClassSig& c2, int n) {
        std::vector<double> s1Inv, s2Inv;
        if (!invertMatrix(c1.covMatrix, s1Inv, n) || !invertMatrix(c2.covMatrix, s2Inv, n))
            return 2000.0;

        // S1 - S2
        std::vector<double> sDiff(n * n);
        for (int i = 0; i < n * n; ++i)
            sDiff[i] = c1.covMatrix[i] - c2.covMatrix[i];

        // S2inv - S1inv
        std::vector<double> sInvDiff(n * n);
        for (int i = 0; i < n * n; ++i)
            sInvDiff[i] = s2Inv[i] - s1Inv[i];

        // S1inv + S2inv
        std::vector<double> sInvSum(n * n);
        for (int i = 0; i < n * n; ++i)
            sInvSum[i] = s1Inv[i] + s2Inv[i];

        double term1 = 0.5 * trace(matMul(sDiff, sInvDiff, n), n);

        // (mu1 - mu2) * (mu1 - mu2)^T
        std::vector<double> diff(n);
        for (int i = 0; i < n; ++i)
            diff[i] = c1.mean[i] - c2.mean[i];

        std::vector<double> diffOuter(n * n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                diffOuter[i * n + j] = diff[i] * diff[j];

        double term2 = 0.5 * trace(matMul(sInvSum, diffOuter, n), n);

        double D = term1 + term2;
        double TD = 2000.0 * (1.0 - std::exp(-D / 8.0));
        return std::min(TD, 2000.0);
    }

public:
    bool execute() override {
        // Read both signature files
        std::vector<ClassSig> classes1, classes2;
        int numBands1 = 0, numBands2 = 0;

        QString sigPath1 = parameter("sig_file1").toString();
        QString sigPath2 = parameter("sig_file2").toString();

        if (!readSigFile(sigPath1, classes1, numBands1)) {
            setError("Failed to read or parse first signature file: " + sigPath1);
            return false;
        }
        if (!readSigFile(sigPath2, classes2, numBands2)) {
            setError("Failed to read or parse second signature file: " + sigPath2);
            return false;
        }

        if (numBands1 != numBands2) {
            setError(QString("Band count mismatch: file 1 has %1 bands, file 2 has %2 bands")
                     .arg(numBands1).arg(numBands2));
            return false;
        }

        int n = numBands1;

        reportProgress(0.2, "Computing pairwise distances...");

        // Open output report
        QString outPath = parameter("output_report").toString();
        std::ofstream out(outPath.toStdString());
        if (!out.is_open()) {
            setError("Failed to open output report: " + outPath);
            return false;
        }

        out << "# Signature Comparison Report\n";
        out << "# File 1: " << sigPath1.toStdString() << "\n";
        out << "# File 2: " << sigPath2.toStdString() << "\n";
        out << "# Number of bands: " << n << "\n";
        out << "# Classes in File 1: " << classes1.size() << "\n";
        out << "# Classes in File 2: " << classes2.size() << "\n";
        out << "#\n";
        out << "# Jeffries-Matusita Distance: range [0, 2], values > 1.9 indicate good separability\n";
        out << "# Transformed Divergence: range [0, 2000], values > 1900 indicate good separability\n";
        out << "#\n\n";

        out << std::fixed << std::setprecision(4);

        out << "Class_File1, Class_File2, Jeffries_Matusita, Transformed_Divergence\n";

        int totalPairs = static_cast<int>(classes1.size() * classes2.size());
        int pairsDone = 0;

        for (const auto& c1 : classes1) {
            for (const auto& c2 : classes2) {
                double jm = jeffriesMatusita(c1, c2, n);
                double td = transformedDivergence(c1, c2, n);

                out << c1.classId << ", " << c2.classId
                    << ", " << jm << ", " << td << "\n";

                pairsDone++;
                if (pairsDone % 10 == 0)
                    reportProgress(0.2 + 0.7 * static_cast<double>(pairsDone) / totalPairs);
            }
        }

        out.close();

        reportProgress(1.0, QString("Comparison report written: %1 class pairs evaluated.")
                       .arg(totalPairs));
        return true;
    }
};

REGISTER_MODULE(SigCompModule)

} // namespace aplaceholder
