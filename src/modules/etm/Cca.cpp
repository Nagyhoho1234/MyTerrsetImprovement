#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>

namespace aplaceholder {

class CcaModule : public Module {
public:
    QString name() const override { return "CCA"; }
    QString description() const override {
        return "Canonical Correlation Analysis. Finds linear combinations of two sets "
               "of variables that are maximally correlated.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("set1_bands", "Set 1 bands (comma-separated)",
                "First set of raster bands for CCA"),
            ParameterDef::file("set2_bands", "Set 2 bands (comma-separated)",
                "Second set of raster bands for CCA"),
            ParameterDef::output("output_prefix", "Output prefix"),
            ParameterDef::integer("num_components", "Number of components", 3, 1, 100,
                "Number of canonical components to extract"),
        };
    }

    bool execute() override {
        // Parse input band lists
        QStringList set1Files = parameter("set1_bands").toString().split(",", Qt::SkipEmptyParts);
        QStringList set2Files = parameter("set2_bands").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : set1Files) f = f.trimmed();
        for (auto& f : set2Files) f = f.trimmed();

        int p = set1Files.size();
        int q = set2Files.size();
        int numComp = parameter("num_components").toInt();
        QString outPrefix = parameter("output_prefix").toString();

        if (p < 1 || q < 1) {
            setError("Both sets must have at least one band");
            return false;
        }
        numComp = std::min(numComp, std::min(p, q));

        reportProgress(0.0, "Reading set 1 bands...");

        // Read set 1
        std::vector<std::unique_ptr<Raster>> set1;
        set1.reserve(p);
        for (int i = 0; i < p; ++i) {
            auto r = GdalIO::read(set1Files[i]);
            if (!r) { setError(QString("Failed to read: %1").arg(set1Files[i])); return false; }
            if (i > 0 && (r->cols() != set1[0]->cols() || r->rows() != set1[0]->rows())) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
            set1.push_back(std::move(r));
        }

        reportProgress(0.1, "Reading set 2 bands...");

        // Read set 2
        std::vector<std::unique_ptr<Raster>> set2;
        set2.reserve(q);
        for (int i = 0; i < q; ++i) {
            auto r = GdalIO::read(set2Files[i]);
            if (!r) { setError(QString("Failed to read: %1").arg(set2Files[i])); return false; }
            if (r->cols() != set1[0]->cols() || r->rows() != set1[0]->rows()) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
            set2.push_back(std::move(r));
        }

        int cols = set1[0]->cols(), rows = set1[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = set1[0]->noDataValue();
        bool hasND = set1[0]->hasNoData();
        int totalVars = p + q;

        reportProgress(0.2, "Computing covariance matrices...");

        // Build valid pixel mask and compute means
        std::vector<bool> validMask(total, true);
        int64_t validCount = 0;

        std::vector<double> means(totalVars, 0.0);
        for (int64_t i = 0; i < total; ++i) {
            bool valid = true;
            for (int v = 0; v < p && valid; ++v) {
                double val = set1[v]->data(0)[i];
                if (hasND && val == noData) valid = false;
            }
            for (int v = 0; v < q && valid; ++v) {
                double val = set2[v]->data(0)[i];
                if (hasND && val == noData) valid = false;
            }
            validMask[i] = valid;
            if (valid) {
                for (int v = 0; v < p; ++v) means[v] += set1[v]->data(0)[i];
                for (int v = 0; v < q; ++v) means[p + v] += set2[v]->data(0)[i];
                validCount++;
            }
        }

        if (validCount < totalVars + 1) {
            setError("Insufficient valid pixels for CCA");
            return false;
        }
        for (int v = 0; v < totalVars; ++v) means[v] /= validCount;

        // Compute covariance matrix (full combined)
        // C = [C11 C12; C21 C22]
        std::vector<std::vector<double>> cov(totalVars, std::vector<double>(totalVars, 0.0));
        for (int64_t i = 0; i < total; ++i) {
            if (!validMask[i]) continue;
            std::vector<double> vals(totalVars);
            for (int v = 0; v < p; ++v) vals[v] = set1[v]->data(0)[i] - means[v];
            for (int v = 0; v < q; ++v) vals[p + v] = set2[v]->data(0)[i] - means[p + v];
            for (int a = 0; a < totalVars; ++a)
                for (int b = a; b < totalVars; ++b) {
                    cov[a][b] += vals[a] * vals[b];
                    if (a != b) cov[b][a] += vals[a] * vals[b];
                }
        }
        for (int a = 0; a < totalVars; ++a)
            for (int b = 0; b < totalVars; ++b)
                cov[a][b] /= (validCount - 1);

        // Extract sub-matrices C11, C22, C12
        auto matMul = [](const std::vector<std::vector<double>>& A,
                         const std::vector<std::vector<double>>& B) {
            int ra = (int)A.size(), ca = (int)A[0].size(), cb = (int)B[0].size();
            std::vector<std::vector<double>> C(ra, std::vector<double>(cb, 0.0));
            for (int i = 0; i < ra; ++i)
                for (int j = 0; j < cb; ++j)
                    for (int k = 0; k < ca; ++k)
                        C[i][j] += A[i][k] * B[k][j];
            return C;
        };

        // Simple symmetric matrix inverse via Gauss-Jordan
        auto matInv = [](std::vector<std::vector<double>> M) {
            int n = (int)M.size();
            std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
            for (int i = 0; i < n; ++i) I[i][i] = 1.0;
            for (int i = 0; i < n; ++i) {
                int pivot = i;
                for (int j = i + 1; j < n; ++j)
                    if (std::abs(M[j][i]) > std::abs(M[pivot][i])) pivot = j;
                std::swap(M[i], M[pivot]);
                std::swap(I[i], I[pivot]);
                double diag = M[i][i];
                if (std::abs(diag) < 1e-15) diag = 1e-15;
                for (int j = 0; j < n; ++j) { M[i][j] /= diag; I[i][j] /= diag; }
                for (int j = 0; j < n; ++j) {
                    if (j == i) continue;
                    double factor = M[j][i];
                    for (int k = 0; k < n; ++k) {
                        M[j][k] -= factor * M[i][k];
                        I[j][k] -= factor * I[i][k];
                    }
                }
            }
            return I;
        };

        // Extract submatrices
        auto subMat = [&](int r0, int c0, int nr, int nc) {
            std::vector<std::vector<double>> S(nr, std::vector<double>(nc));
            for (int i = 0; i < nr; ++i)
                for (int j = 0; j < nc; ++j)
                    S[i][j] = cov[r0 + i][c0 + j];
            return S;
        };

        auto C11 = subMat(0, 0, p, p);
        auto C22 = subMat(p, p, q, q);
        auto C12 = subMat(0, p, p, q);

        // Transpose C12 to get C21
        std::vector<std::vector<double>> C21(q, std::vector<double>(p));
        for (int i = 0; i < q; ++i)
            for (int j = 0; j < p; ++j)
                C21[i][j] = C12[j][i];

        // CCA: solve C11^{-1} C12 C22^{-1} C21 a = lambda^2 a
        // via power iteration for top components
        auto C11inv = matInv(C11);
        auto C22inv = matInv(C22);
        auto M1 = matMul(C11inv, matMul(C12, matMul(C22inv, C21)));

        reportProgress(0.4, "Extracting canonical components via power iteration...");

        // Power iteration with deflation for eigenvalues/vectors of M1
        std::vector<std::vector<double>> eigvecs_a(numComp, std::vector<double>(p, 0.0));
        std::vector<double> eigenvalues(numComp, 0.0);

        auto M1work = M1; // working copy for deflation
        for (int comp = 0; comp < numComp; ++comp) {
            // Initialize eigenvector
            std::vector<double> v(p, 1.0 / std::sqrt(p));
            double eigval = 0.0;

            for (int iter = 0; iter < 200; ++iter) {
                // Multiply: w = M1work * v
                std::vector<double> w(p, 0.0);
                for (int i = 0; i < p; ++i)
                    for (int j = 0; j < p; ++j)
                        w[i] += M1work[i][j] * v[j];

                // Compute eigenvalue estimate
                double norm = 0.0;
                for (int i = 0; i < p; ++i) norm += w[i] * w[i];
                norm = std::sqrt(norm);
                if (norm < 1e-15) break;
                eigval = norm;

                // Normalize
                for (int i = 0; i < p; ++i) v[i] = w[i] / norm;
            }

            eigvecs_a[comp] = v;
            eigenvalues[comp] = eigval;

            // Deflate: M1work = M1work - eigval * v * v^T
            for (int i = 0; i < p; ++i)
                for (int j = 0; j < p; ++j)
                    M1work[i][j] -= eigval * v[i] * v[j];
        }

        // Compute set2 weights: b = C22^{-1} C21 a / ||...||
        std::vector<std::vector<double>> eigvecs_b(numComp, std::vector<double>(q, 0.0));
        auto C22invC21 = matMul(C22inv, C21);
        for (int comp = 0; comp < numComp; ++comp) {
            for (int i = 0; i < q; ++i) {
                for (int j = 0; j < p; ++j)
                    eigvecs_b[comp][i] += C22invC21[i][j] * eigvecs_a[comp][j];
            }
            double norm = 0.0;
            for (int i = 0; i < q; ++i) norm += eigvecs_b[comp][i] * eigvecs_b[comp][i];
            norm = std::sqrt(norm);
            if (norm > 1e-15) {
                for (int i = 0; i < q; ++i) eigvecs_b[comp][i] /= norm;
            }
        }

        reportProgress(0.6, "Computing canonical variate images...");

        // Project each pixel onto canonical variates
        for (int comp = 0; comp < numComp; ++comp) {
            Raster outU(cols, rows, 1, DataType::Float64);
            outU.setGeoTransform(set1[0]->geoTransform());
            outU.setProjection(set1[0]->projection());
            outU.setNoDataValue(noData);

            Raster outV(cols, rows, 1, DataType::Float64);
            outV.setGeoTransform(set1[0]->geoTransform());
            outV.setProjection(set1[0]->projection());
            outV.setNoDataValue(noData);

            auto& uData = outU.data(0);
            auto& vData = outV.data(0);

            for (int64_t i = 0; i < total; ++i) {
                if (!validMask[i]) {
                    uData[i] = noData;
                    vData[i] = noData;
                    continue;
                }
                double u = 0.0, v = 0.0;
                for (int j = 0; j < p; ++j)
                    u += eigvecs_a[comp][j] * (set1[j]->data(0)[i] - means[j]);
                for (int j = 0; j < q; ++j)
                    v += eigvecs_b[comp][j] * (set2[j]->data(0)[i] - means[p + j]);
                uData[i] = u;
                vData[i] = v;
            }

            QString uPath = QString("%1_U%2.tif").arg(outPrefix).arg(comp + 1);
            QString vPath = QString("%1_V%2.tif").arg(outPrefix).arg(comp + 1);
            if (!GdalIO::write(outU, uPath)) {
                setError(QString("Failed to write: %1").arg(uPath));
                return false;
            }
            if (!GdalIO::write(outV, vPath)) {
                setError(QString("Failed to write: %1").arg(vPath));
                return false;
            }

            reportProgress(0.6 + 0.35 * (comp + 1.0) / numComp);
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(CcaModule)

} // namespace aplaceholder
