#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class CanCorModule : public Module {
public:
    QString name() const override { return "CANCOR"; }
    QString description() const override {
        return "Canonical Correlation Analysis between two sets of variables. "
               "Finds pairs of linear combinations with maximum correlation.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_group1", "Input raster group 1 (multi-band image)"),
            ParameterDef::file("input_group2", "Input raster group 2 (multi-band image)"),
            ParameterDef::integer("num_variates", "Number of canonical variates", 3, 1, 256,
                "Number of canonical variate pairs to retain"),
            ParameterDef::output("output_prefix", "Output prefix for canonical variate images"),
        };
    }

    bool execute() override {
        auto raster1 = GdalIO::read(parameter("input_group1").toString());
        if (!raster1) {
            setError("Failed to read input raster group 1");
            return false;
        }

        auto raster2 = GdalIO::read(parameter("input_group2").toString());
        if (!raster2) {
            setError("Failed to read input raster group 2");
            return false;
        }

        if (raster1->cols() != raster2->cols() || raster1->rows() != raster2->rows()) {
            setError("Input raster groups must have the same dimensions");
            return false;
        }

        int cols = raster1->cols(), rows = raster1->rows();
        int p = raster1->bands();  // bands in group 1
        int q = raster2->bands();  // bands in group 2
        int numVar = parameter("num_variates").toInt();
        QString outPrefix = parameter("output_prefix").toString();

        int minPQ = std::min(p, q);
        if (numVar > minPQ) numVar = minPQ;

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = raster1->hasNoData();
        double noData = raster1->noDataValue();

        // Collect band data
        std::vector<const std::vector<double>*> bandsX(p), bandsY(q);
        for (int b = 0; b < p; ++b)
            bandsX[b] = &raster1->data(b);
        for (int b = 0; b < q; ++b)
            bandsY[b] = &raster2->data(b);

        reportProgress(0.1, "Computing means...");

        // Compute means
        std::vector<double> meanX(p, 0.0), meanY(q, 0.0);
        int64_t validCount = 0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bandsX[0])[i] == noData) continue;
            validCount++;
            for (int b = 0; b < p; ++b)
                meanX[b] += (*bandsX[b])[i];
            for (int b = 0; b < q; ++b)
                meanY[b] += (*bandsY[b])[i];
        }
        if (validCount < 2) {
            setError("Insufficient valid pixels for canonical correlation analysis");
            return false;
        }
        for (int b = 0; b < p; ++b) meanX[b] /= validCount;
        for (int b = 0; b < q; ++b) meanY[b] /= validCount;

        reportProgress(0.2, "Computing covariance matrices...");

        // Compute Sxx (p x p), Syy (q x q), Sxy (p x q)
        std::vector<double> Sxx(p * p, 0.0), Syy(q * q, 0.0), Sxy(p * q, 0.0);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bandsX[0])[i] == noData) continue;
            for (int b1 = 0; b1 < p; ++b1) {
                double dx1 = (*bandsX[b1])[i] - meanX[b1];
                for (int b2 = b1; b2 < p; ++b2) {
                    double dx2 = (*bandsX[b2])[i] - meanX[b2];
                    Sxx[b1 * p + b2] += dx1 * dx2;
                }
                for (int b2 = 0; b2 < q; ++b2) {
                    double dy2 = (*bandsY[b2])[i] - meanY[b2];
                    Sxy[b1 * q + b2] += dx1 * dy2;
                }
            }
            for (int b1 = 0; b1 < q; ++b1) {
                double dy1 = (*bandsY[b1])[i] - meanY[b1];
                for (int b2 = b1; b2 < q; ++b2) {
                    double dy2 = (*bandsY[b2])[i] - meanY[b2];
                    Syy[b1 * q + b2] += dy1 * dy2;
                }
            }
        }
        double denom = validCount - 1;
        for (int b1 = 0; b1 < p; ++b1)
            for (int b2 = b1; b2 < p; ++b2) {
                Sxx[b1 * p + b2] /= denom;
                Sxx[b2 * p + b1] = Sxx[b1 * p + b2];
            }
        for (int b1 = 0; b1 < q; ++b1)
            for (int b2 = b1; b2 < q; ++b2) {
                Syy[b1 * q + b2] /= denom;
                Syy[b2 * q + b1] = Syy[b1 * q + b2];
            }
        for (int b1 = 0; b1 < p; ++b1)
            for (int b2 = 0; b2 < q; ++b2)
                Sxy[b1 * q + b2] /= denom;

        reportProgress(0.3, "Inverting covariance matrices...");

        // Gauss-Jordan matrix inversion helper
        auto invertMatrix = [&](const std::vector<double>& mat, int dim) -> std::vector<double> {
            std::vector<double> inv(dim * dim, 0.0);
            std::vector<double> aug(dim * 2 * dim, 0.0);
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j)
                    aug[i * 2 * dim + j] = mat[i * dim + j];
                aug[i * 2 * dim + dim + i] = 1.0;
            }
            for (int i = 0; i < dim; ++i) {
                int maxRow = i;
                for (int k = i + 1; k < dim; ++k)
                    if (std::abs(aug[k * 2 * dim + i]) > std::abs(aug[maxRow * 2 * dim + i]))
                        maxRow = k;
                if (maxRow != i)
                    for (int j = 0; j < 2 * dim; ++j)
                        std::swap(aug[i * 2 * dim + j], aug[maxRow * 2 * dim + j]);
                double pivot = aug[i * 2 * dim + i];
                if (std::abs(pivot) < 1e-15) return {};
                for (int j = 0; j < 2 * dim; ++j)
                    aug[i * 2 * dim + j] /= pivot;
                for (int k = 0; k < dim; ++k) {
                    if (k == i) continue;
                    double factor = aug[k * 2 * dim + i];
                    for (int j = 0; j < 2 * dim; ++j)
                        aug[k * 2 * dim + j] -= factor * aug[i * 2 * dim + j];
                }
            }
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                    inv[i * dim + j] = aug[i * 2 * dim + dim + j];
            return inv;
        };

        auto SxxInv = invertMatrix(Sxx, p);
        if (SxxInv.empty()) {
            setError("Covariance matrix of group 1 is singular");
            return false;
        }
        auto SyyInv = invertMatrix(Syy, q);
        if (SyyInv.empty()) {
            setError("Covariance matrix of group 2 is singular");
            return false;
        }

        reportProgress(0.4, "Computing canonical correlation matrix...");

        // Compute Syx = Sxy^T (q x p)
        std::vector<double> Syx(q * p, 0.0);
        for (int i = 0; i < p; ++i)
            for (int j = 0; j < q; ++j)
                Syx[j * p + i] = Sxy[i * q + j];

        // Matrix multiply helper: C = A(r x s) * B(s x t)
        auto matMul = [](const std::vector<double>& A, const std::vector<double>& B,
                         int r, int s, int t) -> std::vector<double> {
            std::vector<double> C(r * t, 0.0);
            for (int i = 0; i < r; ++i)
                for (int j = 0; j < t; ++j)
                    for (int k = 0; k < s; ++k)
                        C[i * t + j] += A[i * s + k] * B[k * t + j];
            return C;
        };

        // Form M = SxxInv * Sxy * SyyInv * Syx (p x p)
        auto T1 = matMul(SxxInv, Sxy, p, p, q);      // p x q
        auto T2 = matMul(T1, SyyInv, p, q, q);        // p x q
        auto M = matMul(T2, Syx, p, q, p);            // p x p

        // Symmetrize for Jacobi
        std::vector<double> A(p * p, 0.0);
        for (int i = 0; i < p; ++i)
            for (int j = 0; j < p; ++j)
                A[i * p + j] = 0.5 * (M[i * p + j] + M[j * p + i]);

        reportProgress(0.5, "Computing eigenvectors (Jacobi rotation)...");

        // Jacobi eigenvalue decomposition on A (p x p)
        std::vector<double> V(p * p, 0.0);
        for (int i = 0; i < p; ++i) V[i * p + i] = 1.0;

        const int maxSweeps = 100 * p * p;
        for (int sweep = 0; sweep < maxSweeps; ++sweep) {
            double maxOff = 0.0;
            int pi = 0, qi = 1;
            for (int i = 0; i < p; ++i)
                for (int j = i + 1; j < p; ++j)
                    if (std::abs(A[i * p + j]) > maxOff) {
                        maxOff = std::abs(A[i * p + j]);
                        pi = i; qi = j;
                    }
            if (maxOff < 1e-12) break;

            double app = A[pi * p + pi];
            double aqq = A[qi * p + qi];
            double apq = A[pi * p + qi];

            double theta;
            if (std::abs(app - aqq) < 1e-15)
                theta = M_PI / 4.0;
            else
                theta = 0.5 * std::atan2(2.0 * apq, app - aqq);

            double c = std::cos(theta);
            double s = std::sin(theta);

            std::vector<double> ap(p), aq(p);
            for (int j = 0; j < p; ++j) {
                ap[j] = A[pi * p + j];
                aq[j] = A[qi * p + j];
            }
            for (int j = 0; j < p; ++j) {
                A[pi * p + j] =  c * ap[j] + s * aq[j];
                A[qi * p + j] = -s * ap[j] + c * aq[j];
            }
            std::vector<double> colp(p), colq(p);
            for (int i = 0; i < p; ++i) {
                colp[i] = A[i * p + pi];
                colq[i] = A[i * p + qi];
            }
            for (int i = 0; i < p; ++i) {
                A[i * p + pi] =  c * colp[i] + s * colq[i];
                A[i * p + qi] = -s * colp[i] + c * colq[i];
            }
            for (int i = 0; i < p; ++i) {
                double vp = V[i * p + pi];
                double vq = V[i * p + qi];
                V[i * p + pi] =  c * vp + s * vq;
                V[i * p + qi] = -s * vp + c * vq;
            }
        }

        // Eigenvalues (squared canonical correlations) from diagonal
        std::vector<double> eigenvalues(p);
        for (int i = 0; i < p; ++i)
            eigenvalues[i] = A[i * p + i];

        // Sort by eigenvalue descending
        std::vector<int> order(p);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return eigenvalues[a] > eigenvalues[b];
        });

        reportProgress(0.6, "Projecting data onto canonical variates...");

        // Output: U (group 1 variates) and V_out (group 2 variates) interleaved
        // Band layout: U1, V1, U2, V2, ...
        int totalBands = numVar * 2;
        Raster output(cols, rows, totalBands, DataType::Float64);
        output.setGeoTransform(raster1->geoTransform());
        output.setProjection(raster1->projection());
        output.setNoDataValue(-9999.0);

        for (int var = 0; var < numVar; ++var) {
            int eigIdx = order[var];

            // Eigenvector a (weight for X)
            std::vector<double> aVec(p);
            for (int b = 0; b < p; ++b)
                aVec[b] = V[b * p + eigIdx];

            // Compute b = SyyInv * Syx * a (weight for Y)
            std::vector<double> tmp(q, 0.0);
            for (int i = 0; i < q; ++i)
                for (int k = 0; k < p; ++k)
                    tmp[i] += Syx[i * p + k] * aVec[k];
            std::vector<double> bVec(q, 0.0);
            for (int i = 0; i < q; ++i)
                for (int k = 0; k < q; ++k)
                    bVec[i] += SyyInv[i * q + k] * tmp[k];

            // Project X onto a (U variate)
            auto& outU = output.data(var * 2);
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && (*bandsX[0])[i] == noData) {
                    outU[i] = -9999.0;
                    continue;
                }
                double val = 0.0;
                for (int b = 0; b < p; ++b)
                    val += ((*bandsX[b])[i] - meanX[b]) * aVec[b];
                outU[i] = val;
            }

            // Project Y onto b (V variate)
            auto& outV = output.data(var * 2 + 1);
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && (*bandsX[0])[i] == noData) {
                    outV[i] = -9999.0;
                    continue;
                }
                double val = 0.0;
                for (int b = 0; b < q; ++b)
                    val += ((*bandsY[b])[i] - meanY[b]) * bVec[b];
                outV[i] = val;
            }

            reportProgress(0.6 + 0.3 * (var + 1.0) / numVar);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outPrefix);
    }
};

REGISTER_MODULE(CanCorModule)

} // namespace aplaceholder
