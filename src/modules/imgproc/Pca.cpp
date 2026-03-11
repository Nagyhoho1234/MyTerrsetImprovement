#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class PcaModule : public Module {
public:
    QString name() const override { return "PCA"; }
    QString description() const override {
        return "Principal Components Analysis. Transforms a multi-band image into "
               "uncorrelated principal component bands ordered by variance explained.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input multi-band image"),
            ParameterDef::output("output", "Output principal components image"),
            ParameterDef::integer("num_components", "Number of components", 3, 1, 256,
                "Number of principal components to retain"),
            ParameterDef::boolean("standardize", "Standardize input bands", true,
                "Standardize bands to zero mean and unit variance before analysis"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input image");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int n = raster->bands();  // number of bands = matrix dimension
        int numComp = parameter("num_components").toInt();
        bool standardize = parameter("standardize").toBool();

        if (numComp > n) numComp = n;

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        // Collect band data
        std::vector<const std::vector<double>*> bands(n);
        for (int b = 0; b < n; ++b)
            bands[b] = &raster->data(b);

        reportProgress(0.1, "Computing band statistics...");

        // Step 1: Compute band means
        std::vector<double> means(n, 0.0);
        int64_t validCount = 0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) continue;
            validCount++;
            for (int b = 0; b < n; ++b)
                means[b] += (*bands[b])[i];
        }
        if (validCount == 0) {
            setError("No valid pixels found in input image");
            return false;
        }
        for (int b = 0; b < n; ++b)
            means[b] /= validCount;

        reportProgress(0.2, "Computing covariance matrix...");

        // Step 2: Compute covariance matrix (row-major: covMat[i*n + j])
        std::vector<double> covMat(n * n, 0.0);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) continue;
            for (int b1 = 0; b1 < n; ++b1) {
                double d1 = (*bands[b1])[i] - means[b1];
                for (int b2 = b1; b2 < n; ++b2) {
                    double d2 = (*bands[b2])[i] - means[b2];
                    covMat[b1 * n + b2] += d1 * d2;
                }
            }
        }
        // Normalize and symmetrize
        for (int b1 = 0; b1 < n; ++b1) {
            for (int b2 = b1; b2 < n; ++b2) {
                covMat[b1 * n + b2] /= (validCount - 1);
                covMat[b2 * n + b1] = covMat[b1 * n + b2];
            }
        }

        // If standardize: convert covariance to correlation matrix
        std::vector<double> stddevs(n, 1.0);
        if (standardize) {
            for (int b = 0; b < n; ++b) {
                stddevs[b] = std::sqrt(covMat[b * n + b]);
                if (stddevs[b] < 1e-15) stddevs[b] = 1.0;
            }
            for (int b1 = 0; b1 < n; ++b1)
                for (int b2 = 0; b2 < n; ++b2)
                    covMat[b1 * n + b2] /= (stddevs[b1] * stddevs[b2]);
        }

        reportProgress(0.3, "Computing eigenvalues (Jacobi rotation)...");

        // Step 3: Jacobi eigenvalue algorithm on the symmetric matrix covMat
        // A will be diagonalized; V accumulates the rotation matrices (columns = eigenvectors)
        std::vector<double> A(covMat);
        std::vector<double> V(n * n, 0.0);
        for (int i = 0; i < n; ++i) V[i * n + i] = 1.0;

        const int maxSweeps = 100 * n * n;
        for (int sweep = 0; sweep < maxSweeps; ++sweep) {
            // Find largest off-diagonal element
            double maxOff = 0.0;
            int p = 0, q = 1;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j)
                    if (std::abs(A[i * n + j]) > maxOff) {
                        maxOff = std::abs(A[i * n + j]);
                        p = i; q = j;
                    }
            if (maxOff < 1e-12) break;

            // Compute Jacobi rotation parameters
            double app = A[p * n + p];
            double aqq = A[q * n + q];
            double apq = A[p * n + q];

            double theta;
            if (std::abs(app - aqq) < 1e-15)
                theta = M_PI / 4.0;
            else
                theta = 0.5 * std::atan2(2.0 * apq, app - aqq);

            double c = std::cos(theta);
            double s = std::sin(theta);

            // Apply similarity transform: A' = G^T A G
            // Only rows/cols p and q are affected.
            // First save the affected rows
            std::vector<double> ap(n), aq(n);
            for (int j = 0; j < n; ++j) {
                ap[j] = A[p * n + j];
                aq[j] = A[q * n + j];
            }

            // Update row p and row q (left multiply by G^T)
            for (int j = 0; j < n; ++j) {
                A[p * n + j] =  c * ap[j] + s * aq[j];
                A[q * n + j] = -s * ap[j] + c * aq[j];
            }

            // Now update col p and col q (right multiply by G)
            // Save current columns first
            std::vector<double> colp(n), colq(n);
            for (int i = 0; i < n; ++i) {
                colp[i] = A[i * n + p];
                colq[i] = A[i * n + q];
            }
            for (int i = 0; i < n; ++i) {
                A[i * n + p] =  c * colp[i] + s * colq[i];
                A[i * n + q] = -s * colp[i] + c * colq[i];
            }

            // Update eigenvector matrix: V <- V * G
            for (int i = 0; i < n; ++i) {
                double vp = V[i * n + p];
                double vq = V[i * n + q];
                V[i * n + p] =  c * vp + s * vq;
                V[i * n + q] = -s * vp + c * vq;
            }
        }

        // Eigenvalues are diagonal of A
        std::vector<double> eigenvalues(n);
        for (int i = 0; i < n; ++i)
            eigenvalues[i] = A[i * n + i];

        // Step 4: Sort eigenvectors by eigenvalue descending
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return eigenvalues[a] > eigenvalues[b];
        });

        reportProgress(0.5, "Projecting data onto principal components...");

        // Step 5: Project data onto top-K eigenvectors
        Raster output(cols, rows, numComp, DataType::Float64);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(-9999.0);

        for (int comp = 0; comp < numComp; ++comp) {
            int eigIdx = order[comp];
            // Eigenvector is column eigIdx of V
            std::vector<double> eigvec(n);
            for (int b = 0; b < n; ++b)
                eigvec[b] = V[b * n + eigIdx];

            auto& outData = output.data(comp);
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && (*bands[0])[i] == noData) {
                    outData[i] = -9999.0;
                    continue;
                }

                double val = 0.0;
                for (int b = 0; b < n; ++b) {
                    double centered = (*bands[b])[i] - means[b];
                    if (standardize)
                        centered /= stddevs[b];
                    val += centered * eigvec[b];
                }
                outData[i] = val;
            }

            reportProgress(0.5 + 0.4 * (comp + 1.0) / numComp);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PcaModule)

} // namespace aplaceholder
