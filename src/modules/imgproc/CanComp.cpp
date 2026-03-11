#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>

namespace aplaceholder {

class CanCompModule : public Module {
public:
    QString name() const override { return "CANCOMP"; }
    QString description() const override {
        return "Canonical Components Analysis. Finds components that maximize "
               "between-group variance relative to within-group variance.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster group (multi-band image)"),
            ParameterDef::file("training_classes", "Training classes raster"),
            ParameterDef::integer("num_components", "Number of canonical components", 3, 1, 256,
                "Number of canonical components to retain"),
            ParameterDef::output("output_prefix", "Output prefix for component images"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster group");
            return false;
        }

        auto classRaster = GdalIO::read(parameter("training_classes").toString());
        if (!classRaster) {
            setError("Failed to read training classes raster");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int n = raster->bands();
        int numComp = parameter("num_components").toInt();
        QString outPrefix = parameter("output_prefix").toString();

        if (numComp > n) numComp = n;

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        // Collect band data
        std::vector<const std::vector<double>*> bands(n);
        for (int b = 0; b < n; ++b)
            bands[b] = &raster->data(b);

        const auto& classData = classRaster->data(0);

        reportProgress(0.1, "Identifying classes and computing statistics...");

        // Identify unique classes (excluding nodata / background)
        std::map<int, int64_t> classCounts;
        std::map<int, std::vector<double>> classMeans;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) continue;
            int cls = static_cast<int>(classData[i]);
            if (cls <= 0) continue;
            classCounts[cls]++;
        }

        int numClasses = static_cast<int>(classCounts.size());
        if (numClasses < 2) {
            setError("At least two classes are required for canonical components analysis");
            return false;
        }

        // Compute per-class means
        for (auto& kv : classCounts)
            classMeans[kv.first].assign(n, 0.0);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) continue;
            int cls = static_cast<int>(classData[i]);
            if (classCounts.find(cls) == classCounts.end()) continue;
            for (int b = 0; b < n; ++b)
                classMeans[cls][b] += (*bands[b])[i];
        }
        for (auto& kv : classMeans)
            for (int b = 0; b < n; ++b)
                kv.second[b] /= classCounts[kv.first];

        // Compute grand mean (weighted by class counts)
        int64_t totalValid = 0;
        for (auto& kv : classCounts) totalValid += kv.second;

        std::vector<double> grandMean(n, 0.0);
        for (auto& kv : classMeans)
            for (int b = 0; b < n; ++b)
                grandMean[b] += classCounts[kv.first] * kv.second[b];
        for (int b = 0; b < n; ++b)
            grandMean[b] /= totalValid;

        reportProgress(0.2, "Computing within-group covariance matrix...");

        // Within-group covariance matrix W (n x n)
        std::vector<double> W(n * n, 0.0);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) continue;
            int cls = static_cast<int>(classData[i]);
            if (classCounts.find(cls) == classCounts.end()) continue;
            for (int b1 = 0; b1 < n; ++b1) {
                double d1 = (*bands[b1])[i] - classMeans[cls][b1];
                for (int b2 = b1; b2 < n; ++b2) {
                    double d2 = (*bands[b2])[i] - classMeans[cls][b2];
                    W[b1 * n + b2] += d1 * d2;
                }
            }
        }
        for (int b1 = 0; b1 < n; ++b1) {
            for (int b2 = b1; b2 < n; ++b2) {
                W[b1 * n + b2] /= (totalValid - numClasses);
                W[b2 * n + b1] = W[b1 * n + b2];
            }
        }

        reportProgress(0.3, "Computing between-group covariance matrix...");

        // Between-group covariance matrix B (n x n)
        std::vector<double> B(n * n, 0.0);
        for (auto& kv : classMeans) {
            int cls = kv.first;
            int64_t cnt = classCounts[cls];
            for (int b1 = 0; b1 < n; ++b1) {
                double d1 = kv.second[b1] - grandMean[b1];
                for (int b2 = b1; b2 < n; ++b2) {
                    double d2 = kv.second[b2] - grandMean[b2];
                    B[b1 * n + b2] += cnt * d1 * d2;
                }
            }
        }
        for (int b1 = 0; b1 < n; ++b1) {
            for (int b2 = b1; b2 < n; ++b2) {
                B[b1 * n + b2] /= (numClasses - 1);
                B[b2 * n + b1] = B[b1 * n + b2];
            }
        }

        reportProgress(0.4, "Inverting within-group covariance matrix...");

        // Invert W via Gauss-Jordan elimination
        std::vector<double> Winv(n * n, 0.0);
        std::vector<double> aug(n * 2 * n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                aug[i * 2 * n + j] = W[i * n + j];
            aug[i * 2 * n + n + i] = 1.0;
        }
        for (int i = 0; i < n; ++i) {
            int maxRow = i;
            for (int k = i + 1; k < n; ++k)
                if (std::abs(aug[k * 2 * n + i]) > std::abs(aug[maxRow * 2 * n + i]))
                    maxRow = k;
            if (maxRow != i)
                for (int j = 0; j < 2 * n; ++j)
                    std::swap(aug[i * 2 * n + j], aug[maxRow * 2 * n + j]);

            double pivot = aug[i * 2 * n + i];
            if (std::abs(pivot) < 1e-15) {
                setError("Within-group covariance matrix is singular");
                return false;
            }
            for (int j = 0; j < 2 * n; ++j)
                aug[i * 2 * n + j] /= pivot;
            for (int k = 0; k < n; ++k) {
                if (k == i) continue;
                double factor = aug[k * 2 * n + i];
                for (int j = 0; j < 2 * n; ++j)
                    aug[k * 2 * n + j] -= factor * aug[i * 2 * n + j];
            }
        }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                Winv[i * n + j] = aug[i * 2 * n + n + j];

        // Compute M = W^{-1} * B (n x n)
        std::vector<double> M(n * n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k)
                    M[i * n + j] += Winv[i * n + k] * B[k * n + j];

        // Symmetrize for Jacobi: A = (M + M^T) / 2
        std::vector<double> A(n * n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                A[i * n + j] = 0.5 * (M[i * n + j] + M[j * n + i]);

        reportProgress(0.5, "Computing eigenvectors (Jacobi rotation)...");

        // Jacobi eigenvalue decomposition
        std::vector<double> V(n * n, 0.0);
        for (int i = 0; i < n; ++i) V[i * n + i] = 1.0;

        const int maxSweeps = 100 * n * n;
        for (int sweep = 0; sweep < maxSweeps; ++sweep) {
            double maxOff = 0.0;
            int p = 0, q = 1;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j)
                    if (std::abs(A[i * n + j]) > maxOff) {
                        maxOff = std::abs(A[i * n + j]);
                        p = i; q = j;
                    }
            if (maxOff < 1e-12) break;

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

            std::vector<double> ap(n), aq(n);
            for (int j = 0; j < n; ++j) {
                ap[j] = A[p * n + j];
                aq[j] = A[q * n + j];
            }
            for (int j = 0; j < n; ++j) {
                A[p * n + j] =  c * ap[j] + s * aq[j];
                A[q * n + j] = -s * ap[j] + c * aq[j];
            }
            std::vector<double> colp(n), colq(n);
            for (int i = 0; i < n; ++i) {
                colp[i] = A[i * n + p];
                colq[i] = A[i * n + q];
            }
            for (int i = 0; i < n; ++i) {
                A[i * n + p] =  c * colp[i] + s * colq[i];
                A[i * n + q] = -s * colp[i] + c * colq[i];
            }
            for (int i = 0; i < n; ++i) {
                double vp = V[i * n + p];
                double vq = V[i * n + q];
                V[i * n + p] =  c * vp + s * vq;
                V[i * n + q] = -s * vp + c * vq;
            }
        }

        // Eigenvalues from diagonal, sort descending
        std::vector<double> eigenvalues(n);
        for (int i = 0; i < n; ++i)
            eigenvalues[i] = A[i * n + i];

        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return eigenvalues[a] > eigenvalues[b];
        });

        reportProgress(0.6, "Projecting data onto canonical components...");

        // Project data onto top-K eigenvectors
        Raster output(cols, rows, numComp, DataType::Float64);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(-9999.0);

        for (int comp = 0; comp < numComp; ++comp) {
            int eigIdx = order[comp];
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
                    double centered = (*bands[b])[i] - grandMean[b];
                    val += centered * eigvec[b];
                }
                outData[i] = val;
            }

            reportProgress(0.6 + 0.3 * (comp + 1.0) / numComp);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outPrefix);
    }
};

REGISTER_MODULE(CanCompModule)

} // namespace aplaceholder
