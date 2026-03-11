#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>

namespace aplaceholder {

class MssaModule : public Module {
public:
    QString name() const override { return "MSSA"; }
    QString description() const override {
        return "Multichannel Singular Spectrum Analysis. Embeds time series in a trajectory "
               "matrix, computes SVD, and reconstructs signal components.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series_bands", "Time series bands (comma-separated)",
                "Input raster bands representing time steps"),
            ParameterDef::output("output_prefix", "Output prefix"),
            ParameterDef::integer("window_length", "Window (embedding) length", 12, 2, 999,
                "Length of embedding window (L). Should be less than N/2."),
            ParameterDef::integer("num_components", "Number of components", 5, 1, 50,
                "Number of reconstructed components to output"),
        };
    }

    bool execute() override {
        QStringList files = parameter("time_series_bands").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : files) f = f.trimmed();
        int n = files.size();
        int L = parameter("window_length").toInt();
        int numComp = parameter("num_components").toInt();
        QString outPrefix = parameter("output_prefix").toString();

        int K = n - L + 1; // Number of lagged vectors
        if (K < 1 || L < 2) {
            setError(QString("Window length %1 is too large for time series of length %2")
                     .arg(L).arg(n));
            return false;
        }
        numComp = std::min(numComp, std::min(L, K));

        reportProgress(0.0, "Reading time series...");

        std::vector<std::unique_ptr<Raster>> rasters;
        rasters.reserve(n);
        for (int t = 0; t < n; ++t) {
            auto r = GdalIO::read(files[t]);
            if (!r) { setError(QString("Failed to read: %1").arg(files[t])); return false; }
            if (t > 0 && (r->cols() != rasters[0]->cols() || r->rows() != rasters[0]->rows())) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
            rasters.push_back(std::move(r));
        }

        int cols = rasters[0]->cols(), rows = rasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = rasters[0]->noDataValue();
        bool hasND = rasters[0]->hasNoData();

        // Output: one multi-band raster per component (n bands each = reconstructed series)
        std::vector<Raster> compOuts;
        compOuts.reserve(numComp);
        for (int c = 0; c < numComp; ++c) {
            compOuts.emplace_back(cols, rows, n, DataType::Float64);
            compOuts[c].setGeoTransform(rasters[0]->geoTransform());
            compOuts[c].setProjection(rasters[0]->projection());
            compOuts[c].setNoDataValue(noData);
        }

        reportProgress(0.1, "Computing MSSA decomposition per pixel...");

        for (int64_t i = 0; i < total; ++i) {
            std::vector<double> series(n);
            bool valid = true;
            for (int t = 0; t < n; ++t) {
                double val = rasters[t]->data(0)[i];
                if (hasND && val == noData) { valid = false; break; }
                series[t] = val;
            }

            if (!valid) {
                for (int c = 0; c < numComp; ++c)
                    for (int t = 0; t < n; ++t)
                        compOuts[c].data(t)[i] = noData;
                continue;
            }

            // Build trajectory matrix X: L x K
            // X[j][k] = series[j + k]
            // Compute lag-covariance matrix C = X * X^T / K (L x L)
            std::vector<double> covMat(L * L, 0.0);
            for (int r = 0; r < L; ++r) {
                for (int c = r; c < L; ++c) {
                    double sum = 0.0;
                    for (int k = 0; k < K; ++k)
                        sum += series[r + k] * series[c + k];
                    covMat[r * L + c] = sum / K;
                    covMat[c * L + r] = sum / K;
                }
            }

            // Eigendecomposition via power iteration with deflation
            std::vector<std::vector<double>> eigvecs(numComp, std::vector<double>(L));
            std::vector<double> eigvals(numComp, 0.0);

            std::vector<double> covWork(covMat);
            for (int comp = 0; comp < numComp; ++comp) {
                std::vector<double> v(L, 1.0 / std::sqrt((double)L));
                double eigval = 0.0;

                for (int iter = 0; iter < 200; ++iter) {
                    std::vector<double> w(L, 0.0);
                    for (int r = 0; r < L; ++r)
                        for (int c = 0; c < L; ++c)
                            w[r] += covWork[r * L + c] * v[c];

                    double norm = 0.0;
                    for (int r = 0; r < L; ++r) norm += w[r] * w[r];
                    norm = std::sqrt(norm);
                    if (norm < 1e-15) break;
                    eigval = norm;
                    for (int r = 0; r < L; ++r) v[r] = w[r] / norm;
                }

                eigvecs[comp] = v;
                eigvals[comp] = eigval;

                // Deflate
                for (int r = 0; r < L; ++r)
                    for (int c = 0; c < L; ++c)
                        covWork[r * L + c] -= eigval * v[r] * v[c];
            }

            // Reconstruct each component via diagonal averaging
            for (int comp = 0; comp < numComp; ++comp) {
                const auto& ev = eigvecs[comp];

                // Project trajectory matrix onto eigenvector to get principal component
                std::vector<double> pc(K, 0.0);
                for (int k = 0; k < K; ++k) {
                    for (int j = 0; j < L; ++j)
                        pc[k] += ev[j] * series[j + k];
                }

                // Reconstruct rank-1 trajectory matrix: R[j][k] = ev[j] * pc[k]
                // Diagonal averaging to get reconstructed time series
                std::vector<double> recon(n, 0.0);
                std::vector<int> counts(n, 0);
                for (int j = 0; j < L; ++j) {
                    for (int k = 0; k < K; ++k) {
                        recon[j + k] += ev[j] * pc[k];
                        counts[j + k]++;
                    }
                }
                for (int t = 0; t < n; ++t) {
                    recon[t] /= counts[t];
                    compOuts[comp].data(t)[i] = recon[t];
                }
            }

            if (i % 500000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");

        for (int comp = 0; comp < numComp; ++comp) {
            QString path = QString("%1_component_%2.tif").arg(outPrefix).arg(comp + 1);
            if (!GdalIO::write(compOuts[comp], path)) {
                setError(QString("Failed to write: %1").arg(path));
                return false;
            }
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(MssaModule)

} // namespace aplaceholder
