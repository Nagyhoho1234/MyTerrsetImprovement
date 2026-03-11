#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class HantsModule : public Module {
public:
    QString name() const override { return "HANTS"; }
    QString description() const override {
        return "Harmonic Analysis of Time Series. Fits Fourier harmonics to per-pixel "
               "time series for temporal smoothing and gap filling.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series_bands", "Time series bands (comma-separated)",
                "Input raster bands representing time steps"),
            ParameterDef::output("output_prefix", "Output prefix"),
            ParameterDef::integer("num_harmonics", "Number of harmonics", 3, 1, 20,
                "Number of Fourier harmonics to fit"),
            ParameterDef::real("fit_error_tolerance", "Fit error tolerance", 0.05, 0.0, 10.0,
                "Tolerance for iterative outlier removal (fraction of range)"),
        };
    }

    bool execute() override {
        QStringList files = parameter("time_series_bands").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : files) f = f.trimmed();
        int n = files.size();
        int numHarm = parameter("num_harmonics").toInt();
        double fitTol = parameter("fit_error_tolerance").toDouble();
        QString outPrefix = parameter("output_prefix").toString();

        if (n < 2 * numHarm + 1) {
            setError(QString("At least %1 time steps required for %2 harmonics")
                     .arg(2 * numHarm + 1).arg(numHarm));
            return false;
        }

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

        // Output: smoothed time series (n bands)
        Raster smoothOut(cols, rows, n, DataType::Float64);
        smoothOut.setGeoTransform(rasters[0]->geoTransform());
        smoothOut.setProjection(rasters[0]->projection());
        smoothOut.setNoDataValue(noData);

        // Also output amplitude and phase per harmonic
        Raster ampOut(cols, rows, numHarm, DataType::Float64);
        ampOut.setGeoTransform(rasters[0]->geoTransform());
        ampOut.setProjection(rasters[0]->projection());
        ampOut.setNoDataValue(noData);

        Raster phaseOut(cols, rows, numHarm, DataType::Float64);
        phaseOut.setGeoTransform(rasters[0]->geoTransform());
        phaseOut.setProjection(rasters[0]->projection());
        phaseOut.setNoDataValue(noData);

        reportProgress(0.1, "Fitting harmonics per pixel...");

        const double twoPi = 2.0 * M_PI;
        int nCoeffs = 2 * numHarm + 1; // a0 + sum(ak*cos + bk*sin)

        for (int64_t i = 0; i < total; ++i) {
            // Gather time series
            std::vector<double> series(n);
            bool valid = true;
            for (int t = 0; t < n; ++t) {
                double val = rasters[t]->data(0)[i];
                if (hasND && val == noData) { valid = false; break; }
                series[t] = val;
            }

            if (!valid) {
                for (int t = 0; t < n; ++t) smoothOut.data(t)[i] = noData;
                for (int h = 0; h < numHarm; ++h) {
                    ampOut.data(h)[i] = noData;
                    phaseOut.data(h)[i] = noData;
                }
                continue;
            }

            // HANTS iterative fitting
            std::vector<bool> useMask(n, true);
            std::vector<double> coeffs(nCoeffs, 0.0);

            // Compute data range for tolerance
            double dmin = *std::min_element(series.begin(), series.end());
            double dmax = *std::max_element(series.begin(), series.end());
            double threshold = fitTol * (dmax - dmin);

            for (int iteration = 0; iteration < 10; ++iteration) {
                // Build and solve normal equations: A^T W A x = A^T W y
                // where A is the Fourier basis matrix, W is diagonal weight matrix
                std::vector<std::vector<double>> AtA(nCoeffs, std::vector<double>(nCoeffs, 0.0));
                std::vector<double> AtY(nCoeffs, 0.0);

                for (int t = 0; t < n; ++t) {
                    if (!useMask[t]) continue;
                    double phase = twoPi * t / n;
                    std::vector<double> basis(nCoeffs);
                    basis[0] = 1.0;
                    for (int h = 0; h < numHarm; ++h) {
                        basis[2 * h + 1] = std::cos((h + 1) * phase);
                        basis[2 * h + 2] = std::sin((h + 1) * phase);
                    }
                    for (int a = 0; a < nCoeffs; ++a) {
                        AtY[a] += basis[a] * series[t];
                        for (int b = a; b < nCoeffs; ++b) {
                            AtA[a][b] += basis[a] * basis[b];
                            if (a != b) AtA[b][a] = AtA[a][b];
                        }
                    }
                }

                // Solve via Gauss elimination
                coeffs = solveLinearSystem(AtA, AtY);

                // Check fit and remove outliers
                bool changed = false;
                for (int t = 0; t < n; ++t) {
                    double phase = twoPi * t / n;
                    double fitted = coeffs[0];
                    for (int h = 0; h < numHarm; ++h) {
                        fitted += coeffs[2 * h + 1] * std::cos((h + 1) * phase);
                        fitted += coeffs[2 * h + 2] * std::sin((h + 1) * phase);
                    }
                    double residual = std::abs(series[t] - fitted);
                    if (residual > threshold && useMask[t]) {
                        useMask[t] = false;
                        changed = true;
                    }
                }
                if (!changed) break;
            }

            // Reconstruct smoothed series and extract amplitude/phase
            for (int t = 0; t < n; ++t) {
                double phase = twoPi * t / n;
                double fitted = coeffs[0];
                for (int h = 0; h < numHarm; ++h) {
                    fitted += coeffs[2 * h + 1] * std::cos((h + 1) * phase);
                    fitted += coeffs[2 * h + 2] * std::sin((h + 1) * phase);
                }
                smoothOut.data(t)[i] = fitted;
            }

            for (int h = 0; h < numHarm; ++h) {
                double ak = coeffs[2 * h + 1];
                double bk = coeffs[2 * h + 2];
                ampOut.data(h)[i] = std::sqrt(ak * ak + bk * bk);
                phaseOut.data(h)[i] = std::atan2(bk, ak);
            }

            if (i % 500000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");

        if (!GdalIO::write(smoothOut, QString("%1_smoothed.tif").arg(outPrefix))) {
            setError("Failed to write smoothed output");
            return false;
        }
        if (!GdalIO::write(ampOut, QString("%1_amplitude.tif").arg(outPrefix))) {
            setError("Failed to write amplitude output");
            return false;
        }
        if (!GdalIO::write(phaseOut, QString("%1_phase.tif").arg(outPrefix))) {
            setError("Failed to write phase output");
            return false;
        }

        reportProgress(1.0);
        return true;
    }

private:
    static std::vector<double> solveLinearSystem(
        std::vector<std::vector<double>> A, std::vector<double> b) {
        int n = (int)b.size();
        // Forward elimination with partial pivoting
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i + 1; j < n; ++j)
                if (std::abs(A[j][i]) > std::abs(A[pivot][i])) pivot = j;
            std::swap(A[i], A[pivot]);
            std::swap(b[i], b[pivot]);
            double diag = A[i][i];
            if (std::abs(diag) < 1e-15) continue;
            for (int j = i + 1; j < n; ++j) {
                double factor = A[j][i] / diag;
                for (int k = i; k < n; ++k) A[j][k] -= factor * A[i][k];
                b[j] -= factor * b[i];
            }
        }
        // Back substitution
        std::vector<double> x(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            double sum = b[i];
            for (int j = i + 1; j < n; ++j) sum -= A[i][j] * x[j];
            x[i] = (std::abs(A[i][i]) > 1e-15) ? sum / A[i][i] : 0.0;
        }
        return x;
    }
};

REGISTER_MODULE(HantsModule)

} // namespace aplaceholder
