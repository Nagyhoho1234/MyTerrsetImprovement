#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class StaModule : public Module {
public:
    QString name() const override { return "STA"; }
    QString description() const override {
        return "Seasonal Trend Analysis. Decomposes time series into trend, seasonal amplitude, "
               "and seasonal phase using harmonic regression with time-varying coefficients.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series_bands", "Time series bands (comma-separated)",
                "Input raster bands representing time steps"),
            ParameterDef::output("output_trend", "Output trend raster"),
            ParameterDef::output("output_amplitude", "Output seasonal amplitude raster"),
            ParameterDef::output("output_phase", "Output seasonal phase raster"),
            ParameterDef::integer("period", "Seasonal period", 12, 2, 365,
                "Number of observations per seasonal cycle (e.g. 12 for monthly)"),
        };
    }

    bool execute() override {
        QStringList files = parameter("time_series_bands").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : files) f = f.trimmed();
        int n = files.size();
        int period = parameter("period").toInt();
        QString outTrend = parameter("output_trend").toString();
        QString outAmplitude = parameter("output_amplitude").toString();
        QString outPhase = parameter("output_phase").toString();

        if (n < period) {
            setError(QString("At least %1 time steps required for period %2").arg(period).arg(period));
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

        // Output rasters: trend slope, amplitude trend, phase trend
        Raster trendOut(cols, rows, 1, DataType::Float64);
        trendOut.setGeoTransform(rasters[0]->geoTransform());
        trendOut.setProjection(rasters[0]->projection());
        trendOut.setNoDataValue(noData);

        Raster ampOut(cols, rows, 1, DataType::Float64);
        ampOut.setGeoTransform(rasters[0]->geoTransform());
        ampOut.setProjection(rasters[0]->projection());
        ampOut.setNoDataValue(noData);

        Raster phaseOut(cols, rows, 1, DataType::Float64);
        phaseOut.setGeoTransform(rasters[0]->geoTransform());
        phaseOut.setProjection(rasters[0]->projection());
        phaseOut.setNoDataValue(noData);

        auto& trendData = trendOut.data(0);
        auto& ampData = ampOut.data(0);
        auto& phaseData = phaseOut.data(0);

        reportProgress(0.1, "Computing Seasonal Trend Analysis...");

        const double twoPi = 2.0 * M_PI;

        // STA model: y(t) = a0 + a1*t + a2*cos(2pi*t/P) + a3*sin(2pi*t/P)
        //                  + a4*t*cos(2pi*t/P) + a5*t*sin(2pi*t/P) + e
        // This allows trend in mean (a1), seasonal amplitude and phase,
        // and time-varying seasonal amplitude/phase (a4, a5).
        const int nCoeffs = 6;

        // Precompute design matrix
        std::vector<std::vector<double>> A(n, std::vector<double>(nCoeffs));
        for (int t = 0; t < n; ++t) {
            double omega = twoPi * t / period;
            A[t][0] = 1.0;
            A[t][1] = t;
            A[t][2] = std::cos(omega);
            A[t][3] = std::sin(omega);
            A[t][4] = t * std::cos(omega);
            A[t][5] = t * std::sin(omega);
        }

        // Precompute A^T A and prepare for solving
        std::vector<std::vector<double>> AtA(nCoeffs, std::vector<double>(nCoeffs, 0.0));
        for (int i = 0; i < nCoeffs; ++i)
            for (int j = i; j < nCoeffs; ++j) {
                double sum = 0.0;
                for (int t = 0; t < n; ++t) sum += A[t][i] * A[t][j];
                AtA[i][j] = sum;
                AtA[j][i] = sum;
            }

        // Invert AtA once (same for all pixels since no missing data handling per pixel)
        // For pixels with nodata, skip entirely
        auto AtAinv = invertMatrix(AtA);

        for (int64_t idx = 0; idx < total; ++idx) {
            std::vector<double> series(n);
            bool valid = true;
            for (int t = 0; t < n; ++t) {
                double val = rasters[t]->data(0)[idx];
                if (hasND && val == noData) { valid = false; break; }
                series[t] = val;
            }

            if (!valid) {
                trendData[idx] = noData;
                ampData[idx] = noData;
                phaseData[idx] = noData;
                continue;
            }

            // Compute A^T y
            std::vector<double> AtY(nCoeffs, 0.0);
            for (int i = 0; i < nCoeffs; ++i)
                for (int t = 0; t < n; ++t)
                    AtY[i] += A[t][i] * series[t];

            // Coefficients = (AtA)^{-1} AtY
            std::vector<double> coeffs(nCoeffs, 0.0);
            for (int i = 0; i < nCoeffs; ++i)
                for (int j = 0; j < nCoeffs; ++j)
                    coeffs[i] += AtAinv[i][j] * AtY[j];

            // Trend: linear trend slope = a1
            trendData[idx] = coeffs[1];

            // Mean seasonal amplitude = sqrt(a2^2 + a3^2)
            ampData[idx] = std::sqrt(coeffs[2] * coeffs[2] + coeffs[3] * coeffs[3]);

            // Mean seasonal phase = atan2(a3, a2)
            phaseData[idx] = std::atan2(coeffs[3], coeffs[2]);

            if (idx % 500000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(idx) / total);
        }

        reportProgress(0.9, "Writing outputs...");
        if (!GdalIO::write(trendOut, outTrend)) {
            setError("Failed to write trend output");
            return false;
        }
        if (!GdalIO::write(ampOut, outAmplitude)) {
            setError("Failed to write amplitude output");
            return false;
        }
        if (!GdalIO::write(phaseOut, outPhase)) {
            setError("Failed to write phase output");
            return false;
        }

        reportProgress(1.0);
        return true;
    }

private:
    static std::vector<std::vector<double>> invertMatrix(std::vector<std::vector<double>> M) {
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
    }
};

REGISTER_MODULE(StaModule)

} // namespace aplaceholder
