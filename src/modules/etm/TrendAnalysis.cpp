#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class TrendAnalysisModule : public Module {
public:
    QString name() const override { return "TREND_ANALYSIS"; }
    QString description() const override {
        return "Performs linear trend analysis on a time series of raster images. "
               "Computes per-pixel slope and statistical significance of trends.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_series", "Input time series (multi-file)"),
            ParameterDef::output("output_trend", "Output trend image"),
            ParameterDef::output("output_significance", "Output significance image"),
        };
    }

    bool execute() override {
        // Parse comma-separated input file list
        QString inputStr = parameter("input_series").toString();
        QStringList files = inputStr.split(",", Qt::SkipEmptyParts);
        for (auto& f : files) f = f.trimmed();

        int n = files.size();
        if (n < 3) {
            setError("At least 3 time steps are required for trend analysis");
            return false;
        }

        reportProgress(0.0, "Reading input series...");

        // Read all rasters
        std::vector<std::unique_ptr<Raster>> rasters;
        rasters.reserve(n);
        for (int t = 0; t < n; ++t) {
            auto r = GdalIO::read(files[t]);
            if (!r) {
                setError(QString("Failed to read raster: %1").arg(files[t]));
                return false;
            }
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

        Raster trendOut(cols, rows, 1, DataType::Float64);
        trendOut.setGeoTransform(rasters[0]->geoTransform());
        trendOut.setProjection(rasters[0]->projection());
        trendOut.setNoDataValue(noData);

        Raster sigOut(cols, rows, 1, DataType::Float64);
        sigOut.setGeoTransform(rasters[0]->geoTransform());
        sigOut.setProjection(rasters[0]->projection());
        sigOut.setNoDataValue(noData);

        auto& trendData = trendOut.data(0);
        auto& sigData = sigOut.data(0);

        // Precompute time variable sums: t = 0, 1, ..., n-1
        double sumT = 0.0, sumT2 = 0.0;
        for (int t = 0; t < n; ++t) {
            sumT += t;
            sumT2 += static_cast<double>(t) * t;
        }
        double meanT = sumT / n;

        reportProgress(0.1, "Computing per-pixel trends...");

        for (int64_t i = 0; i < total; ++i) {
            // Gather pixel time series
            bool valid = true;
            double sumY = 0.0, sumTY = 0.0, sumY2 = 0.0;

            for (int t = 0; t < n; ++t) {
                double val = rasters[t]->data(0)[i];
                if (hasND && val == noData) {
                    valid = false;
                    break;
                }
                sumY += val;
                sumTY += t * val;
                sumY2 += val * val;
            }

            if (!valid) {
                trendData[i] = noData;
                sigData[i] = noData;
                continue;
            }

            // OLS: y = a + b*t
            // b = (n*sumTY - sumT*sumY) / (n*sumT2 - sumT*sumT)
            double denom = n * sumT2 - sumT * sumT;
            if (denom == 0.0) {
                trendData[i] = 0.0;
                sigData[i] = 1.0;
                continue;
            }

            double b = (n * sumTY - sumT * sumY) / denom;
            double a = (sumY - b * sumT) / n;
            trendData[i] = b;

            // Compute residual sum of squares for t-test
            double sse = 0.0;
            for (int t = 0; t < n; ++t) {
                double predicted = a + b * t;
                double residual = rasters[t]->data(0)[i] - predicted;
                sse += residual * residual;
            }

            if (n <= 2) {
                sigData[i] = 1.0;
                continue;
            }

            double mse = sse / (n - 2);
            double seTDenom = sumT2 - sumT * sumT / n;
            if (seTDenom <= 0.0 || mse <= 0.0) {
                sigData[i] = (b == 0.0) ? 1.0 : 0.0;
                continue;
            }

            double seB = std::sqrt(mse / seTDenom);
            double tStat = b / seB;

            // Approximate two-tailed p-value from t-distribution using
            // the Beta incomplete function approximation
            int df = n - 2;
            double pValue = approxTTestPValue(std::abs(tStat), df);
            sigData[i] = pValue;

            if (i % 500000 == 0)
                reportProgress(0.1 + 0.85 * static_cast<double>(i) / total);
        }

        reportProgress(0.95, "Writing outputs...");
        if (!GdalIO::write(trendOut, parameter("output_trend").toString())) {
            setError("Failed to write trend output");
            return false;
        }
        if (!GdalIO::write(sigOut, parameter("output_significance").toString())) {
            setError("Failed to write significance output");
            return false;
        }

        reportProgress(1.0);
        return true;
    }

private:
    // Approximate two-tailed p-value for t-distribution using the
    // relationship: p = 1 - I_x(df/2, 1/2) where x = df/(df + t^2)
    // Uses a simple continued-fraction / series approximation of the
    // regularized incomplete beta function.
    static double approxTTestPValue(double tAbs, int df) {
        if (df <= 0) return 1.0;
        double x = static_cast<double>(df) / (df + tAbs * tAbs);
        double a = df / 2.0;
        double b = 0.5;

        // Regularized incomplete beta function I_x(a, b) via series expansion
        double betaIx = regIncBeta(x, a, b);
        double pValue = betaIx; // two-tailed
        return std::clamp(pValue, 0.0, 1.0);
    }

    // Regularized incomplete beta function I_x(a, b) using a continued
    // fraction representation (Lentz's method)
    static double regIncBeta(double x, double a, double b) {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;

        // Use symmetry relation if needed for convergence
        if (x > (a + 1.0) / (a + b + 2.0)) {
            return 1.0 - regIncBeta(1.0 - x, b, a);
        }

        double lnBeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
        double front = std::exp(std::log(x) * a + std::log(1.0 - x) * b - lnBeta) / a;

        // Lentz's continued fraction
        const double eps = 1e-14;
        const double tiny = 1e-30;
        double f = 1.0, c = 1.0, d = 1.0 - (a + b) * x / (a + 1.0);
        if (std::abs(d) < tiny) d = tiny;
        d = 1.0 / d;
        f = d;

        for (int m = 1; m <= 200; ++m) {
            // Even step
            double num = m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
            d = 1.0 + num * d; if (std::abs(d) < tiny) d = tiny; d = 1.0 / d;
            c = 1.0 + num / c; if (std::abs(c) < tiny) c = tiny;
            f *= d * c;

            // Odd step
            num = -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0));
            d = 1.0 + num * d; if (std::abs(d) < tiny) d = tiny; d = 1.0 / d;
            c = 1.0 + num / c; if (std::abs(c) < tiny) c = tiny;
            double delta = d * c;
            f *= delta;

            if (std::abs(delta - 1.0) < eps) break;
        }

        return front * f;
    }
};

REGISTER_MODULE(TrendAnalysisModule)

} // namespace aplaceholder
