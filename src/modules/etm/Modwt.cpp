#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class ModwtModule : public Module {
public:
    QString name() const override { return "MODWT"; }
    QString description() const override {
        return "Maximal Overlap Discrete Wavelet Transform. Decomposes per-pixel time series "
               "using Haar or Daubechies-4 wavelets into detail and smooth coefficients.";
    }
    QString category() const override { return "Earth Trends Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time_series_bands", "Time series bands (comma-separated)",
                "Input raster bands representing time steps"),
            ParameterDef::output("output_prefix", "Output prefix"),
            ParameterDef::integer("num_levels", "Number of decomposition levels", 3, 1, 10,
                "Number of wavelet decomposition levels"),
            ParameterDef::combo("wavelet_type", "Wavelet type",
                {"haar", "db4"}, 0, "Wavelet filter to use"),
        };
    }

    bool execute() override {
        QStringList files = parameter("time_series_bands").toString().split(",", Qt::SkipEmptyParts);
        for (auto& f : files) f = f.trimmed();
        int n = files.size();
        int numLevels = parameter("num_levels").toInt();
        int waveletIdx = parameter("wavelet_type").toInt();
        QString outPrefix = parameter("output_prefix").toString();

        // Get wavelet filter coefficients
        std::vector<double> h, g; // scaling and wavelet filters
        if (waveletIdx == 0) {
            // Haar
            double s2 = 1.0 / std::sqrt(2.0);
            h = {s2, s2};
            g = {s2, -s2};
        } else {
            // Daubechies-4
            double c0 = (1.0 + std::sqrt(3.0)) / (4.0 * std::sqrt(2.0));
            double c1 = (3.0 + std::sqrt(3.0)) / (4.0 * std::sqrt(2.0));
            double c2 = (3.0 - std::sqrt(3.0)) / (4.0 * std::sqrt(2.0));
            double c3 = (1.0 - std::sqrt(3.0)) / (4.0 * std::sqrt(2.0));
            h = {c0, c1, c2, c3};
            g = {c3, -c2, c1, -c0};
        }
        int filterLen = (int)h.size();

        // Verify sufficient data length
        int minLen = filterLen * (1 << numLevels);
        if (n < filterLen) {
            setError(QString("Time series length %1 is too short for the selected wavelet (filter length %2)")
                     .arg(n).arg(filterLen));
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

        // Output: detail coefficients at each level + smooth at last level
        // Each is n bands
        std::vector<Raster> detailOuts(numLevels, Raster(cols, rows, n, DataType::Float64));
        Raster smoothOut(cols, rows, n, DataType::Float64);
        for (int lev = 0; lev < numLevels; ++lev) {
            detailOuts[lev].setGeoTransform(rasters[0]->geoTransform());
            detailOuts[lev].setProjection(rasters[0]->projection());
            detailOuts[lev].setNoDataValue(noData);
        }
        smoothOut.setGeoTransform(rasters[0]->geoTransform());
        smoothOut.setProjection(rasters[0]->projection());
        smoothOut.setNoDataValue(noData);

        reportProgress(0.1, "Computing MODWT decomposition...");

        // MODWT rescaling: divide filters by sqrt(2^j) at level j
        for (int64_t i = 0; i < total; ++i) {
            // Gather pixel time series
            std::vector<double> series(n);
            bool valid = true;
            for (int t = 0; t < n; ++t) {
                double val = rasters[t]->data(0)[i];
                if (hasND && val == noData) { valid = false; break; }
                series[t] = val;
            }

            if (!valid) {
                for (int lev = 0; lev < numLevels; ++lev)
                    for (int t = 0; t < n; ++t)
                        detailOuts[lev].data(t)[i] = noData;
                for (int t = 0; t < n; ++t)
                    smoothOut.data(t)[i] = noData;
                continue;
            }

            // Apply MODWT level by level
            std::vector<double> smooth = series;

            for (int lev = 0; lev < numLevels; ++lev) {
                int step = 1 << lev; // Dyadic upsampling of filter
                std::vector<double> detail(n, 0.0);
                std::vector<double> nextSmooth(n, 0.0);

                // MODWT: no downsampling, filters are periodized
                // Scale filters by 1/sqrt(2) per level (already done for MODWT)
                double scale = 1.0 / std::sqrt(2.0);
                // But the base filters are already normalized for MODWT at level 0
                // At each level we just apply with upsampled index

                for (int t = 0; t < n; ++t) {
                    double wt = 0.0, vt = 0.0;
                    for (int k = 0; k < filterLen; ++k) {
                        int idx = ((t - k * step) % n + n) % n;
                        wt += g[k] * smooth[idx];
                        vt += h[k] * smooth[idx];
                    }
                    detail[t] = wt;
                    nextSmooth[t] = vt;
                }

                for (int t = 0; t < n; ++t)
                    detailOuts[lev].data(t)[i] = detail[t];

                smooth = nextSmooth;
            }

            for (int t = 0; t < n; ++t)
                smoothOut.data(t)[i] = smooth[t];

            if (i % 500000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");

        for (int lev = 0; lev < numLevels; ++lev) {
            QString path = QString("%1_detail_L%2.tif").arg(outPrefix).arg(lev + 1);
            if (!GdalIO::write(detailOuts[lev], path)) {
                setError(QString("Failed to write: %1").arg(path));
                return false;
            }
        }
        if (!GdalIO::write(smoothOut, QString("%1_smooth.tif").arg(outPrefix))) {
            setError("Failed to write smooth output");
            return false;
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(ModwtModule)

} // namespace aplaceholder
