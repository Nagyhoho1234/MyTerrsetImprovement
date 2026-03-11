#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class TfaModule : public Module {
public:
    QString name() const override { return "TFA"; }
    QString description() const override {
        return "Time series Fourier Analysis. Decomposes an image time series into "
               "harmonic components (amplitude and phase for each frequency).";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster group file (multi-band time series)"),
            ParameterDef::integer("num_harmonics", "Number of harmonics", 3, 1, 128,
                "Number of harmonic components to extract"),
            ParameterDef::output("output_prefix", "Output prefix for harmonic images"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster group file");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int N = raster->bands();  // number of time steps
        int numHarmonics = parameter("num_harmonics").toInt();
        QString outPrefix = parameter("output_prefix").toString();

        if (N < 2) {
            setError("Input must have at least 2 bands (time steps)");
            return false;
        }

        // Maximum harmonics is N/2 (Nyquist limit)
        int maxHarmonics = N / 2;
        if (numHarmonics > maxHarmonics) numHarmonics = maxHarmonics;

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        // Collect band data
        std::vector<const std::vector<double>*> bands(N);
        for (int b = 0; b < N; ++b)
            bands[b] = &raster->data(b);

        reportProgress(0.1, "Precomputing Fourier basis functions...");

        // Output layout: band 0 = mean (DC component),
        // then for each harmonic k=1..numHarmonics: amplitude band, phase band
        int outBands = 1 + numHarmonics * 2;
        Raster output(cols, rows, outBands, DataType::Float64);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(-9999.0);

        // Precompute cosine and sine tables
        std::vector<std::vector<double>> cosTable(numHarmonics, std::vector<double>(N));
        std::vector<std::vector<double>> sinTable(numHarmonics, std::vector<double>(N));
        for (int h = 0; h < numHarmonics; ++h) {
            int freq = h + 1;
            for (int t = 0; t < N; ++t) {
                double angle = 2.0 * M_PI * freq * t / N;
                cosTable[h][t] = std::cos(angle);
                sinTable[h][t] = std::sin(angle);
            }
        }

        reportProgress(0.2, "Decomposing time series per pixel...");

        auto& meanBand = output.data(0);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*bands[0])[i] == noData) {
                meanBand[i] = -9999.0;
                for (int h = 0; h < numHarmonics; ++h) {
                    output.data(1 + h * 2)[i] = -9999.0;
                    output.data(2 + h * 2)[i] = -9999.0;
                }
                continue;
            }

            // Compute mean (DC component)
            double mean = 0.0;
            for (int t = 0; t < N; ++t)
                mean += (*bands[t])[i];
            mean /= N;
            meanBand[i] = mean;

            // Compute Fourier coefficients for each harmonic
            for (int h = 0; h < numHarmonics; ++h) {
                double realPart = 0.0;
                double imagPart = 0.0;
                for (int t = 0; t < N; ++t) {
                    double val = (*bands[t])[i];
                    realPart += val * cosTable[h][t];
                    imagPart += val * sinTable[h][t];
                }
                realPart *= 2.0 / N;
                imagPart *= 2.0 / N;

                double amplitude = std::sqrt(realPart * realPart + imagPart * imagPart);
                double phase = std::atan2(imagPart, realPart);

                output.data(1 + h * 2)[i] = amplitude;
                output.data(2 + h * 2)[i] = phase;
            }

            // Report progress periodically
            if (i % (total / 10 + 1) == 0)
                reportProgress(0.2 + 0.7 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outPrefix);
    }
};

REGISTER_MODULE(TfaModule)

} // namespace aplaceholder
