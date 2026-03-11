#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

namespace aplaceholder {

class FiltFqModule : public Module {
public:
    QString name() const override { return "FILTFQ"; }
    QString description() const override {
        return "Fourier frequency domain filter. Applies low-pass, high-pass, band-pass, "
               "or band-reject filtering directly to a spatial domain image by performing "
               "a forward FFT, multiplying by a Butterworth filter kernel, and computing "
               "the inverse FFT to produce the filtered result.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output filtered image"),
            ParameterDef::combo("filter_type", "Filter type",
                                {"Low-pass", "High-pass", "Band-pass", "Band-reject"}, 0,
                                "Type of frequency domain filter to apply."),
            ParameterDef::real("cutoff_frequency", "Cutoff frequency", 30.0, 0.001, 999999,
                               "Cutoff frequency in pixels from the center of the frequency domain. "
                               "For band-pass and band-reject, this is the center frequency of the band."),
            ParameterDef::real("bandwidth", "Bandwidth", 15.0, 0.001, 999999,
                               "Bandwidth of the filter in pixels (used for band-pass and band-reject only)."),
            ParameterDef::integer("order", "Butterworth order", 2, 1, 20,
                                  "Order of the Butterworth filter. Higher values produce a sharper cutoff."),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString outputPath = parameter("output").toString();
        int filterIdx = parameter("filter_type").toInt();
        double cutoff = parameter("cutoff_frequency").toDouble();
        double bandwidth = parameter("bandwidth").toDouble();
        int order = parameter("order").toInt();

        auto input = GdalIO::read(inputPath);
        if (!input) {
            setError("Failed to read input image: " + inputPath);
            return false;
        }

        int origCols = input->cols();
        int origRows = input->rows();
        int cols = nextPow2(origCols);
        int rows = nextPow2(origRows);

        bool hasND = input->hasNoData();
        double inputND = input->noDataValue();
        const auto& inData = input->data(0);

        reportProgress(0.05, "Preparing data...");

        // Build complex data array (zero-padded to next power of 2)
        std::vector<std::vector<std::complex<double>>> data(
            rows, std::vector<std::complex<double>>(cols, {0.0, 0.0}));

        for (int r = 0; r < origRows; ++r) {
            for (int c = 0; c < origCols; ++c) {
                double val = inData[static_cast<int64_t>(r) * origCols + c];
                if (hasND && val == inputND)
                    val = 0.0;
                data[r][c] = {val, 0.0};
            }
        }

        reportProgress(0.15, "Computing forward FFT...");

        // Forward 2D FFT
        fft2d(data, rows, cols, false);

        reportProgress(0.45, "Applying frequency filter...");

        // Apply Butterworth filter in the frequency domain
        double centerX = cols / 2.0;
        double centerY = rows / 2.0;

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                // Map to centered frequency coordinates
                // FFT output has DC at (0,0); shift conceptually to center
                int freqR = (r + rows / 2) % rows;
                int freqC = (c + cols / 2) % cols;
                double dx = freqC - centerX;
                double dy = freqR - centerY;
                double dist = std::sqrt(dx * dx + dy * dy);

                double H = computeFilter(filterIdx, dist, cutoff, bandwidth, order);

                data[r][c] *= H;
            }
        }

        reportProgress(0.55, "Computing inverse FFT...");

        // Inverse 2D FFT
        fft2d(data, rows, cols, true);

        reportProgress(0.85, "Writing output...");

        // Extract real part back to raster (original dimensions)
        double noData = input->hasNoData() ? input->noDataValue() : -9999.0;

        Raster output(origCols, origRows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);

        auto& dst = output.data(0);
        for (int r = 0; r < origRows; ++r) {
            for (int c = 0; c < origCols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * origCols + c;
                if (hasND && inData[idx] == inputND) {
                    dst[idx] = noData;
                } else {
                    dst[idx] = data[r][c].real();
                }
            }
        }

        if (!GdalIO::write(output, outputPath)) {
            setError("Failed to write output image");
            return false;
        }

        reportProgress(1.0, "Frequency filtering complete.");
        return true;
    }

private:
    // Compute Butterworth filter gain for a given distance from center
    static double computeFilter(int filterIdx, double dist, double cutoff,
                                double bandwidth, int order) {
        double H = 0.0;

        switch (filterIdx) {
        case 0: // Low-pass
            if (cutoff > 0.0) {
                H = 1.0 / (1.0 + std::pow(dist / cutoff, 2.0 * order));
            }
            break;

        case 1: // High-pass
            if (dist > 0.0) {
                H = 1.0 / (1.0 + std::pow(cutoff / dist, 2.0 * order));
            } else {
                H = 0.0; // DC component fully suppressed
            }
            break;

        case 2: // Band-pass
            if (dist > 0.0 && bandwidth > 0.0) {
                double term = (dist * dist - cutoff * cutoff) / (dist * bandwidth);
                H = 1.0 / (1.0 + std::pow(term, 2.0 * order));
            } else {
                H = 0.0;
            }
            break;

        case 3: // Band-reject (notch)
            if (dist > 0.0 && bandwidth > 0.0) {
                double term = (dist * bandwidth) / (dist * dist - cutoff * cutoff);
                H = 1.0 / (1.0 + std::pow(term, 2.0 * order));
            } else {
                H = 1.0; // DC component passes through
            }
            break;
        }

        return H;
    }

    // Next power of 2 >= n
    static int nextPow2(int n) {
        int p = 1;
        while (p < n) p <<= 1;
        return p;
    }

    // 1D FFT in-place (Cooley-Tukey radix-2 DIT)
    static void fft1d(std::vector<std::complex<double>>& x, bool inverse) {
        int N = static_cast<int>(x.size());
        if (N <= 1) return;

        // Bit-reversal permutation
        for (int i = 1, j = 0; i < N; ++i) {
            int bit = N >> 1;
            for (; j & bit; bit >>= 1) {
                j ^= bit;
            }
            j ^= bit;
            if (i < j) std::swap(x[i], x[j]);
        }

        // Cooley-Tukey butterfly
        for (int len = 2; len <= N; len <<= 1) {
            double angle = 2.0 * M_PI / len * (inverse ? -1.0 : 1.0);
            std::complex<double> wlen(std::cos(angle), std::sin(angle));
            for (int i = 0; i < N; i += len) {
                std::complex<double> w(1.0, 0.0);
                for (int j = 0; j < len / 2; ++j) {
                    std::complex<double> u = x[i + j];
                    std::complex<double> v = x[i + j + len / 2] * w;
                    x[i + j] = u + v;
                    x[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }

        if (inverse) {
            for (auto& val : x)
                val /= static_cast<double>(N);
        }
    }

    // 2D FFT via row-column decomposition
    static void fft2d(std::vector<std::vector<std::complex<double>>>& data,
                      int rows, int cols, bool inverse) {
        // Transform rows
        for (int r = 0; r < rows; ++r) {
            fft1d(data[r], inverse);
        }

        // Transform columns
        std::vector<std::complex<double>> col(rows);
        for (int c = 0; c < cols; ++c) {
            for (int r = 0; r < rows; ++r)
                col[r] = data[r][c];
            fft1d(col, inverse);
            for (int r = 0; r < rows; ++r)
                data[r][c] = col[r];
        }
    }
};

REGISTER_MODULE(FiltFqModule)

} // namespace aplaceholder
