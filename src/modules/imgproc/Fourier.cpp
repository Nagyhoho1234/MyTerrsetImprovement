#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

namespace aplaceholder {

class FourierModule : public Module {
public:
    QString name() const override { return "FOURIER"; }
    QString description() const override {
        return "Forward and Inverse 2D Fast Fourier Transform. Computes the 2D DFT "
               "using the Cooley-Tukey algorithm (dimensions padded to next power of 2). "
               "Forward outputs magnitude and phase images; inverse reconstructs from them.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image (forward) or magnitude input (inverse)"),
            ParameterDef::output("output_magnitude", "Output magnitude image"),
            ParameterDef::output("output_phase", "Output phase image"),
            ParameterDef::combo("direction", "Transform direction",
                                {"forward", "inverse"}, 0,
                                "Forward: spatial to frequency domain. "
                                "Inverse: frequency domain to spatial domain."),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString magPath = parameter("output_magnitude").toString();
        QString phasePath = parameter("output_phase").toString();
        int dirIdx = parameter("direction").toInt();
        bool forward = (dirIdx == 0);

        if (forward) {
            return executeForward(inputPath, magPath, phasePath);
        } else {
            return executeInverse(inputPath, magPath, phasePath);
        }
    }

private:
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

    bool executeForward(const QString& inputPath, const QString& magPath,
                        const QString& phasePath) {
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

        reportProgress(0.1, "Preparing data...");

        // Build complex data array (zero-padded to power of 2)
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

        reportProgress(0.2, "Computing FFT...");

        // Forward 2D FFT
        fft2d(data, rows, cols, false);

        reportProgress(0.8, "Writing magnitude and phase...");

        // Extract magnitude and phase
        double noData = -9999.0;

        Raster magRaster(cols, rows, 1, DataType::Float64);
        magRaster.setGeoTransform(input->geoTransform());
        magRaster.setProjection(input->projection());
        magRaster.setNoDataValue(noData);

        Raster phaseRaster(cols, rows, 1, DataType::Float64);
        phaseRaster.setGeoTransform(input->geoTransform());
        phaseRaster.setProjection(input->projection());
        phaseRaster.setNoDataValue(noData);

        auto& magData = magRaster.data(0);
        auto& phaseData = phaseRaster.data(0);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                magData[idx] = std::abs(data[r][c]);
                phaseData[idx] = std::arg(data[r][c]);
            }
        }

        if (!GdalIO::write(magRaster, magPath)) {
            setError("Failed to write magnitude output");
            return false;
        }
        if (!GdalIO::write(phaseRaster, phasePath)) {
            setError("Failed to write phase output");
            return false;
        }

        reportProgress(1.0, "Forward FFT complete.");
        return true;
    }

    bool executeInverse(const QString& magPath, const QString& outMagPath,
                        const QString& phasePath) {
        // For inverse: "input" holds magnitude, "output_phase" holds phase input
        // The reconstructed spatial image goes to "output_magnitude"
        auto magRaster = GdalIO::read(magPath);
        auto phaseRaster = GdalIO::read(phasePath);

        if (!magRaster) {
            setError("Failed to read magnitude input: " + magPath);
            return false;
        }
        if (!phaseRaster) {
            setError("Failed to read phase input: " + phasePath);
            return false;
        }

        int cols = magRaster->cols();
        int rows = magRaster->rows();

        if (phaseRaster->cols() != cols || phaseRaster->rows() != rows) {
            setError("Magnitude and phase images must have the same dimensions");
            return false;
        }

        const auto& magData = magRaster->data(0);
        const auto& phaseData = phaseRaster->data(0);

        reportProgress(0.1, "Reconstructing complex data...");

        // Reconstruct complex data from magnitude and phase
        std::vector<std::vector<std::complex<double>>> data(
            rows, std::vector<std::complex<double>>(cols));

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                double mag = magData[idx];
                double phase = phaseData[idx];
                data[r][c] = std::polar(mag, phase);
            }
        }

        reportProgress(0.2, "Computing inverse FFT...");

        // Inverse 2D FFT
        fft2d(data, rows, cols, true);

        reportProgress(0.8, "Writing reconstructed image...");

        double noData = -9999.0;
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(magRaster->geoTransform());
        output.setProjection(magRaster->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;
                outData[idx] = data[r][c].real();
            }
        }

        if (!GdalIO::write(output, outMagPath)) {
            setError("Failed to write reconstructed output");
            return false;
        }

        reportProgress(1.0, "Inverse FFT complete.");
        return true;
    }
};

REGISTER_MODULE(FourierModule)

} // namespace aplaceholder
