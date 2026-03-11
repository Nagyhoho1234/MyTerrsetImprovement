#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class FiltFreqModule : public Module {
public:
    QString name() const override { return "FILTFREQ"; }
    QString description() const override {
        return "Frequency domain filtering using Butterworth filters. Applies low-pass, "
               "high-pass, or band-pass filtering to magnitude and phase images from "
               "a Fourier transform.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_magnitude", "Input magnitude image"),
            ParameterDef::file("input_phase", "Input phase image"),
            ParameterDef::output("output_magnitude", "Output filtered magnitude image"),
            ParameterDef::output("output_phase", "Output filtered phase image"),
            ParameterDef::combo("filter_type", "Filter type",
                                {"low_pass", "high_pass", "band_pass"}, 0,
                                "Type of frequency domain filter to apply."),
            ParameterDef::real("cutoff_frequency", "Cutoff frequency", 30.0, 0.001, 999999,
                               "Cutoff frequency (distance from center in pixels). "
                               "For band_pass, this is the center frequency."),
            ParameterDef::integer("order", "Butterworth order", 2, 1, 20,
                                  "Order of the Butterworth filter. Higher = sharper cutoff."),
        };
    }

    bool execute() override {
        QString magInPath = parameter("input_magnitude").toString();
        QString phaseInPath = parameter("input_phase").toString();
        QString magOutPath = parameter("output_magnitude").toString();
        QString phaseOutPath = parameter("output_phase").toString();
        int filterIdx = parameter("filter_type").toInt();
        double cutoff = parameter("cutoff_frequency").toDouble();
        int order = parameter("order").toInt();

        auto magIn = GdalIO::read(magInPath);
        auto phaseIn = GdalIO::read(phaseInPath);

        if (!magIn) {
            setError("Failed to read magnitude input: " + magInPath);
            return false;
        }
        if (!phaseIn) {
            setError("Failed to read phase input: " + phaseInPath);
            return false;
        }

        int cols = magIn->cols();
        int rows = magIn->rows();

        if (phaseIn->cols() != cols || phaseIn->rows() != rows) {
            setError("Magnitude and phase images must have the same dimensions");
            return false;
        }

        const auto& magData = magIn->data(0);
        const auto& phaseData = phaseIn->data(0);

        double noData = -9999.0;

        Raster magOut(cols, rows, 1, DataType::Float64);
        magOut.setGeoTransform(magIn->geoTransform());
        magOut.setProjection(magIn->projection());
        magOut.setNoDataValue(noData);

        Raster phaseOut(cols, rows, 1, DataType::Float64);
        phaseOut.setGeoTransform(phaseIn->geoTransform());
        phaseOut.setProjection(phaseIn->projection());
        phaseOut.setNoDataValue(noData);

        auto& magOutData = magOut.data(0);
        auto& phaseOutData = phaseOut.data(0);

        // Center of frequency domain
        double centerX = cols / 2.0;
        double centerY = rows / 2.0;

        // Bandwidth for band-pass (use cutoff/2 as bandwidth)
        double bandwidth = cutoff * 0.5;
        if (bandwidth < 1.0) bandwidth = 1.0;

        int64_t total = static_cast<int64_t>(cols) * rows;
        int64_t processed = 0;

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;

                // Distance from center of frequency domain
                double dx = c - centerX;
                double dy = r - centerY;
                double dist = std::sqrt(dx * dx + dy * dy);

                // Compute Butterworth filter gain H
                double H = 0.0;

                switch (filterIdx) {
                case 0: // low_pass
                    // H = 1 / (1 + (dist/cutoff)^(2*order))
                    if (cutoff > 0.0) {
                        H = 1.0 / (1.0 + std::pow(dist / cutoff, 2.0 * order));
                    }
                    break;

                case 1: // high_pass
                    // H = 1 / (1 + (cutoff/dist)^(2*order))
                    if (dist > 0.0) {
                        H = 1.0 / (1.0 + std::pow(cutoff / dist, 2.0 * order));
                    } else {
                        H = 0.0; // DC component fully suppressed
                    }
                    break;

                case 2: // band_pass
                    // Butterworth band-pass: H = 1 / (1 + ((dist^2 - cutoff^2) / (dist * bandwidth))^(2*order))
                    if (dist > 0.0 && bandwidth > 0.0) {
                        double term = (dist * dist - cutoff * cutoff) / (dist * bandwidth);
                        H = 1.0 / (1.0 + std::pow(term, 2.0 * order));
                    } else {
                        H = 0.0;
                    }
                    break;
                }

                // Apply filter: multiply magnitude by H, phase passes through
                magOutData[idx] = magData[idx] * H;
                phaseOutData[idx] = phaseData[idx];

                ++processed;
                if (processed % 1000000 == 0)
                    reportProgress(static_cast<double>(processed) / total);
            }
        }

        reportProgress(0.9, "Writing filtered outputs...");

        if (!GdalIO::write(magOut, magOutPath)) {
            setError("Failed to write filtered magnitude output");
            return false;
        }
        if (!GdalIO::write(phaseOut, phaseOutPath)) {
            setError("Failed to write filtered phase output");
            return false;
        }

        reportProgress(1.0, "Frequency filtering complete.");
        return true;
    }
};

REGISTER_MODULE(FiltFreqModule)

} // namespace aplaceholder
