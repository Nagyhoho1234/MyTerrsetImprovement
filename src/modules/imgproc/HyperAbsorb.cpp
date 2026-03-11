#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class HyperAbsorbModule : public Module {
public:
    QString name() const override { return "HYPERABSORB"; }
    QString description() const override {
        return "Absorption feature analysis for hyperspectral data. Performs continuum removal "
               "and measures absorption depth, width, and area at each pixel for a specified "
               "absorption center wavelength.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::file("wavelengths_file", "Wavelengths file (CSV)",
                "CSV with band center wavelengths in order matching the input bands"),
            ParameterDef::output("output_depth", "Output absorption depth image",
                "Output raster with absorption depth at the specified center wavelength"),
            ParameterDef::output("output_width", "Output absorption width image",
                "Output raster with absorption feature width (wavelength units)"),
            ParameterDef::real("absorption_center", "Absorption center wavelength",
                2.2, 0.0, 100000.0,
                "Wavelength of the absorption feature center (same units as wavelengths file)"),
        };
    }

    bool execute() override {
        // --------------------------------------------------------------------
        // 1. Read input bands
        // --------------------------------------------------------------------
        QString bandsParam = parameter("bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths) p = p.trimmed();

        if (bandPaths.isEmpty()) {
            setError("No input bands specified");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> bandRasters;
        for (const auto& path : bandPaths) {
            auto r = GdalIO::read(path);
            if (!r) {
                setError("Failed to read band image: " + path);
                return false;
            }
            bandRasters.push_back(std::move(r));
        }

        int numBands = static_cast<int>(bandRasters.size());
        int cols = bandRasters[0]->cols();
        int rows = bandRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = bandRasters[0]->hasNoData();
        double noData = bandRasters[0]->noDataValue();

        for (int b = 1; b < numBands; ++b) {
            if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                setError("All input bands must have the same dimensions");
                return false;
            }
        }

        // --------------------------------------------------------------------
        // 2. Read wavelengths file (CSV, single row or column)
        // --------------------------------------------------------------------
        QString wavPath = parameter("wavelengths_file").toString();
        std::vector<double> wavelengths;

        {
            std::ifstream file(wavPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open wavelengths file: " + wavPath);
                return false;
            }

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;

                std::istringstream iss(line);
                std::string token;
                while (std::getline(iss, token, ',')) {
                    size_t start = token.find_first_not_of(" \t");
                    size_t end = token.find_last_not_of(" \t\r\n");
                    if (start != std::string::npos && end != std::string::npos) {
                        try {
                            wavelengths.push_back(std::stod(token.substr(start, end - start + 1)));
                        } catch (...) {}
                    }
                }
            }
        }

        if (static_cast<int>(wavelengths.size()) < numBands) {
            setError(QString("Wavelengths file has %1 values but %2 bands provided")
                     .arg(wavelengths.size()).arg(numBands));
            return false;
        }
        wavelengths.resize(numBands);

        double absCenter = parameter("absorption_center").toDouble();

        // Find the band index closest to the absorption center
        int centerBand = 0;
        double minDiff = std::abs(wavelengths[0] - absCenter);
        for (int b = 1; b < numBands; ++b) {
            double diff = std::abs(wavelengths[b] - absCenter);
            if (diff < minDiff) {
                minDiff = diff;
                centerBand = b;
            }
        }

        reportProgress(0.05, "Analyzing absorption features...");

        // --------------------------------------------------------------------
        // 3. Collect band data pointers
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // --------------------------------------------------------------------
        // 4. Create output rasters
        // --------------------------------------------------------------------
        Raster outputDepth(cols, rows, 1, DataType::Float64);
        outputDepth.setGeoTransform(bandRasters[0]->geoTransform());
        outputDepth.setProjection(bandRasters[0]->projection());
        outputDepth.setNoDataValue(noData);
        auto& depthData = outputDepth.data(0);

        Raster outputWidth(cols, rows, 1, DataType::Float64);
        outputWidth.setGeoTransform(bandRasters[0]->geoTransform());
        outputWidth.setProjection(bandRasters[0]->projection());
        outputWidth.setNoDataValue(noData);
        auto& widthData = outputWidth.data(0);

        // --------------------------------------------------------------------
        // 5. For each pixel: continuum removal + absorption depth/width
        // --------------------------------------------------------------------
        std::vector<double> spectrum(numBands);

        for (int64_t idx = 0; idx < total; ++idx) {
            if (hasND && (*bands[0])[idx] == noData) {
                depthData[idx] = noData;
                widthData[idx] = noData;
                continue;
            }

            for (int b = 0; b < numBands; ++b)
                spectrum[b] = (*bands[b])[idx];

            // ----------------------------------------------------------------
            // Continuum removal: compute convex hull upper envelope
            // Use a simple approach: build a piecewise-linear continuum
            // connecting local maxima that form the convex hull of the spectrum
            // ----------------------------------------------------------------
            std::vector<double> continuum(numBands);

            // Build upper convex hull over (wavelength, reflectance) points
            // Start from the first point, always connect to the point that
            // maintains the hull (steepest descent from current point)
            std::vector<int> hullPts;
            hullPts.push_back(0);

            int current = 0;
            while (current < numBands - 1) {
                // Find next hull point: the point that makes the largest slope
                // from current point (or least negative slope)
                double bestSlope = -std::numeric_limits<double>::max();
                int bestIdx = current + 1;

                for (int b = current + 1; b < numBands; ++b) {
                    double dw = wavelengths[b] - wavelengths[current];
                    if (std::abs(dw) < 1e-15) continue;
                    double slope = (spectrum[b] - spectrum[current]) / dw;
                    if (slope > bestSlope) {
                        bestSlope = slope;
                        bestIdx = b;
                    }
                }
                hullPts.push_back(bestIdx);
                current = bestIdx;
            }

            // Ensure last band is included
            if (hullPts.back() != numBands - 1)
                hullPts.push_back(numBands - 1);

            // Interpolate continuum between hull points
            for (size_t h = 0; h + 1 < hullPts.size(); ++h) {
                int b1 = hullPts[h];
                int b2 = hullPts[h + 1];
                double w1 = wavelengths[b1], w2 = wavelengths[b2];
                double v1 = spectrum[b1], v2 = spectrum[b2];

                for (int b = b1; b <= b2; ++b) {
                    double t = (std::abs(w2 - w1) > 1e-15)
                               ? (wavelengths[b] - w1) / (w2 - w1)
                               : 0.0;
                    continuum[b] = v1 + t * (v2 - v1);
                }
            }

            // Continuum-removed spectrum
            // CR = spectrum / continuum (1.0 = on continuum, < 1.0 = absorption)
            // Absorption depth = 1 - CR at center

            double crCenter = (std::abs(continuum[centerBand]) > 1e-15)
                              ? spectrum[centerBand] / continuum[centerBand]
                              : 1.0;
            double depth = 1.0 - crCenter;
            if (depth < 0.0) depth = 0.0;

            // Absorption width: find the range of bands around center where
            // continuum-removed values are below a threshold (e.g., half-depth)
            double halfDepthThreshold = 1.0 - depth * 0.5;
            double leftWav = wavelengths[centerBand];
            double rightWav = wavelengths[centerBand];

            // Search left
            for (int b = centerBand - 1; b >= 0; --b) {
                double cr = (std::abs(continuum[b]) > 1e-15)
                            ? spectrum[b] / continuum[b] : 1.0;
                if (cr >= halfDepthThreshold) {
                    // Interpolate the crossing point
                    double crNext = (std::abs(continuum[b + 1]) > 1e-15)
                                    ? spectrum[b + 1] / continuum[b + 1] : 1.0;
                    double frac = (halfDepthThreshold - cr) / (crNext - cr);
                    leftWav = wavelengths[b] + frac * (wavelengths[b + 1] - wavelengths[b]);
                    break;
                }
                leftWav = wavelengths[b];
            }

            // Search right
            for (int b = centerBand + 1; b < numBands; ++b) {
                double cr = (std::abs(continuum[b]) > 1e-15)
                            ? spectrum[b] / continuum[b] : 1.0;
                if (cr >= halfDepthThreshold) {
                    double crPrev = (std::abs(continuum[b - 1]) > 1e-15)
                                    ? spectrum[b - 1] / continuum[b - 1] : 1.0;
                    double frac = (halfDepthThreshold - crPrev) / (cr - crPrev);
                    rightWav = wavelengths[b - 1] + frac * (wavelengths[b] - wavelengths[b - 1]);
                    break;
                }
                rightWav = wavelengths[b];
            }

            double width = rightWav - leftWav;
            if (width < 0.0) width = 0.0;

            depthData[idx] = depth;
            widthData[idx] = width;

            if (idx % 1000000 == 0)
                reportProgress(0.05 + 0.85 * static_cast<double>(idx) / total);
        }

        // --------------------------------------------------------------------
        // 6. Write outputs
        // --------------------------------------------------------------------
        reportProgress(0.9, "Writing output rasters...");

        QString depthPath = parameter("output_depth").toString();
        if (!GdalIO::write(outputDepth, depthPath)) {
            setError("Failed to write depth output: " + depthPath);
            return false;
        }

        QString widthPath = parameter("output_width").toString();
        if (!GdalIO::write(outputWidth, widthPath)) {
            setError("Failed to write width output: " + widthPath);
            return false;
        }

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(HyperAbsorbModule)

} // namespace aplaceholder
