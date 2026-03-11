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

class HyperSamModule : public Module {
public:
    QString name() const override { return "HYPERSAM"; }
    QString description() const override {
        return "Spectral Angle Mapper for hyperspectral image classification. "
               "Computes the angle between each pixel spectrum and reference endmember spectra. "
               "Assigns each pixel to the endmember with the smallest spectral angle.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::file("endmember_file", "Endmember file (CSV)",
                "CSV: each row is an endmember, columns are band values"),
            ParameterDef::output("output", "Output classified image",
                "Output raster with class assignments (0 = unclassified)"),
            ParameterDef::real("threshold_angle", "Threshold angle (radians)",
                0.1, 0.0, 3.14159,
                "Maximum spectral angle for classification; pixels exceeding this are unclassified"),
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
        // 2. Read endmember file (CSV)
        // --------------------------------------------------------------------
        QString emPath = parameter("endmember_file").toString();
        std::vector<std::vector<double>> endmembers;

        {
            std::ifstream file(emPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open endmember file: " + emPath);
                return false;
            }

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;

                std::istringstream iss(line);
                std::vector<double> values;
                std::string token;
                while (std::getline(iss, token, ',')) {
                    size_t start = token.find_first_not_of(" \t");
                    size_t end = token.find_last_not_of(" \t\r\n");
                    if (start != std::string::npos && end != std::string::npos)
                        values.push_back(std::stod(token.substr(start, end - start + 1)));
                }

                if (static_cast<int>(values.size()) >= numBands) {
                    values.resize(numBands);
                    endmembers.push_back(std::move(values));
                }
            }
        }

        int numEM = static_cast<int>(endmembers.size());
        if (numEM == 0) {
            setError("No valid endmembers found in endmember file");
            return false;
        }

        double thresholdAngle = parameter("threshold_angle").toDouble();

        // --------------------------------------------------------------------
        // 3. Precompute endmember magnitudes
        // --------------------------------------------------------------------
        std::vector<double> emMagnitude(numEM, 0.0);
        for (int e = 0; e < numEM; ++e) {
            double sum = 0.0;
            for (int b = 0; b < numBands; ++b)
                sum += endmembers[e][b] * endmembers[e][b];
            emMagnitude[e] = std::sqrt(sum);
        }

        reportProgress(0.05, "Classifying pixels via Spectral Angle Mapper...");

        // --------------------------------------------------------------------
        // 4. Collect band data pointers
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // --------------------------------------------------------------------
        // 5. Create output raster (class assignment)
        // --------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Float32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);

        // --------------------------------------------------------------------
        // 6. For each pixel, compute spectral angle to each endmember
        // --------------------------------------------------------------------
        std::vector<double> pixelVals(numBands);

        for (int64_t idx = 0; idx < total; ++idx) {
            if (hasND && (*bands[0])[idx] == noData) {
                outData[idx] = noData;
                continue;
            }

            // Get pixel spectral values and magnitude
            double pixMag = 0.0;
            for (int b = 0; b < numBands; ++b) {
                pixelVals[b] = (*bands[b])[idx];
                pixMag += pixelVals[b] * pixelVals[b];
            }
            pixMag = std::sqrt(pixMag);

            if (pixMag < 1e-12) {
                outData[idx] = 0.0; // unclassified
                continue;
            }

            // Find endmember with minimum spectral angle
            double minAngle = std::numeric_limits<double>::max();
            int bestClass = 0; // 0 = unclassified

            for (int e = 0; e < numEM; ++e) {
                if (emMagnitude[e] < 1e-12) continue;

                double dotProduct = 0.0;
                for (int b = 0; b < numBands; ++b)
                    dotProduct += pixelVals[b] * endmembers[e][b];

                double cosAngle = dotProduct / (pixMag * emMagnitude[e]);
                // Clamp to valid range for acos
                cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                double angle = std::acos(cosAngle);

                if (angle < minAngle) {
                    minAngle = angle;
                    bestClass = e + 1; // 1-based class assignment
                }
            }

            // Apply threshold: if minimum angle exceeds threshold, leave unclassified
            if (minAngle > thresholdAngle)
                bestClass = 0;

            outData[idx] = static_cast<double>(bestClass);

            if (idx % 1000000 == 0)
                reportProgress(0.05 + 0.85 * static_cast<double>(idx) / total);
        }

        // --------------------------------------------------------------------
        // 7. Write output
        // --------------------------------------------------------------------
        reportProgress(0.9, "Writing output...");
        QString outputPath = parameter("output").toString();
        if (!GdalIO::write(output, outputPath)) {
            setError("Failed to write output: " + outputPath);
            return false;
        }

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(HyperSamModule)

} // namespace aplaceholder
