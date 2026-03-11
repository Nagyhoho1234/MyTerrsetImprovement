#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>

namespace aplaceholder {

class CalibrateModule : public Module {
public:
    QString name() const override { return "CALIBRATE"; }
    QString description() const override {
        return "DN to Radiance Calibration. Applies per-band offset and gain "
               "calibration coefficients from a metadata file to convert raw "
               "Digital Number values to calibrated radiance or reflectance. "
               "Formula: output = Offset + (Gain * DN).";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated paths)",
                "Comma-separated list of band raster file paths to calibrate"),
            ParameterDef::file("calibration_file", "Calibration coefficients file",
                "Text file with per-band offset and gain. Each non-comment line: "
                "band_index, offset, gain (CSV format). Lines starting with # are comments."),
            ParameterDef::output("output_prefix", "Output prefix for calibrated images"),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Parse band file paths
        // ------------------------------------------------------------------
        QString bandsParam = parameter("input_bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths)
            p = p.trimmed();

        if (bandPaths.isEmpty()) {
            setError("No band images specified");
            return false;
        }

        int numBands = bandPaths.size();

        // ------------------------------------------------------------------
        // 2. Read calibration coefficients
        // ------------------------------------------------------------------
        QString calPath = parameter("calibration_file").toString();

        struct CalCoeff {
            double offset;
            double gain;
        };

        std::vector<CalCoeff> coeffs(numBands, {0.0, 1.0}); // defaults

        {
            std::ifstream file(calPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open calibration file: " + calPath);
                return false;
            }

            std::string line;
            int lineIdx = 0;
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

                if (values.size() >= 3) {
                    // Format: band_index, offset, gain
                    int bandIdx = static_cast<int>(values[0]);
                    if (bandIdx >= 0 && bandIdx < numBands) {
                        coeffs[bandIdx].offset = values[1];
                        coeffs[bandIdx].gain = values[2];
                    }
                } else if (values.size() >= 2) {
                    // Format: offset, gain (sequential by line order)
                    if (lineIdx < numBands) {
                        coeffs[lineIdx].offset = values[0];
                        coeffs[lineIdx].gain = values[1];
                    }
                }

                lineIdx++;
            }
        }

        // ------------------------------------------------------------------
        // 3. Process each band
        // ------------------------------------------------------------------
        QString prefix = parameter("output_prefix").toString();

        for (int b = 0; b < numBands; ++b) {
            reportProgress(static_cast<double>(b) / numBands,
                           QString("Calibrating band %1 of %2...").arg(b + 1).arg(numBands));

            auto raster = GdalIO::read(bandPaths[b]);
            if (!raster) {
                setError("Failed to read band image: " + bandPaths[b]);
                return false;
            }

            int cols = raster->cols(), rows = raster->rows();
            int64_t total = static_cast<int64_t>(cols) * rows;
            bool hasND = raster->hasNoData();
            double nd = raster->noDataValue();

            // Create output raster
            Raster output(cols, rows, 1, DataType::Float32);
            output.setGeoTransform(raster->geoTransform());
            output.setProjection(raster->projection());
            if (hasND) output.setNoDataValue(nd);

            const auto& inData = raster->data(0);
            auto& outData = output.data(0);

            double offset = coeffs[b].offset;
            double gain = coeffs[b].gain;

            for (int64_t i = 0; i < total; ++i) {
                if (hasND && inData[i] == nd) {
                    outData[i] = nd;
                    continue;
                }
                // output = Offset + (Gain * DN)
                outData[i] = offset + gain * inData[i];
            }

            // Write output
            QString outPath = QString("%1_band%2.tif").arg(prefix).arg(b + 1);
            if (!GdalIO::write(output, outPath)) {
                setError("Failed to write calibrated image: " + outPath);
                return false;
            }
        }

        reportProgress(1.0, QString("Calibration complete: %1 bands processed.").arg(numBands));
        return true;
    }
};

REGISTER_MODULE(CalibrateModule)

} // namespace aplaceholder
