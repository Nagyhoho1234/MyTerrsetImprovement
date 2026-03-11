#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <QFile>
#include <QTextStream>

namespace aplaceholder {

class RadianceModule : public Module {
public:
    QString name() const override { return "RADIANCE"; }
    QString description() const override {
        return "Convert raw Digital Number (DN) values to at-sensor spectral radiance. "
               "Uses L = gain * DN + offset, or reads gain/offset per band from a "
               "calibration file.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image (DN values)"),
            ParameterDef::output("output", "Output radiance image"),
            ParameterDef::real("gain", "Gain (multiplicative factor)", 1.0, -999999, 999999,
                               "Radiance gain. Ignored if calibration file is provided."),
            ParameterDef::real("offset", "Offset (additive constant)", 0.0, -999999, 999999,
                               "Radiance offset. Ignored if calibration file is provided."),
            ParameterDef::file("calibration_file", "Calibration file (optional)",
                               "Text file with gain and offset per band, one line per band: "
                               "gain offset"),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString outputPath = parameter("output").toString();
        QString calFile = parameter("calibration_file").toString();

        auto input = GdalIO::read(inputPath);
        if (!input) {
            setError("Failed to read input image: " + inputPath);
            return false;
        }

        int numBands = input->bands();
        int cols = input->cols();
        int rows = input->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Build per-band gain/offset arrays
        std::vector<double> gains(numBands, parameter("gain").toDouble());
        std::vector<double> offsets(numBands, parameter("offset").toDouble());

        // If a calibration file is provided, read gain/offset per band
        if (!calFile.isEmpty()) {
            QFile file(calFile);
            if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
                setError("Failed to open calibration file: " + calFile);
                return false;
            }
            QTextStream in(&file);
            int band = 0;
            while (!in.atEnd() && band < numBands) {
                QString line = in.readLine().trimmed();
                if (line.isEmpty() || line.startsWith('#'))
                    continue;
                QStringList parts = line.simplified().split(' ');
                if (parts.size() < 2) {
                    setError("Calibration file line " + QString::number(band + 1) +
                             " must have at least gain and offset values");
                    return false;
                }
                bool okG = false, okO = false;
                gains[band] = parts[0].toDouble(&okG);
                offsets[band] = parts[1].toDouble(&okO);
                if (!okG || !okO) {
                    setError("Invalid numeric values in calibration file at band " +
                             QString::number(band + 1));
                    return false;
                }
                ++band;
            }
            if (band < numBands) {
                setError("Calibration file has fewer entries (" + QString::number(band) +
                         ") than image bands (" + QString::number(numBands) + ")");
                return false;
            }
        }

        // Create output raster
        Raster output(cols, rows, numBands, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        double noData = -9999.0;
        output.setNoDataValue(noData);

        bool hasND = input->hasNoData();
        double inputND = input->noDataValue();

        int64_t totalPixels = total * numBands;
        int64_t processed = 0;

        for (int b = 0; b < numBands; ++b) {
            double g = gains[b];
            double o = offsets[b];
            const auto& inData = input->data(b);
            auto& outData = output.data(b);

            for (int64_t i = 0; i < total; ++i) {
                double dn = inData[i];
                if (hasND && dn == inputND) {
                    outData[i] = noData;
                } else {
                    // L = gain * DN + offset
                    outData[i] = g * dn + o;
                }

                ++processed;
                if (processed % 1000000 == 0)
                    reportProgress(static_cast<double>(processed) / totalPixels);
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
    }
};

REGISTER_MODULE(RadianceModule)

} // namespace aplaceholder
