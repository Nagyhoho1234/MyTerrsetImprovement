#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <QStringList>

namespace aplaceholder {

class AtmosCModule : public Module {
public:
    QString name() const override { return "ATMOSC"; }
    QString description() const override {
        return "Atmospheric correction of satellite imagery. Supports Dark Object "
               "Subtraction (DOS1) and Cos(t) (COST) correction methods. DOS subtracts "
               "the minimum DN per band; COST additionally applies cos(theta) correction "
               "for solar zenith angle.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input bands (comma-separated file paths)"),
            ParameterDef::output("output_prefix", "Output file prefix"),
            ParameterDef::combo("method", "Correction method",
                                {"DOS1", "COST"}, 0,
                                "DOS1: Dark Object Subtraction. "
                                "COST: Cos(t) model with solar zenith correction."),
            ParameterDef::real("sun_elevation", "Sun elevation (degrees)", 45.0, 0.0, 90.0,
                               "Required for COST method. Solar elevation angle in degrees."),
        };
    }

    bool execute() override {
        QString inputBandsStr = parameter("input_bands").toString();
        QString outputPrefix = parameter("output_prefix").toString();
        int methodIdx = parameter("method").toInt();
        double sunElevation = parameter("sun_elevation").toDouble();

        QStringList bandPaths = inputBandsStr.split(',', Qt::SkipEmptyParts);
        for (auto& p : bandPaths)
            p = p.trimmed();

        if (bandPaths.isEmpty()) {
            setError("No input bands specified");
            return false;
        }

        // Read all input bands
        std::vector<std::unique_ptr<Raster>> bands;
        for (const auto& path : bandPaths) {
            auto raster = GdalIO::read(path);
            if (!raster) {
                setError("Failed to read input band: " + path);
                return false;
            }
            bands.push_back(std::move(raster));
        }

        // Validate dimensions
        int cols = bands[0]->cols();
        int rows = bands[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (size_t b = 1; b < bands.size(); ++b) {
            if (bands[b]->cols() != cols || bands[b]->rows() != rows) {
                setError("All input bands must have the same dimensions");
                return false;
            }
        }

        // Compute cos(theta) for COST method
        // theta = solar zenith angle = 90 - sun_elevation
        double cosTheta = 1.0;
        if (methodIdx == 1) { // COST
            double zenithDeg = 90.0 - sunElevation;
            double zenithRad = zenithDeg * M_PI / 180.0;
            cosTheta = std::cos(zenithRad);
            if (cosTheta <= 0.0) {
                setError("Invalid sun elevation: cos(zenith) <= 0");
                return false;
            }
        }

        int totalBands = static_cast<int>(bands.size());
        int64_t totalWork = total * totalBands;
        int64_t processed = 0;

        for (int b = 0; b < totalBands; ++b) {
            const auto& inData = bands[b]->data(0);
            bool hasND = bands[b]->hasNoData();
            double inputND = bands[b]->noDataValue();

            // Find the minimum DN (dark object value), excluding NoData
            double minDN = std::numeric_limits<double>::max();
            for (int64_t i = 0; i < total; ++i) {
                double val = inData[i];
                if (hasND && val == inputND)
                    continue;
                if (val < minDN)
                    minDN = val;
            }

            if (minDN == std::numeric_limits<double>::max())
                minDN = 0.0; // all NoData -- fallback

            // Create output raster
            Raster output(cols, rows, 1, DataType::Float64);
            output.setGeoTransform(bands[b]->geoTransform());
            output.setProjection(bands[b]->projection());
            double noData = -9999.0;
            output.setNoDataValue(noData);
            auto& outData = output.data(0);

            for (int64_t i = 0; i < total; ++i) {
                double val = inData[i];
                if (hasND && val == inputND) {
                    outData[i] = noData;
                } else {
                    // DOS: subtract dark object value
                    double corrected = val - minDN;
                    if (corrected < 0.0)
                        corrected = 0.0;

                    // COST: divide by cos(theta) to correct for solar zenith
                    if (methodIdx == 1) {
                        corrected = corrected / cosTheta;
                    }

                    outData[i] = corrected;
                }

                ++processed;
                if (processed % 1000000 == 0)
                    reportProgress(static_cast<double>(processed) / totalWork);
            }

            // Write output with band suffix
            QString outPath = outputPrefix + "_band" + QString::number(b + 1);
            if (!GdalIO::write(output, outPath)) {
                setError("Failed to write output for band " + QString::number(b + 1));
                return false;
            }
        }

        reportProgress(1.0, "Atmospheric correction complete.");
        return true;
    }
};

REGISTER_MODULE(AtmosCModule)

} // namespace aplaceholder
