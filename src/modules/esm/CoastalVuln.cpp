#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <memory>

namespace aplaceholder {

class CoastalVulnModule : public Module {
public:
    QString name() const override { return "COASTAL_VULNERABILITY"; }
    QString description() const override {
        return "Coastal Vulnerability Index (CVI). Combines six factor rasters — "
               "geomorphology, coastal slope, sea level rise, wave height, tidal range, "
               "and shoreline change — using the formula: CVI = sqrt((a*b*c*d*e*f)/6).";
    }
    QString category() const override { return "Ecosystem Services"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("factors", "Factor rasters (6 comma-separated paths)",
                "Six comma-separated raster paths in order: geomorphology, coastal slope, "
                "sea level rise, wave height, tidal range, shoreline change. "
                "Each should be ranked 1-5 (very low to very high vulnerability)."),
            ParameterDef::output("output", "Output CVI raster",
                "Continuous raster of Coastal Vulnerability Index values"),
        };
    }

    bool execute() override {
        // 1. Parse comma-separated factor paths
        QString factorsStr = parameter("factors").toString();
        QStringList factorPaths = factorsStr.split(',', Qt::SkipEmptyParts);

        for (int i = 0; i < factorPaths.size(); ++i)
            factorPaths[i] = factorPaths[i].trimmed();

        if (factorPaths.size() != 6) {
            setError(QString("Expected exactly 6 factor rasters, got %1. "
                             "Provide paths separated by commas.")
                         .arg(factorPaths.size()));
            return false;
        }

        // 2. Load all six factor rasters
        std::vector<std::unique_ptr<Raster>> factors;
        QStringList factorNames = {"Geomorphology", "Coastal Slope", "Sea Level Rise",
                                   "Wave Height", "Tidal Range", "Shoreline Change"};

        for (int i = 0; i < 6; ++i) {
            auto raster = GdalIO::read(factorPaths[i]);
            if (!raster) {
                setError("Failed to read " + factorNames[i] + " raster: " + factorPaths[i]);
                return false;
            }
            factors.push_back(std::move(raster));
        }

        reportProgress(0.10, "All 6 factor rasters loaded.");

        // 3. Validate dimensions
        int cols = factors[0]->cols();
        int rows = factors[0]->rows();
        int64_t total = factors[0]->cellCount();

        for (int i = 1; i < 6; ++i) {
            if (factors[i]->cols() != cols || factors[i]->rows() != rows) {
                setError(QString("Dimension mismatch: %1 is %2x%3 but expected %4x%5")
                             .arg(factorNames[i])
                             .arg(factors[i]->cols()).arg(factors[i]->rows())
                             .arg(cols).arg(rows));
                return false;
            }
        }

        // 4. Compute CVI = sqrt((a*b*c*d*e*f) / 6)
        double noData = factors[0]->noDataValue();
        bool hasND = factors[0]->hasNoData();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(factors[0]->geoTransform());
        output.setProjection(factors[0]->projection());
        output.setNoDataValue(noData);

        auto& outData = output.data(0);

        // Get references to all factor data arrays
        std::vector<const std::vector<double>*> fData;
        for (int i = 0; i < 6; ++i)
            fData.push_back(&factors[i]->data(0));

        std::vector<bool> fHasND(6);
        std::vector<double> fND(6);
        for (int i = 0; i < 6; ++i) {
            fHasND[i] = factors[i]->hasNoData();
            fND[i] = factors[i]->noDataValue();
        }

        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            for (int f = 0; f < 6; ++f) {
                if (fHasND[f] && (*fData[f])[i] == fND[f]) {
                    skip = true;
                    break;
                }
            }

            if (skip) {
                outData[i] = noData;
                continue;
            }

            double product = 1.0;
            for (int f = 0; f < 6; ++f)
                product *= (*fData[f])[i];

            outData[i] = std::sqrt(product / 6.0);
            ++validCount;

            if (i % 500000 == 0)
                reportProgress(0.10 + 0.80 * static_cast<double>(i) / total);
        }

        reportProgress(0.90, "Writing output...");

        // 5. Write output
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        auto stats = output.computeStats(0);
        reportProgress(1.0,
            QString("Done. CVI — min: %1, max: %2, mean: %3, valid pixels: %4")
                .arg(stats.min, 0, 'f', 4)
                .arg(stats.max, 0, 'f', 4)
                .arg(stats.mean, 0, 'f', 4)
                .arg(validCount));

        return true;
    }
};

REGISTER_MODULE(CoastalVulnModule)

} // namespace aplaceholder
