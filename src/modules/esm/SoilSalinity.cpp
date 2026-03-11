#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class SoilSalinityModule : public Module {
public:
    QString name() const override { return "SOIL_SALINITY"; }
    QString description() const override {
        return "Soil salinity index from satellite bands. Supports three index types: "
               "NDSI = (B1 - B2)/(B1 + B2), SI1 = sqrt(B1 * B2), "
               "SI2 = (B1 - B2)/(B1 + B2). Band assignment depends on sensor — "
               "typically band1=Green/Red and band2=Red/NIR.";
    }
    QString category() const override { return "Ecosystem Services"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("band1", "Band 1 raster",
                "First spectral band raster (e.g., Green or Red band)"),
            ParameterDef::file("band2", "Band 2 raster",
                "Second spectral band raster (e.g., Red or NIR band)"),
            ParameterDef::output("output", "Output salinity index raster",
                "Continuous raster of computed salinity index values"),
            ParameterDef::combo("index_type", "Index type",
                {"ndsi", "si1", "si2"}, 0,
                "ndsi: Normalized Difference Salinity Index = (B1-B2)/(B1+B2); "
                "si1: Salinity Index 1 = sqrt(B1*B2); "
                "si2: Salinity Index 2 = (B1-B2)/(B1+B2)"),
        };
    }

    bool execute() override {
        // 1. Load input bands
        auto band1 = GdalIO::read(parameter("band1").toString());
        if (!band1) {
            setError("Failed to read band 1 raster.");
            return false;
        }

        auto band2 = GdalIO::read(parameter("band2").toString());
        if (!band2) {
            setError("Failed to read band 2 raster.");
            return false;
        }

        reportProgress(0.10, "Input bands loaded.");

        // 2. Validate dimensions
        int cols = band1->cols();
        int rows = band1->rows();
        int64_t total = band1->cellCount();

        if (band2->cols() != cols || band2->rows() != rows) {
            setError(QString("Dimension mismatch: band1 is %1x%2 but band2 is %3x%4")
                         .arg(cols).arg(rows)
                         .arg(band2->cols()).arg(band2->rows()));
            return false;
        }

        // 3. Determine index type
        int indexType = parameter("index_type").toInt();
        // 0 = ndsi, 1 = si1, 2 = si2

        // 4. NoData handling
        double noData = band1->noDataValue();
        bool hasND_1 = band1->hasNoData();
        bool hasND_2 = band2->hasNoData();
        double nd1 = band1->noDataValue();
        double nd2 = band2->noDataValue();

        // 5. Compute salinity index
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(band1->geoTransform());
        output.setProjection(band1->projection());
        output.setNoDataValue(noData);

        const auto& b1Data = band1->data(0);
        const auto& b2Data = band2->data(0);
        auto& outData = output.data(0);

        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            if (hasND_1 && b1Data[i] == nd1) skip = true;
            if (hasND_2 && b2Data[i] == nd2) skip = true;

            if (skip) {
                outData[i] = noData;
                continue;
            }

            double b1 = b1Data[i];
            double b2 = b2Data[i];

            switch (indexType) {
            case 0: // NDSI: (B1 - B2) / (B1 + B2)
            case 2: // SI2:  same formula as NDSI
            {
                double sum = b1 + b2;
                if (std::abs(sum) < 1e-10) {
                    outData[i] = 0.0;
                } else {
                    outData[i] = (b1 - b2) / sum;
                }
                break;
            }
            case 1: // SI1: sqrt(B1 * B2)
            {
                double product = b1 * b2;
                if (product < 0.0) {
                    outData[i] = 0.0;
                } else {
                    outData[i] = std::sqrt(product);
                }
                break;
            }
            default:
                outData[i] = noData;
                break;
            }

            ++validCount;

            if (i % 500000 == 0)
                reportProgress(0.10 + 0.80 * static_cast<double>(i) / total);
        }

        reportProgress(0.90, "Writing output...");

        // 6. Write output
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        QStringList indexNames = {"NDSI", "SI1", "SI2"};
        auto stats = output.computeStats(0);
        reportProgress(1.0,
            QString("Done. %1 — min: %2, max: %3, mean: %4, valid pixels: %5")
                .arg(indexNames.value(indexType, "Unknown"))
                .arg(stats.min, 0, 'f', 4)
                .arg(stats.max, 0, 'f', 4)
                .arg(stats.mean, 0, 'f', 4)
                .arg(validCount));

        return true;
    }
};

REGISTER_MODULE(SoilSalinityModule)

} // namespace aplaceholder
