#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ScenGenModule : public Module {
public:
    QString name() const override { return "SCENGEN"; }
    QString description() const override {
        return "Climate scenario generation using the delta change method. Applies a "
               "delta (change factor) to a baseline raster. Supports additive "
               "(future = baseline + delta) and multiplicative (future = baseline * delta) "
               "methods. Delta can be a single value or a spatially varying raster.";
    }
    QString category() const override { return "Climate Adaptation"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("baseline", "Baseline climate raster",
                "Baseline (current) climate raster (e.g., temperature, precipitation)"),
            ParameterDef::file("delta", "Delta value or raster",
                "Path to a delta raster, or a single numeric value as a string. "
                "For additive: change in same units as baseline. "
                "For multiplicative: ratio (e.g., 1.1 = 10% increase)."),
            ParameterDef::output("output", "Output future scenario raster",
                "Raster representing the projected future climate variable"),
            ParameterDef::combo("method", "Method",
                {"additive", "multiplicative"}, 0,
                "additive: future = baseline + delta; "
                "multiplicative: future = baseline * delta"),
        };
    }

    bool execute() override {
        // 1. Load baseline raster
        auto baseline = GdalIO::read(parameter("baseline").toString());
        if (!baseline) {
            setError("Failed to read baseline raster.");
            return false;
        }

        reportProgress(0.05, "Baseline raster loaded.");

        int cols = baseline->cols();
        int rows = baseline->rows();
        int64_t total = baseline->cellCount();

        // 2. Determine if delta is a raster or a scalar
        QString deltaStr = parameter("delta").toString();
        bool deltaIsScalar = false;
        double deltaScalar = 0.0;
        std::unique_ptr<Raster> deltaRaster;

        // Try to parse as a number first
        bool ok = false;
        deltaScalar = deltaStr.toDouble(&ok);
        if (ok) {
            deltaIsScalar = true;
            reportProgress(0.10, QString("Delta is scalar: %1").arg(deltaScalar));
        } else {
            // Try to load as a raster
            deltaRaster = GdalIO::read(deltaStr);
            if (!deltaRaster) {
                setError("Delta parameter is neither a valid number nor a readable raster: " + deltaStr);
                return false;
            }

            if (deltaRaster->cols() != cols || deltaRaster->rows() != rows) {
                setError(QString("Dimension mismatch: baseline is %1x%2 but delta raster is %3x%4")
                             .arg(cols).arg(rows)
                             .arg(deltaRaster->cols()).arg(deltaRaster->rows()));
                return false;
            }

            reportProgress(0.10, "Delta raster loaded.");
        }

        // 3. Determine method
        int methodIdx = parameter("method").toInt();
        bool isAdditive = (methodIdx == 0);

        // 4. NoData handling
        double noData = baseline->noDataValue();
        bool hasND_B = baseline->hasNoData();
        double ndB = baseline->noDataValue();

        bool hasND_D = false;
        double ndD = 0.0;
        const std::vector<double>* dData = nullptr;
        if (!deltaIsScalar) {
            hasND_D = deltaRaster->hasNoData();
            ndD = deltaRaster->noDataValue();
            dData = &deltaRaster->data(0);
        }

        // 5. Compute future = baseline op delta
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(baseline->geoTransform());
        output.setProjection(baseline->projection());
        output.setNoDataValue(noData);

        const auto& bData = baseline->data(0);
        auto& outData = output.data(0);

        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            if (hasND_B && bData[i] == ndB) skip = true;
            if (!deltaIsScalar && hasND_D && (*dData)[i] == ndD) skip = true;

            if (skip) {
                outData[i] = noData;
                continue;
            }

            double baseVal = bData[i];
            double deltaVal = deltaIsScalar ? deltaScalar : (*dData)[i];

            if (isAdditive) {
                outData[i] = baseVal + deltaVal;
            } else {
                outData[i] = baseVal * deltaVal;
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

        auto stats = output.computeStats(0);
        QString methodName = isAdditive ? "additive" : "multiplicative";
        reportProgress(1.0,
            QString("Done. Method: %1, output — min: %2, max: %3, mean: %4, valid pixels: %5")
                .arg(methodName)
                .arg(stats.min, 0, 'f', 4)
                .arg(stats.max, 0, 'f', 4)
                .arg(stats.mean, 0, 'f', 4)
                .arg(validCount));

        return true;
    }
};

REGISTER_MODULE(ScenGenModule)

} // namespace aplaceholder
