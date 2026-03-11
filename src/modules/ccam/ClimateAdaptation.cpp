#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace aplaceholder {

class ClimateAdaptationModule : public Module {
public:
    QString name() const override { return "CLIMATE_ADAPTATION"; }
    QString description() const override {
        return "Climate change vulnerability and adaptation assessment. Computes a "
               "composite Vulnerability Index from three input layers — exposure, "
               "sensitivity, and adaptive capacity — using the formula: "
               "V = (E_norm + S_norm) / 2 - weight * AC_norm, clamped to [0, 1]. "
               "Each input is min-max normalized to [0, 1] before combination.";
    }
    QString category() const override { return "Climate Adaptation"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("exposure", "Exposure raster",
                "Continuous raster representing climate exposure "
                "(e.g., projected temperature change in degrees)"),
            ParameterDef::file("sensitivity", "Sensitivity raster",
                "Continuous raster representing system sensitivity "
                "(e.g., crop-dependence index, population density)"),
            ParameterDef::file("adaptive_capacity", "Adaptive capacity raster",
                "Continuous raster representing adaptive capacity "
                "(e.g., infrastructure index, GDP per capita index)"),
            ParameterDef::output("output", "Output vulnerability raster",
                "Continuous raster [0-1] where 1 = highest vulnerability"),
            ParameterDef::real("ac_weight", "Adaptive capacity weight",
                0.5, 0.0, 2.0,
                "Weight applied to the normalized adaptive capacity layer "
                "when computing vulnerability. Default 0.5. Higher values "
                "give more credit to adaptive capacity, reducing vulnerability."),
        };
    }

    bool execute() override {
        // ----------------------------------------------------------
        // 1. Load the three input rasters
        // ----------------------------------------------------------
        auto exposure    = GdalIO::read(parameter("exposure").toString());
        auto sensitivity = GdalIO::read(parameter("sensitivity").toString());
        auto adapCap     = GdalIO::read(parameter("adaptive_capacity").toString());

        if (!exposure) {
            setError("Failed to read exposure raster.");
            return false;
        }
        if (!sensitivity) {
            setError("Failed to read sensitivity raster.");
            return false;
        }
        if (!adapCap) {
            setError("Failed to read adaptive capacity raster.");
            return false;
        }

        reportProgress(0.05, "Input rasters loaded.");

        // ----------------------------------------------------------
        // 2. Validate dimensions
        // ----------------------------------------------------------
        int cols = exposure->cols();
        int rows = exposure->rows();
        int64_t total = exposure->cellCount();

        if (sensitivity->cols() != cols || sensitivity->rows() != rows ||
            adapCap->cols() != cols    || adapCap->rows() != rows) {
            setError("All three input rasters must have the same dimensions. "
                     "Exposure: " + QString("%1x%2").arg(cols).arg(rows) +
                     ", Sensitivity: " + QString("%1x%2").arg(sensitivity->cols())
                                                           .arg(sensitivity->rows()) +
                     ", Adaptive Capacity: " + QString("%1x%2").arg(adapCap->cols())
                                                                 .arg(adapCap->rows()));
            return false;
        }

        // ----------------------------------------------------------
        // 3. Determine NoData handling
        //    Use exposure's NoData as the master value; treat any
        //    pixel that is NoData in ANY input as NoData in output.
        // ----------------------------------------------------------
        double noData = exposure->noDataValue();
        bool hasND_E  = exposure->hasNoData();
        bool hasND_S  = sensitivity->hasNoData();
        bool hasND_A  = adapCap->hasNoData();
        double ndE = exposure->noDataValue();
        double ndS = sensitivity->noDataValue();
        double ndA = adapCap->noDataValue();

        // ----------------------------------------------------------
        // 4. First pass — compute min/max of each layer for
        //    min-max normalization, skipping NoData.
        // ----------------------------------------------------------
        const auto& dE = exposure->data(0);
        const auto& dS = sensitivity->data(0);
        const auto& dA = adapCap->data(0);

        double minE =  std::numeric_limits<double>::max();
        double maxE = -std::numeric_limits<double>::max();
        double minS =  std::numeric_limits<double>::max();
        double maxS = -std::numeric_limits<double>::max();
        double minA =  std::numeric_limits<double>::max();
        double maxA = -std::numeric_limits<double>::max();

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            if (hasND_E && dE[i] == ndE) skip = true;
            if (hasND_S && dS[i] == ndS) skip = true;
            if (hasND_A && dA[i] == ndA) skip = true;
            if (skip) continue;

            if (dE[i] < minE) minE = dE[i];
            if (dE[i] > maxE) maxE = dE[i];
            if (dS[i] < minS) minS = dS[i];
            if (dS[i] > maxS) maxS = dS[i];
            if (dA[i] < minA) minA = dA[i];
            if (dA[i] > maxA) maxA = dA[i];

            if (i % 2000000 == 0)
                reportProgress(0.05 + 0.25 * static_cast<double>(i) / total,
                               "Computing statistics...");
        }

        // Guard against constant-valued layers (range = 0)
        double rangeE = (maxE > minE) ? (maxE - minE) : 1.0;
        double rangeS = (maxS > minS) ? (maxS - minS) : 1.0;
        double rangeA = (maxA > minA) ? (maxA - minA) : 1.0;

        reportProgress(0.30, QString("Stats — E:[%1,%2] S:[%3,%4] AC:[%5,%6]")
                        .arg(minE).arg(maxE).arg(minS).arg(maxS)
                        .arg(minA).arg(maxA));

        // ----------------------------------------------------------
        // 5. Second pass — normalise and compute vulnerability
        //    V = (E_n + S_n) / 2  -  weight * AC_n
        //    Clamp to [0, 1]
        // ----------------------------------------------------------
        double acWeight = parameter("ac_weight").toDouble();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(exposure->geoTransform());
        output.setProjection(exposure->projection());
        output.setNoDataValue(noData);

        auto& outData = output.data(0);

        double sumV = 0.0;
        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            if (hasND_E && dE[i] == ndE) skip = true;
            if (hasND_S && dS[i] == ndS) skip = true;
            if (hasND_A && dA[i] == ndA) skip = true;

            if (skip) {
                outData[i] = noData;
                continue;
            }

            // Min-max normalize to [0, 1]
            double eNorm = (dE[i] - minE) / rangeE;
            double sNorm = (dS[i] - minS) / rangeS;
            double aNorm = (dA[i] - minA) / rangeA;

            // Composite vulnerability index
            double v = (eNorm + sNorm) / 2.0 - acWeight * aNorm;

            // Clamp to [0, 1]
            v = std::max(0.0, std::min(1.0, v));

            outData[i] = v;
            sumV += v;
            ++validCount;

            if (i % 500000 == 0)
                reportProgress(0.30 + 0.60 * static_cast<double>(i) / total);
        }

        reportProgress(0.90, "Writing output...");

        // ----------------------------------------------------------
        // 6. Write output raster
        // ----------------------------------------------------------
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        double meanV = (validCount > 0) ? sumV / validCount : 0.0;
        reportProgress(1.0,
            QString("Done. Vulnerability index — mean: %1, valid pixels: %2")
                .arg(meanV, 0, 'f', 4)
                .arg(validCount));

        return true;
    }
};

REGISTER_MODULE(ClimateAdaptationModule)

} // namespace aplaceholder
