#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class EcoCropModule : public Module {
public:
    QString name() const override { return "ECOCROP"; }
    QString description() const override {
        return "Ecological crop suitability assessment. Evaluates temperature and "
               "precipitation suitability for a given crop based on parameter ranges "
               "(Tmin, Topt_low, Topt_high, Tmax, Pmin, Popt_low, Popt_high, Pmax). "
               "Suitability is 0 outside absolute limits, linearly interpolated in "
               "marginal zones, and 1 in the optimal range. Final suitability is the "
               "minimum of temperature and precipitation suitability.";
    }
    QString category() const override { return "Climate Adaptation"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("temperature", "Temperature raster (degrees C)",
                "Mean annual temperature or growing season temperature in degrees Celsius"),
            ParameterDef::file("precipitation", "Precipitation raster (mm)",
                "Mean annual precipitation or growing season precipitation in mm"),
            ParameterDef::output("output", "Output crop suitability raster",
                "Continuous raster [0-1] where 1 = fully suitable"),
            ParameterDef::file("crop_params_file", "Crop parameters file (CSV)",
                "CSV with one data row containing: Tmin, Topt_low, Topt_high, Tmax, "
                "Pmin, Popt_low, Popt_high, Pmax"),
        };
    }

    bool execute() override {
        // 1. Load input rasters
        auto tempRaster = GdalIO::read(parameter("temperature").toString());
        if (!tempRaster) {
            setError("Failed to read temperature raster.");
            return false;
        }

        auto precRaster = GdalIO::read(parameter("precipitation").toString());
        if (!precRaster) {
            setError("Failed to read precipitation raster.");
            return false;
        }

        reportProgress(0.05, "Input rasters loaded.");

        // 2. Validate dimensions
        int cols = tempRaster->cols();
        int rows = tempRaster->rows();
        int64_t total = tempRaster->cellCount();

        if (precRaster->cols() != cols || precRaster->rows() != rows) {
            setError(QString("Dimension mismatch: temperature is %1x%2 but precipitation is %3x%4")
                         .arg(cols).arg(rows)
                         .arg(precRaster->cols()).arg(precRaster->rows()));
            return false;
        }

        // 3. Parse crop parameters CSV
        //    Expected: header row, then one data row with 8 values
        //    Tmin, Topt_low, Topt_high, Tmax, Pmin, Popt_low, Popt_high, Pmax
        QString csvPath = parameter("crop_params_file").toString();
        double Tmin, ToptLow, ToptHigh, Tmax;
        double Pmin, PoptLow, PoptHigh, Pmax;

        {
            QFile csvFile(csvPath);
            if (!csvFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
                setError("Failed to open crop parameters file: " + csvPath);
                return false;
            }

            QTextStream in(&csvFile);
            bool headerSkipped = false;
            bool parsed = false;

            while (!in.atEnd()) {
                QString line = in.readLine().trimmed();
                if (line.isEmpty() || line.startsWith('#'))
                    continue;

                if (!headerSkipped) {
                    headerSkipped = true;
                    continue;
                }

                QStringList fields = line.split(',');
                if (fields.size() < 8) {
                    setError("Crop parameters CSV: expected 8 values "
                             "(Tmin, Topt_low, Topt_high, Tmax, Pmin, Popt_low, Popt_high, Pmax)");
                    return false;
                }

                Tmin     = fields[0].trimmed().toDouble();
                ToptLow  = fields[1].trimmed().toDouble();
                ToptHigh = fields[2].trimmed().toDouble();
                Tmax     = fields[3].trimmed().toDouble();
                Pmin     = fields[4].trimmed().toDouble();
                PoptLow  = fields[5].trimmed().toDouble();
                PoptHigh = fields[6].trimmed().toDouble();
                Pmax     = fields[7].trimmed().toDouble();
                parsed = true;
                break;
            }

            if (!parsed) {
                setError("Crop parameters CSV: no data row found.");
                return false;
            }
        }

        reportProgress(0.10,
            QString("Crop params — T:[%1,%2,%3,%4] P:[%5,%6,%7,%8]")
                .arg(Tmin).arg(ToptLow).arg(ToptHigh).arg(Tmax)
                .arg(Pmin).arg(PoptLow).arg(PoptHigh).arg(Pmax));

        // 4. NoData handling
        double noData = tempRaster->noDataValue();
        bool hasND_T = tempRaster->hasNoData();
        bool hasND_P = precRaster->hasNoData();
        double ndT = tempRaster->noDataValue();
        double ndP = precRaster->noDataValue();

        // 5. Compute suitability
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(tempRaster->geoTransform());
        output.setProjection(tempRaster->projection());
        output.setNoDataValue(noData);

        const auto& tData = tempRaster->data(0);
        const auto& pData = precRaster->data(0);
        auto& outData = output.data(0);

        int64_t validCount = 0;
        double sumSuit = 0.0;

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            if (hasND_T && tData[i] == ndT) skip = true;
            if (hasND_P && pData[i] == ndP) skip = true;

            if (skip) {
                outData[i] = noData;
                continue;
            }

            double temp = tData[i];
            double prec = pData[i];

            // Temperature suitability (trapezoidal function)
            double tSuit = trapezoidalSuitability(temp, Tmin, ToptLow, ToptHigh, Tmax);

            // Precipitation suitability (trapezoidal function)
            double pSuit = trapezoidalSuitability(prec, Pmin, PoptLow, PoptHigh, Pmax);

            // Overall suitability is the minimum (limiting factor approach)
            double suit = std::min(tSuit, pSuit);

            outData[i] = suit;
            sumSuit += suit;
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

        double meanSuit = (validCount > 0) ? sumSuit / validCount : 0.0;
        reportProgress(1.0,
            QString("Done. Mean suitability: %1, valid pixels: %2")
                .arg(meanSuit, 0, 'f', 4)
                .arg(validCount));

        return true;
    }

private:
    // Trapezoidal suitability function:
    //   0 below absMin or above absMax
    //   linear ramp from absMin to optLow and from optHigh to absMax
    //   1 between optLow and optHigh
    static double trapezoidalSuitability(double val, double absMin, double optLow,
                                          double optHigh, double absMax) {
        if (val <= absMin || val >= absMax)
            return 0.0;
        if (val >= optLow && val <= optHigh)
            return 1.0;
        if (val < optLow) {
            // Linear ramp up from absMin to optLow
            double range = optLow - absMin;
            return (range > 0.0) ? (val - absMin) / range : 0.0;
        }
        // val > optHigh: linear ramp down from optHigh to absMax
        double range = absMax - optHigh;
        return (range > 0.0) ? (absMax - val) / range : 0.0;
    }
};

REGISTER_MODULE(EcoCropModule)

} // namespace aplaceholder
