#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <QMap>
#include <cmath>
#include <limits>

namespace aplaceholder {

class ExtractModule : public Module {
public:
    QString name() const override { return "EXTRACT"; }
    QString description() const override {
        return "Extract raster values using a feature/mask raster. For each zone in the "
               "feature raster, computes statistics (mean/sum/min/max/count) from a value raster.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("feature_raster", "Feature/zone raster",
                "Raster with integer zone identifiers"),
            ParameterDef::file("value_raster", "Value raster",
                "Raster from which to extract statistics"),
            ParameterDef::output("output_file", "Output CSV file"),
            ParameterDef::combo("statistic", "Statistic",
                {"mean", "sum", "min", "max", "count"}, 0,
                "Statistic to compute per zone"),
        };
    }

    bool execute() override {
        auto features = GdalIO::read(parameter("feature_raster").toString());
        auto values = GdalIO::read(parameter("value_raster").toString());
        if (!features || !values) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = features->cols(), rows = features->rows();
        if (values->cols() != cols || values->rows() != rows) {
            setError("Feature and value rasters must have the same dimensions");
            return false;
        }

        int statIdx = parameter("statistic").toInt();
        const auto& dFeat = features->data(0);
        const auto& dVal = values->data(0);
        double noData = values->noDataValue();
        double featND = features->noDataValue();
        bool hasND = values->hasNoData();
        bool hasFeatND = features->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Accumulate statistics per zone
        struct ZoneStats {
            double sum = 0.0;
            double min = std::numeric_limits<double>::max();
            double max = std::numeric_limits<double>::lowest();
            int64_t count = 0;
        };
        QMap<int, ZoneStats> zones;

        for (int64_t i = 0; i < total; ++i) {
            if (hasFeatND && dFeat[i] == featND) continue;
            if (hasND && dVal[i] == noData) continue;

            int zone = static_cast<int>(dFeat[i]);
            auto& zs = zones[zone];
            zs.sum += dVal[i];
            if (dVal[i] < zs.min) zs.min = dVal[i];
            if (dVal[i] > zs.max) zs.max = dVal[i];
            zs.count++;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total * 0.8);
        }

        // Write output CSV
        QFile outFile(parameter("output_file").toString());
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file for writing");
            return false;
        }

        QStringList statNames = {"mean", "sum", "min", "max", "count"};
        QTextStream out(&outFile);
        out << "zone," << statNames[statIdx] << "\n";

        for (auto it = zones.constBegin(); it != zones.constEnd(); ++it) {
            const auto& zs = it.value();
            double result;
            switch (statIdx) {
                case 0: result = (zs.count > 0) ? zs.sum / zs.count : 0.0; break;
                case 1: result = zs.sum; break;
                case 2: result = zs.min; break;
                case 3: result = zs.max; break;
                case 4: result = static_cast<double>(zs.count); break;
                default: result = 0.0;
            }
            out << it.key() << "," << result << "\n";
        }

        outFile.close();
        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(ExtractModule)

} // namespace aplaceholder
