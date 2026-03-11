#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <set>
#include <cmath>

namespace aplaceholder {

class SplitModule : public Module {
public:
    QString name() const override { return "SPLIT"; }
    QString description() const override {
        return "Split raster by unique values. For each unique integer value, "
               "creates a binary raster.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef({
                "output_prefix", "Output prefix",
                ParameterDef::String, {}, {}, 0, 0,
                "Creates prefix_value<N>.tif for each unique value", true
            }),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& src = raster->data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        // Collect unique integer values
        reportProgress(0.0, "Finding unique values...");
        std::set<int> uniqueValues;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && src[i] == noData)
                continue;
            int val = static_cast<int>(std::round(src[i]));
            uniqueValues.insert(val);
        }

        if (uniqueValues.empty()) {
            setError("No non-NoData values found in input raster");
            return false;
        }

        QString prefix = parameter("output_prefix").toString();
        int count = 0;
        int numValues = static_cast<int>(uniqueValues.size());

        for (int val : uniqueValues) {
            Raster binary(cols, rows, 1, DataType::Byte);
            binary.setGeoTransform(raster->geoTransform());
            binary.setProjection(raster->projection());

            auto& dst = binary.data(0);
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && src[i] == noData) {
                    dst[i] = 0.0;
                } else {
                    dst[i] = (static_cast<int>(std::round(src[i])) == val) ? 1.0 : 0.0;
                }
            }

            QString outPath = QString("%1_value%2.tif").arg(prefix).arg(val);
            if (!GdalIO::write(binary, outPath)) {
                setError("Failed to write output for value " + QString::number(val));
                return false;
            }

            ++count;
            reportProgress(static_cast<double>(count) / numValues);
        }

        reportProgress(1.0, "Done");
        return true;
    }
};

REGISTER_MODULE(SplitModule)

} // namespace aplaceholder
