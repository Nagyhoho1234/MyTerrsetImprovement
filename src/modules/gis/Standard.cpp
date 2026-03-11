#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class StandardModule : public Module {
public:
    QString name() const override { return "STANDARD"; }
    QString description() const override {
        return "Standardize raster to z-scores: (x - mean) / stddev.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output image"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input").toString());
        if (!r1) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        const auto& d1 = r1->data(0);
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Pass 1: compute mean
        double sum = 0.0;
        int64_t count = 0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && d1[i] == noData) continue;
            sum += d1[i];
            ++count;
        }

        if (count == 0) {
            setError("No valid pixels found");
            return false;
        }

        double mean = sum / count;

        // Pass 1b: compute stddev
        double sumSq = 0.0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && d1[i] == noData) continue;
            double diff = d1[i] - mean;
            sumSq += diff * diff;
        }
        double stddev = std::sqrt(sumSq / count);
        if (stddev == 0.0) stddev = 1.0; // avoid division by zero

        // Pass 2: compute z-scores
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(noData);

        auto& out = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && d1[i] == noData) {
                out[i] = noData;
                continue;
            }
            out[i] = (d1[i] - mean) / stddev;

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(StandardModule)

} // namespace aplaceholder
