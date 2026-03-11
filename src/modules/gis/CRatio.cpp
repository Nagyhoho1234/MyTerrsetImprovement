#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class CRatioModule : public Module {
public:
    QString name() const override { return "CRATIO"; }
    QString description() const override {
        return "Computes a Gini-like concentration ratio for spatial data. "
               "Measures the degree of inequality in the distribution of values "
               "across the input raster. Outputs a single-value result raster "
               "containing the concentration ratio (0 = perfect equality, 1 = maximum concentration).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::output("output", "Output single-value result raster"),
        };
    }

    bool execute() override {
        auto rIn = GdalIO::read(parameter("input").toString());
        if (!rIn) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = rIn->cols(), rows = rIn->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& data = rIn->data(0);
        double noData = rIn->noDataValue();
        bool hasND = rIn->hasNoData();

        reportProgress(0.0, "Collecting valid pixel values...");

        // Collect all valid (non-NoData) values
        std::vector<double> values;
        values.reserve(total);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            if (data[i] < 0.0) {
                setError("Concentration ratio requires non-negative values");
                return false;
            }
            values.push_back(data[i]);
        }

        int64_t N = static_cast<int64_t>(values.size());
        if (N < 2) {
            setError("Insufficient valid pixels for concentration ratio (need at least 2)");
            return false;
        }

        reportProgress(0.2, "Sorting values...");

        std::sort(values.begin(), values.end());

        reportProgress(0.5, "Computing concentration ratio...");

        // Gini coefficient via the sorted-values formula:
        // G = (2 * sum_i(i * y_i)) / (N * sum(y_i)) - (N + 1) / N
        double Nd = static_cast<double>(N);
        double sumValues = 0.0;
        double weightedSum = 0.0;

        for (int64_t i = 0; i < N; ++i) {
            sumValues += values[i];
            weightedSum += static_cast<double>(i + 1) * values[i];
        }

        double concentrationRatio = 0.0;
        if (sumValues > 0.0) {
            concentrationRatio = (2.0 * weightedSum) / (Nd * sumValues) - (Nd + 1.0) / Nd;
        }

        reportProgress(0.8, "Writing output raster...");

        // Output a single-value 1x1 raster with the concentration ratio
        double outNoData = -9999.0;

        Raster result(1, 1, 1, DataType::Float64);
        result.setGeoTransform(rIn->geoTransform());
        result.setProjection(rIn->projection());
        result.setNoDataValue(outNoData);
        auto& resData = result.data(0);
        resData[0] = concentrationRatio;

        if (!GdalIO::write(result, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0,
            QString("Concentration ratio = %1 (N = %2)")
                .arg(concentrationRatio, 0, 'f', 6)
                .arg(N));

        return true;
    }
};

REGISTER_MODULE(CRatioModule)

} // namespace aplaceholder
