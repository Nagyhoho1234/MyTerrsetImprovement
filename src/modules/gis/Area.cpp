#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <map>

namespace aplaceholder {

class AreaModule : public Module {
public:
    QString name() const override { return "AREA"; }
    QString description() const override {
        return "Area calculation for categorical raster images. "
               "Computes the area occupied by each category and writes "
               "an output raster where each pixel holds its class total area.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input categorical raster image"),
            ParameterDef::output("output", "Output area results"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();

        // Compute pixel area from geo-transform
        double pixelWidth = std::abs(input->geoTransform().pixelWidth);
        double pixelHeight = std::abs(input->geoTransform().pixelHeight);
        double pixelArea = pixelWidth * pixelHeight;

        // Count pixels per class
        reportProgress(0.0, "Counting pixels per class...");
        std::map<int, int64_t> classCounts;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            int classId = static_cast<int>(data[i]);
            classCounts[classId]++;

            if (i % 1000000 == 0)
                reportProgress(0.3 * static_cast<double>(i) / total);
        }

        // Compute total area per class
        std::map<int, double> classArea;
        for (const auto& kv : classCounts) {
            classArea[kv.first] = static_cast<double>(kv.second) * pixelArea;
        }

        // Create output raster: each pixel gets its class's total area
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        reportProgress(0.3, "Writing area values...");
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) {
                out[i] = noData;
                continue;
            }
            int classId = static_cast<int>(data[i]);
            out[i] = classArea[classId];

            if (i % 1000000 == 0)
                reportProgress(0.3 + 0.6 * static_cast<double>(i) / total);
        }

        // Report summary via progress message
        QString summary = "Area summary:\n";
        for (const auto& kv : classArea) {
            summary += QString("  Class %1: %2 pixels, area = %3\n")
                       .arg(kv.first)
                       .arg(classCounts[kv.first])
                       .arg(kv.second, 0, 'f', 4);
        }

        // Build pie chart of category areas
        ChartResult chart;
        chart.type = ChartResult::Pie;
        chart.title = "Area by Category";
        ChartSeries s;
        s.label = "Area";
        for (const auto& kv : classArea) {
            s.x.push_back(static_cast<double>(kv.first));
            s.y.push_back(kv.second);
            chart.categoryLabels.push_back(
                QString("Class %1").arg(kv.first));
        }
        chart.series.push_back(s);
        setChartResult(chart);

        reportProgress(1.0, summary + "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(AreaModule)

} // namespace aplaceholder
