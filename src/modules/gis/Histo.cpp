#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <cmath>
#include <vector>

namespace aplaceholder {

class HistoModule : public Module {
public:
    QString name() const override { return "HISTO"; }
    QString description() const override {
        return "Computes a frequency distribution histogram of raster image values. "
               "Outputs a text file with value, count, and percentage columns.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::output("output_file", "Output text file"),
            ParameterDef::integer("num_bins", "Number of bins", 256, 2, 10000,
                "Number of histogram bins"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int numBins = parameter("num_bins").toInt();
        if (numBins < 2) numBins = 256;

        const auto& data = raster->data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();
        int64_t total = raster->cellCount();

        // First pass: find min/max of valid values
        double minVal = std::numeric_limits<double>::max();
        double maxVal = std::numeric_limits<double>::lowest();
        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            if (std::isnan(data[i])) continue;
            if (data[i] < minVal) minVal = data[i];
            if (data[i] > maxVal) maxVal = data[i];
            validCount++;
        }

        if (validCount == 0) {
            setError("No valid data values found in raster");
            return false;
        }

        reportProgress(0.3, "Computing histogram...");

        // Compute bin width
        double range = maxVal - minVal;
        double binWidth = (range > 0.0) ? range / numBins : 1.0;

        // If all values are the same, use a single bin
        if (range == 0.0) {
            numBins = 1;
        }

        std::vector<int64_t> counts(numBins, 0);

        // Second pass: bin the values
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            if (std::isnan(data[i])) continue;

            int bin;
            if (range == 0.0) {
                bin = 0;
            } else {
                bin = static_cast<int>((data[i] - minVal) / binWidth);
                if (bin >= numBins) bin = numBins - 1;  // Handle max value edge case
            }
            counts[bin]++;

            if (i % 1000000 == 0)
                reportProgress(0.3 + 0.5 * static_cast<double>(i) / total);
        }

        // Write output
        QFile outFile(parameter("output_file").toString());
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file for writing");
            return false;
        }

        QTextStream out(&outFile);
        out << "value,count,percentage\n";

        for (int b = 0; b < numBins; ++b) {
            double binCenter = minVal + (b + 0.5) * binWidth;
            double pct = (validCount > 0)
                ? (static_cast<double>(counts[b]) / validCount) * 100.0
                : 0.0;
            out << binCenter << "," << counts[b] << "," << pct << "\n";
        }

        outFile.close();

        // Build chart result
        ChartResult chart;
        chart.type = ChartResult::Histogram;
        chart.title = "Frequency Histogram";
        chart.xLabel = "Value";
        chart.yLabel = "Count";
        ChartSeries s;
        s.label = "Frequency";
        for (int b = 0; b < numBins; ++b) {
            double binCenter = minVal + (b + 0.5) * binWidth;
            s.x.push_back(binCenter);
            s.y.push_back(static_cast<double>(counts[b]));
        }
        chart.series.push_back(s);
        setChartResult(chart);

        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(HistoModule)

} // namespace aplaceholder
