#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <map>
#include <set>
#include <vector>

namespace aplaceholder {

class CrosstabModule : public Module {
public:
    QString name() const override { return "CROSSTAB"; }
    QString description() const override {
        return "Cross-tabulation of two categorical raster images. "
               "Produces a matrix of joint frequency distributions "
               "with Chi-square and Cramer's V statistics.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First input image"),
            ParameterDef::file("input2", "Second input image"),
            ParameterDef::output("output", "Output cross-tabulation image"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input1").toString());
        auto r2 = GdalIO::read(parameter("input2").toString());
        if (!r1 || !r2) {
            setError("Failed to read input rasters");
            return false;
        }

        if (r1->cols() != r2->cols() || r1->rows() != r2->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& d1 = r1->data(0);
        const auto& d2 = r2->data(0);
        double noData1 = r1->noDataValue();
        double noData2 = r2->noDataValue();
        bool hasND1 = r1->hasNoData();
        bool hasND2 = r2->hasNoData();

        // Collect unique classes and build frequency table
        std::set<int> classes1, classes2;
        // Key: (class_from_image1, class_from_image2) -> count
        std::map<std::pair<int,int>, int64_t> freqTable;

        reportProgress(0.0, "Building frequency table...");
        int64_t validCount = 0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND1 && d1[i] == noData1) continue;
            if (hasND2 && d2[i] == noData2) continue;

            int c1 = static_cast<int>(d1[i]);
            int c2 = static_cast<int>(d2[i]);
            classes1.insert(c1);
            classes2.insert(c2);
            freqTable[{c1, c2}]++;
            validCount++;

            if (i % 1000000 == 0)
                reportProgress(0.3 * static_cast<double>(i) / total);
        }

        if (validCount == 0) {
            setError("No valid pixels found");
            return false;
        }

        // Assign a unique combo ID to each (c1, c2) pair
        std::map<std::pair<int,int>, int> comboId;
        int nextId = 1;
        for (int c1 : classes1) {
            for (int c2 : classes2) {
                if (freqTable.count({c1, c2}) && freqTable[{c1, c2}] > 0) {
                    comboId[{c1, c2}] = nextId++;
                }
            }
        }

        // Create output cross-tab raster: each pixel gets its combo ID
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        double outNoData = -9999;
        output.setNoDataValue(outNoData);
        auto& out = output.data(0);

        reportProgress(0.3, "Writing combo IDs...");
        for (int64_t i = 0; i < total; ++i) {
            if ((hasND1 && d1[i] == noData1) || (hasND2 && d2[i] == noData2)) {
                out[i] = outNoData;
                continue;
            }
            int c1 = static_cast<int>(d1[i]);
            int c2 = static_cast<int>(d2[i]);
            auto it = comboId.find({c1, c2});
            out[i] = (it != comboId.end()) ? it->second : outNoData;

            if (i % 1000000 == 0)
                reportProgress(0.3 + 0.4 * static_cast<double>(i) / total);
        }

        // Compute Chi-square and Cramer's V
        reportProgress(0.7, "Computing statistics...");

        // Row and column marginals
        std::map<int, int64_t> rowMargin, colMargin;
        for (const auto& kv : freqTable) {
            rowMargin[kv.first.first] += kv.second;
            colMargin[kv.first.second] += kv.second;
        }

        double chiSquare = 0.0;
        for (const auto& kv : freqTable) {
            int c1 = kv.first.first;
            int c2 = kv.first.second;
            int64_t observed = kv.second;
            double expected = static_cast<double>(rowMargin[c1]) *
                              static_cast<double>(colMargin[c2]) /
                              static_cast<double>(validCount);
            if (expected > 0) {
                double diff = static_cast<double>(observed) - expected;
                chiSquare += (diff * diff) / expected;
            }
        }

        int r = static_cast<int>(classes1.size());
        int k = static_cast<int>(classes2.size());
        int minDim = std::min(r, k);
        double cramersV = 0.0;
        if (minDim > 1 && validCount > 0) {
            cramersV = std::sqrt(chiSquare / (static_cast<double>(validCount) * (minDim - 1)));
        }

        // Report the frequency table and statistics via progress messages
        QString statsMsg = QString("Cross-tabulation complete. "
                                   "Classes in image 1: %1, Classes in image 2: %2, "
                                   "Chi-square: %3, Cramer's V: %4")
                           .arg(r).arg(k)
                           .arg(chiSquare, 0, 'f', 4)
                           .arg(cramersV, 0, 'f', 4);

        // Print frequency table header
        QString tableMsg = "Frequency table:\n";
        tableMsg += QString("Image1\\Image2");
        for (int c2 : classes2)
            tableMsg += QString("\t%1").arg(c2);
        tableMsg += "\n";

        for (int c1 : classes1) {
            tableMsg += QString("%1").arg(c1);
            for (int c2 : classes2) {
                auto it = freqTable.find({c1, c2});
                int64_t count = (it != freqTable.end()) ? it->second : 0;
                tableMsg += QString("\t%1").arg(count);
            }
            tableMsg += "\n";
        }

        // Build bar chart of category combination counts
        ChartResult chart;
        chart.type = ChartResult::Bar;
        chart.title = QString("Cross-tabulation Counts (Chi² = %1, V = %2)")
                      .arg(chiSquare, 0, 'f', 2).arg(cramersV, 0, 'f', 4);
        chart.xLabel = "Category Combination";
        chart.yLabel = "Count";
        ChartSeries s;
        s.label = "Frequency";
        int barIdx = 0;
        for (int c1 : classes1) {
            for (int c2 : classes2) {
                auto it = freqTable.find({c1, c2});
                if (it != freqTable.end() && it->second > 0) {
                    s.x.push_back(static_cast<double>(barIdx));
                    s.y.push_back(static_cast<double>(it->second));
                    chart.categoryLabels.push_back(
                        QString("%1x%2").arg(c1).arg(c2));
                    barIdx++;
                }
            }
        }
        chart.series.push_back(s);
        setChartResult(chart);

        reportProgress(1.0, statsMsg + "\n" + tableMsg);
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(CrosstabModule)

} // namespace aplaceholder
