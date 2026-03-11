#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>

namespace aplaceholder {

class WindowModule : public Module {
public:
    QString name() const override { return "WINDOW"; }
    QString description() const override {
        return "Moving window (focal) statistics. Computes a statistic within "
               "a square window centered on each pixel. Supports mean, median, "
               "min, max, range, standard deviation, sum, and variety.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::combo("statistic", "Focal statistic",
                {"mean", "median", "min", "max", "range", "std", "sum", "variety"}, 0,
                "Statistic to compute within the moving window"),
            ParameterDef::integer("window_size", "Window size (odd number)",
                3, 3, 101,
                "Side length of the square window in pixels (must be odd)"),
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
        const auto& data = raster->data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        int statIdx = parameter("statistic").toInt();
        int winSize = parameter("window_size").toInt();
        if (winSize < 3) winSize = 3;
        if (winSize % 2 == 0) winSize++; // ensure odd
        int halfWin = winSize / 2;

        reportProgress(0.0, "Computing moving window statistics...");

        double outNoData = -9999.0;
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(outNoData);
        auto& out = output.data(0);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int64_t idx = static_cast<int64_t>(r) * cols + c;

                // Skip if center is NoData
                if (hasND && data[idx] == noData) {
                    out[idx] = outNoData;
                    continue;
                }

                // Collect valid values in window
                std::vector<double> vals;
                vals.reserve(winSize * winSize);
                for (int dr = -halfWin; dr <= halfWin; ++dr) {
                    for (int dc = -halfWin; dc <= halfWin; ++dc) {
                        int nr = r + dr, nc = c + dc;
                        if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                        int64_t ni = static_cast<int64_t>(nr) * cols + nc;
                        if (hasND && data[ni] == noData) continue;
                        vals.push_back(data[ni]);
                    }
                }

                if (vals.empty()) {
                    out[idx] = outNoData;
                    continue;
                }

                double result = 0.0;
                switch (statIdx) {
                    case 0: { // mean
                        double s = 0.0;
                        for (double v : vals) s += v;
                        result = s / vals.size();
                        break;
                    }
                    case 1: { // median
                        std::sort(vals.begin(), vals.end());
                        size_t n = vals.size();
                        result = (n % 2 == 0) ?
                            (vals[n / 2 - 1] + vals[n / 2]) / 2.0 : vals[n / 2];
                        break;
                    }
                    case 2: { // min
                        result = *std::min_element(vals.begin(), vals.end());
                        break;
                    }
                    case 3: { // max
                        result = *std::max_element(vals.begin(), vals.end());
                        break;
                    }
                    case 4: { // range
                        auto mm = std::minmax_element(vals.begin(), vals.end());
                        result = *mm.second - *mm.first;
                        break;
                    }
                    case 5: { // std
                        double s = 0.0, ss = 0.0;
                        for (double v : vals) { s += v; ss += v * v; }
                        double mean = s / vals.size();
                        double var = ss / vals.size() - mean * mean;
                        result = std::sqrt(std::max(var, 0.0));
                        break;
                    }
                    case 6: { // sum
                        double s = 0.0;
                        for (double v : vals) s += v;
                        result = s;
                        break;
                    }
                    case 7: { // variety (number of unique values)
                        std::set<double> unique(vals.begin(), vals.end());
                        result = static_cast<double>(unique.size());
                        break;
                    }
                }

                out[idx] = result;
            }

            if (r % 100 == 0)
                reportProgress(0.9 * static_cast<double>(r) / rows);
        }

        reportProgress(0.92, "Writing output raster...");

        if (!GdalIO::write(output, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        static const char* statNames[] = {
            "mean", "median", "min", "max", "range", "std", "sum", "variety"
        };
        reportProgress(1.0,
            QString("Window %1 (%2x%2) complete")
                .arg(statNames[statIdx])
                .arg(winSize));

        return true;
    }
};

REGISTER_MODULE(WindowModule)

} // namespace aplaceholder
