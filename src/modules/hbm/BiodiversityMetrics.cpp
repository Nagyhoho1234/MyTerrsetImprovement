#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

#include <cmath>
#include <map>
#include <set>
#include <algorithm>

namespace aplaceholder {

class BiodiversityMetricsModule : public Module {
public:
    QString name() const override { return "BIODIVERSITY_METRICS"; }
    QString description() const override {
        return "Computes landscape biodiversity metrics from a land cover classification "
               "raster using a moving window. Supports Shannon Diversity Index, Simpson's "
               "Diversity Index, and Species Richness (unique class count).";
    }
    QString category() const override { return "Habitat & Biodiversity"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input land cover raster",
                "Categorical integer land cover classification raster"),
            ParameterDef::output("output", "Output diversity raster",
                "Output raster containing computed diversity metric values"),
            ParameterDef::combo("metric", "Diversity metric",
                {"Shannon", "Simpson", "Richness"}, 0,
                "Shannon: H = -sum(pi*ln(pi)); "
                "Simpson: D = 1 - sum(pi^2); "
                "Richness: count of unique classes"),
            ParameterDef::integer("window_size", "Window size (cells)",
                3, 3, 101,
                "Side length of square moving window (must be odd). "
                "Corresponds to 3x3, 5x5, 7x7, etc. neighborhood"),
        };
    }

    bool execute() override {
        // --- Read parameters ---------------------------------------------------
        const QString inputPath  = parameter("input").toString();
        const QString outputPath = parameter("output").toString();
        const int metricIndex    = parameter("metric").toInt();
        int windowSize           = parameter("window_size").toInt();

        if (inputPath.isEmpty() || outputPath.isEmpty()) {
            setError("Input raster and output path are required.");
            return false;
        }

        // Force window size to be odd
        if (windowSize % 2 == 0) {
            windowSize += 1;
        }
        if (windowSize < 3) windowSize = 3;

        const int halfWin = windowSize / 2;

        // --- Load input raster -------------------------------------------------
        auto inputRaster = GdalIO::read(inputPath);
        if (!inputRaster) {
            setError("Failed to read input raster: " + inputPath);
            return false;
        }

        const int cols = inputRaster->cols();
        const int rows = inputRaster->rows();
        const bool hasNoData = inputRaster->hasNoData();
        const double noData  = inputRaster->noDataValue();

        // --- Create output raster (same dimensions and georef) -----------------
        auto output = std::make_unique<Raster>(cols, rows, 1, DataType::Float32);
        output->setGeoTransform(inputRaster->geoTransform());
        output->setProjection(inputRaster->projection());
        output->setNoDataValue(noData);
        output->allocate();

        const auto& inData  = inputRaster->data(0);
        auto&       outData = output->data(0);

        const int64_t totalCells = static_cast<int64_t>(rows) * cols;

        // --- Moving window analysis -------------------------------------------
        for (int r = 0; r < rows; ++r) {
            // Report progress every row
            reportProgress(static_cast<double>(r) / rows,
                           QString("Processing row %1 / %2").arg(r + 1).arg(rows));

            for (int c = 0; c < cols; ++c) {
                const int64_t idx = static_cast<int64_t>(r) * cols + c;

                // If the centre pixel is nodata, output nodata
                if (hasNoData && inData[idx] == noData) {
                    outData[idx] = noData;
                    continue;
                }

                // Gather class counts within the window
                std::map<int, int> classCounts;
                int validPixels = 0;

                const int rStart = std::max(0, r - halfWin);
                const int rEnd   = std::min(rows - 1, r + halfWin);
                const int cStart = std::max(0, c - halfWin);
                const int cEnd   = std::min(cols - 1, c + halfWin);

                for (int wr = rStart; wr <= rEnd; ++wr) {
                    for (int wc = cStart; wc <= cEnd; ++wc) {
                        const int64_t wIdx = static_cast<int64_t>(wr) * cols + wc;
                        const double val = inData[wIdx];

                        if (hasNoData && val == noData) continue;

                        int classVal = static_cast<int>(std::round(val));
                        classCounts[classVal]++;
                        validPixels++;
                    }
                }

                if (validPixels == 0) {
                    outData[idx] = noData;
                    continue;
                }

                // Compute the selected metric
                double result = 0.0;

                if (metricIndex == 0) {
                    // Shannon Diversity Index: H = -sum(pi * ln(pi))
                    for (const auto& kv : classCounts) {
                        double pi = static_cast<double>(kv.second) / validPixels;
                        if (pi > 0.0) {
                            result -= pi * std::log(pi);
                        }
                    }
                } else if (metricIndex == 1) {
                    // Simpson's Diversity Index: D = 1 - sum(pi^2)
                    double sumPi2 = 0.0;
                    for (const auto& kv : classCounts) {
                        double pi = static_cast<double>(kv.second) / validPixels;
                        sumPi2 += pi * pi;
                    }
                    result = 1.0 - sumPi2;
                } else {
                    // Species Richness: count of unique classes
                    result = static_cast<double>(classCounts.size());
                }

                outData[idx] = result;
            }
        }

        // --- Write output ------------------------------------------------------
        if (!GdalIO::write(*output, outputPath)) {
            setError("Failed to write output raster: " + outputPath);
            return false;
        }

        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(BiodiversityMetricsModule)

} // namespace aplaceholder
