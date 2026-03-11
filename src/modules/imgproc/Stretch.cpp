#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>

namespace aplaceholder {

class StretchModule : public Module {
public:
    QString name() const override { return "STRETCH"; }
    QString description() const override {
        return "Contrast stretch. Enhances image contrast by remapping pixel values "
               "to the full display range using the selected method.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output stretched image"),
            ParameterDef::combo("method", "Stretch method",
                {"Linear", "Saturated", "Histogram equalization"}, 0,
                "Contrast stretch method to apply"),
            ParameterDef::real("saturation_percent", "Saturation percent", 2.0, 0.0, 50.0,
                "Percentage of pixels to saturate at each tail (used with Saturated method)"),
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
        int method = parameter("method").toInt();
        double satPct = parameter("saturation_percent").toDouble();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(raster->noDataValue());

        const auto& src = raster->data(0);
        auto& dst = output.data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        // Collect valid pixel values and compute statistics
        std::vector<double> validVals;
        validVals.reserve(total);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && src[i] == noData) continue;
            validVals.push_back(src[i]);
        }

        if (validVals.empty()) {
            setError("No valid pixels in input raster");
            return false;
        }

        if (method == 0) {
            // Linear: map [min,max] to [0,255]
            auto stats = raster->computeStats(0);
            double srcMin = stats.min;
            double srcMax = stats.max;
            double range = srcMax - srcMin;

            for (int64_t i = 0; i < total; ++i) {
                if (hasND && src[i] == noData) {
                    dst[i] = noData;
                    continue;
                }
                if (range == 0.0) {
                    dst[i] = 127.0;
                } else {
                    dst[i] = std::clamp((src[i] - srcMin) / range * 255.0, 0.0, 255.0);
                }
                if (i % 1000000 == 0)
                    reportProgress(static_cast<double>(i) / total);
            }
        } else if (method == 1) {
            // Saturated: map [mean - N*stddev, mean + N*stddev] to [0,255]
            auto stats = raster->computeStats(0);
            double srcMin = stats.mean - satPct * stats.stddev;
            double srcMax = stats.mean + satPct * stats.stddev;
            double range = srcMax - srcMin;

            for (int64_t i = 0; i < total; ++i) {
                if (hasND && src[i] == noData) {
                    dst[i] = noData;
                    continue;
                }
                if (range == 0.0) {
                    dst[i] = 127.0;
                } else {
                    double val = (src[i] - srcMin) / range * 255.0;
                    dst[i] = std::clamp(val, 0.0, 255.0);
                }
                if (i % 1000000 == 0)
                    reportProgress(static_cast<double>(i) / total);
            }
        } else if (method == 2) {
            // Histogram equalization: compute CDF, map values for uniform output
            reportProgress(0.0, "Computing histogram...");

            std::sort(validVals.begin(), validVals.end());
            int64_t nValid = static_cast<int64_t>(validVals.size());

            for (int64_t i = 0; i < total; ++i) {
                if (hasND && src[i] == noData) {
                    dst[i] = noData;
                    continue;
                }
                // Find rank of this pixel value via binary search (CDF position)
                auto lb = std::lower_bound(validVals.begin(), validVals.end(), src[i]);
                auto ub = std::upper_bound(validVals.begin(), validVals.end(), src[i]);
                // Use midpoint of the range of equal values for the CDF value
                int64_t rankLow = std::distance(validVals.begin(), lb);
                int64_t rankHigh = std::distance(validVals.begin(), ub);
                double cdf = static_cast<double>(rankLow + rankHigh) / (2.0 * nValid);
                dst[i] = std::clamp(cdf * 255.0, 0.0, 255.0);

                if (i % 1000000 == 0)
                    reportProgress(static_cast<double>(i) / total);
            }
        }

        // Build before/after histogram chart
        {
            const int numBins = 256;
            std::vector<double> beforeHist(numBins, 0.0);
            std::vector<double> afterHist(numBins, 0.0);

            // Determine bin edges for source values
            double srcMin = *std::min_element(validVals.begin(), validVals.end());
            double srcMax = *std::max_element(validVals.begin(), validVals.end());
            double srcRange = srcMax - srcMin;

            for (int64_t i = 0; i < total; ++i) {
                if (hasND && src[i] == noData) continue;

                // Before histogram
                int srcBin = (srcRange > 0.0)
                    ? static_cast<int>((src[i] - srcMin) / srcRange * (numBins - 1))
                    : numBins / 2;
                srcBin = std::clamp(srcBin, 0, numBins - 1);
                beforeHist[srcBin]++;

                // After histogram
                int dstBin = std::clamp(static_cast<int>(dst[i]), 0, numBins - 1);
                afterHist[dstBin]++;
            }

            ChartResult chart;
            chart.type = ChartResult::Histogram;
            chart.title = "Contrast Stretch: Before vs After";
            chart.xLabel = "Pixel Value";
            chart.yLabel = "Frequency";

            ChartSeries beforeSeries;
            beforeSeries.label = "Before";
            beforeSeries.color = ChartColor(70, 130, 180); // steel blue
            for (int i = 0; i < numBins; ++i) {
                beforeSeries.x.push_back(srcMin + i * srcRange / (numBins - 1));
                beforeSeries.y.push_back(beforeHist[i]);
            }
            chart.series.push_back(std::move(beforeSeries));

            ChartSeries afterSeries;
            afterSeries.label = "After";
            afterSeries.color = ChartColor(220, 80, 60); // crimson
            for (int i = 0; i < numBins; ++i) {
                afterSeries.x.push_back(static_cast<double>(i));
                afterSeries.y.push_back(afterHist[i]);
            }
            chart.series.push_back(std::move(afterSeries));

            setChartResult(std::move(chart));
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(StretchModule)

} // namespace aplaceholder
