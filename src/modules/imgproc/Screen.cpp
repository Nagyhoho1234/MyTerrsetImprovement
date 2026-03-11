#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

namespace aplaceholder {

class ScreenModule : public Module {
public:
    QString name() const override { return "SCREEN"; }
    QString description() const override {
        return "Hyperspectral Band Screening. Identifies and flags bands degraded by "
               "atmospheric absorption or noise using signal-to-noise ratio (SNR) and "
               "inter-band correlation analysis. Optionally outputs a cleaned image "
               "with bad bands removed.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input multi-band image or comma-separated band files",
                "Multi-band raster or comma-separated list of single-band files"),
            ParameterDef::output("output_report", "Output screening report",
                "Text file listing good and bad bands with statistics"),
            ParameterDef::output("output_clean", "Output cleaned image (optional)",
                "Output multi-band raster with bad bands removed"),
            ParameterDef::real("snr_threshold", "SNR threshold", 10.0, 0.0, 999999.0,
                "Bands with SNR below this threshold are flagged as bad"),
            ParameterDef::real("correlation_threshold", "Correlation threshold", 0.5, 0.0, 1.0,
                "Bands with neighbor correlation below this threshold are flagged"),
        };
    }

    bool execute() override {
        QString bandsInput = parameter("bands").toString();
        QString reportPath = parameter("output_report").toString();
        QString cleanPath = parameter("output_clean").toString();
        double snrThresh = parameter("snr_threshold").toDouble();
        double corrThresh = parameter("correlation_threshold").toDouble();

        // Load bands - either a single multi-band file or comma-separated files
        std::vector<std::unique_ptr<Raster>> loadedRasters;
        std::unique_ptr<Raster> multiBand;
        int nBands = 0;
        int cols = 0, rows = 0;

        QStringList bandFiles = bandsInput.split(",", Qt::SkipEmptyParts);
        for (auto& f : bandFiles) f = f.trimmed();

        if (bandFiles.size() == 1) {
            // Single multi-band file
            multiBand = GdalIO::read(bandFiles[0]);
            if (!multiBand) {
                setError("Failed to read input: " + bandFiles[0]);
                return false;
            }
            nBands = multiBand->bands();
            cols = multiBand->cols();
            rows = multiBand->rows();
        } else {
            // Multiple single-band files
            for (const auto& f : bandFiles) {
                auto r = GdalIO::read(f);
                if (!r) {
                    setError("Failed to read band file: " + f);
                    return false;
                }
                if (loadedRasters.empty()) {
                    cols = r->cols();
                    rows = r->rows();
                } else if (r->cols() != cols || r->rows() != rows) {
                    setError("Band dimensions mismatch for: " + f);
                    return false;
                }
                loadedRasters.push_back(std::move(r));
            }
            nBands = static_cast<int>(loadedRasters.size());
        }

        if (nBands < 2) {
            setError("At least 2 bands are required for screening");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;

        // Helper to get band data
        auto getBandData = [&](int b) -> const std::vector<double>& {
            if (multiBand)
                return multiBand->data(b);
            else
                return loadedRasters[b]->data(0);
        };

        bool hasND = multiBand ? multiBand->hasNoData() :
                     (!loadedRasters.empty() && loadedRasters[0]->hasNoData());
        double noData = multiBand ? multiBand->noDataValue() :
                        (!loadedRasters.empty() ? loadedRasters[0]->noDataValue() : -9999.0);

        reportProgress(0.1, "Computing band statistics...");

        // Compute per-band statistics: mean, stddev, SNR
        struct BandStats {
            double mean = 0.0;
            double stddev = 0.0;
            double snr = 0.0;
            double corrPrev = 1.0;  // correlation with previous band
            double corrNext = 1.0;  // correlation with next band
            double minCorr = 1.0;   // minimum neighbor correlation
            bool badSnr = false;
            bool badCorr = false;
            bool bad = false;
        };
        std::vector<BandStats> stats(nBands);

        for (int b = 0; b < nBands; ++b) {
            const auto& data = getBandData(b);
            double sum = 0.0, sumSq = 0.0;
            int64_t count = 0;

            for (int64_t i = 0; i < total; ++i) {
                if (hasND && data[i] == noData) continue;
                sum += data[i];
                sumSq += data[i] * data[i];
                count++;
            }

            if (count > 1) {
                stats[b].mean = sum / count;
                double variance = (sumSq - sum * sum / count) / (count - 1);
                stats[b].stddev = std::sqrt(std::max(0.0, variance));
                stats[b].snr = (stats[b].stddev > 1e-15) ?
                               (stats[b].mean / stats[b].stddev) : 0.0;
            }

            stats[b].badSnr = (std::abs(stats[b].snr) < snrThresh);

            reportProgress(0.1 + 0.3 * (b + 1.0) / nBands);
        }

        reportProgress(0.4, "Computing inter-band correlations...");

        // Compute correlations with neighboring bands
        auto computeCorrelation = [&](int b1, int b2) -> double {
            const auto& d1 = getBandData(b1);
            const auto& d2 = getBandData(b2);
            double m1 = stats[b1].mean, m2 = stats[b2].mean;

            double sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
            int64_t count = 0;
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && (d1[i] == noData || d2[i] == noData)) continue;
                double dx = d1[i] - m1;
                double dy = d2[i] - m2;
                sumXY += dx * dy;
                sumX2 += dx * dx;
                sumY2 += dy * dy;
                count++;
            }

            double denom = std::sqrt(sumX2 * sumY2);
            if (denom < 1e-15 || count < 2) return 0.0;
            return sumXY / denom;
        };

        for (int b = 0; b < nBands; ++b) {
            double minCorr = 1.0;
            if (b > 0) {
                stats[b].corrPrev = computeCorrelation(b, b - 1);
                minCorr = std::min(minCorr, stats[b].corrPrev);
            }
            if (b < nBands - 1) {
                stats[b].corrNext = computeCorrelation(b, b + 1);
                minCorr = std::min(minCorr, stats[b].corrNext);
            }
            stats[b].minCorr = minCorr;
            stats[b].badCorr = (minCorr < corrThresh);
            stats[b].bad = stats[b].badSnr || stats[b].badCorr;

            reportProgress(0.4 + 0.3 * (b + 1.0) / nBands);
        }

        reportProgress(0.7, "Writing report...");

        // Write screening report
        {
            std::ofstream rpt(reportPath.toStdString());
            if (!rpt.is_open()) {
                setError("Failed to open report file: " + reportPath);
                return false;
            }

            rpt << "Hyperspectral Band Screening Report\n";
            rpt << "====================================\n\n";
            rpt << "Number of bands: " << nBands << "\n";
            rpt << "Image dimensions: " << cols << " x " << rows << "\n";
            rpt << "SNR threshold: " << snrThresh << "\n";
            rpt << "Correlation threshold: " << corrThresh << "\n\n";

            int goodCount = 0, badCount = 0;
            for (int b = 0; b < nBands; ++b) {
                if (stats[b].bad) badCount++; else goodCount++;
            }
            rpt << "Good bands: " << goodCount << "\n";
            rpt << "Bad bands: " << badCount << "\n\n";

            rpt << "Band Details:\n";
            rpt << "-------------------------------------------------------------\n";
            rpt << "Band  Status  Mean        StdDev      SNR         MinCorr\n";
            rpt << "-------------------------------------------------------------\n";
            for (int b = 0; b < nBands; ++b) {
                char line[256];
                std::snprintf(line, sizeof(line),
                    "%4d  %-6s  %10.4f  %10.4f  %10.4f  %10.4f",
                    b + 1,
                    stats[b].bad ? "BAD" : "GOOD",
                    stats[b].mean, stats[b].stddev,
                    stats[b].snr, stats[b].minCorr);
                rpt << line;
                if (stats[b].badSnr) rpt << "  [low SNR]";
                if (stats[b].badCorr) rpt << "  [low corr]";
                rpt << "\n";
            }

            rpt << "\nGood band indices (1-based): ";
            bool first = true;
            for (int b = 0; b < nBands; ++b) {
                if (!stats[b].bad) {
                    if (!first) rpt << ", ";
                    rpt << (b + 1);
                    first = false;
                }
            }
            rpt << "\n";

            rpt << "Bad band indices (1-based): ";
            first = true;
            for (int b = 0; b < nBands; ++b) {
                if (stats[b].bad) {
                    if (!first) rpt << ", ";
                    rpt << (b + 1);
                    first = false;
                }
            }
            rpt << "\n";
        }

        // Optionally write cleaned multi-band raster
        if (!cleanPath.isEmpty()) {
            reportProgress(0.8, "Writing cleaned image...");

            std::vector<int> goodBands;
            for (int b = 0; b < nBands; ++b) {
                if (!stats[b].bad) goodBands.push_back(b);
            }

            if (goodBands.empty()) {
                setError("All bands flagged as bad - cannot produce cleaned output");
                return false;
            }

            int nGood = static_cast<int>(goodBands.size());

            // Get geo info from source
            const Raster* srcRef = multiBand ? multiBand.get() :
                                   loadedRasters[0].get();

            Raster cleaned(cols, rows, nGood, DataType::Float64);
            cleaned.setGeoTransform(srcRef->geoTransform());
            cleaned.setProjection(srcRef->projection());
            cleaned.setNoDataValue(noData);

            for (int g = 0; g < nGood; ++g) {
                const auto& srcData = getBandData(goodBands[g]);
                auto& dstData = cleaned.data(g);
                for (int64_t i = 0; i < total; ++i)
                    dstData[i] = srcData[i];

                reportProgress(0.8 + 0.15 * (g + 1.0) / nGood);
            }

            if (!GdalIO::write(cleaned, cleanPath)) {
                setError("Failed to write cleaned output");
                return false;
            }
        }

        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(ScreenModule)

} // namespace aplaceholder
