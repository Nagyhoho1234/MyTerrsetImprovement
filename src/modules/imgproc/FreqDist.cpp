#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <QFile>
#include <QTextStream>

namespace aplaceholder {

class FreqDistModule : public Module {
public:
    QString name() const override { return "FREQDIST"; }
    QString description() const override {
        return "Frequency spectrum distance. Compares two magnitude spectra (from FFT) "
               "by computing Euclidean distance and correlation between them.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("magnitude1", "First magnitude spectrum",
                "Magnitude image from FOURIER forward transform"),
            ParameterDef::file("magnitude2", "Second magnitude spectrum",
                "Magnitude image from FOURIER forward transform"),
            ParameterDef::output("output_report", "Output report file (.txt)"),
        };
    }

    bool execute() override {
        auto mag1 = GdalIO::read(parameter("magnitude1").toString());
        if (!mag1) {
            setError("Failed to read first magnitude spectrum");
            return false;
        }

        auto mag2 = GdalIO::read(parameter("magnitude2").toString());
        if (!mag2) {
            setError("Failed to read second magnitude spectrum");
            return false;
        }

        int cols1 = mag1->cols(), rows1 = mag1->rows();
        int cols2 = mag2->cols(), rows2 = mag2->rows();

        if (cols1 != cols2 || rows1 != rows2) {
            setError("Magnitude spectra must have the same dimensions");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols1) * rows1;
        const auto& data1 = mag1->data(0);
        const auto& data2 = mag2->data(0);
        bool hasND1 = mag1->hasNoData();
        bool hasND2 = mag2->hasNoData();
        double nd1 = mag1->noDataValue();
        double nd2 = mag2->noDataValue();

        reportProgress(0.1, "Computing spectral distance metrics...");

        // Compute Euclidean distance, correlation, and basic stats
        double sumDiffSq = 0.0;
        double sum1 = 0.0, sum2 = 0.0;
        double sumSq1 = 0.0, sumSq2 = 0.0;
        double sumProd = 0.0;
        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            double v1 = data1[i];
            double v2 = data2[i];

            if (hasND1 && v1 == nd1) continue;
            if (hasND2 && v2 == nd2) continue;

            double diff = v1 - v2;
            sumDiffSq += diff * diff;
            sum1 += v1;
            sum2 += v2;
            sumSq1 += v1 * v1;
            sumSq2 += v2 * v2;
            sumProd += v1 * v2;
            validCount++;

            if (i % (total / 20 + 1) == 0)
                reportProgress(0.1 + 0.7 * static_cast<double>(i) / total);
        }

        if (validCount == 0) {
            setError("No valid pixels found for comparison");
            return false;
        }

        double euclideanDist = std::sqrt(sumDiffSq);
        double rmsDist = std::sqrt(sumDiffSq / validCount);
        double mean1 = sum1 / validCount;
        double mean2 = sum2 / validCount;

        // Pearson correlation coefficient
        double numerator = sumProd - (sum1 * sum2 / validCount);
        double denom1 = sumSq1 - (sum1 * sum1 / validCount);
        double denom2 = sumSq2 - (sum2 * sum2 / validCount);
        double denomProduct = denom1 * denom2;
        double correlation = 0.0;
        if (denomProduct > 0.0) {
            correlation = numerator / std::sqrt(denomProduct);
        }

        // Build histogram chart of per-pixel magnitude differences
        {
            const int numBins = 100;
            // First pass: find min/max difference
            double diffMin = std::numeric_limits<double>::max();
            double diffMax = std::numeric_limits<double>::lowest();
            for (int64_t i = 0; i < total; ++i) {
                double v1 = data1[i], v2 = data2[i];
                if (hasND1 && v1 == nd1) continue;
                if (hasND2 && v2 == nd2) continue;
                double d = v1 - v2;
                if (d < diffMin) diffMin = d;
                if (d > diffMax) diffMax = d;
            }

            double binWidth = (diffMax - diffMin) / numBins;
            if (binWidth <= 0.0) binWidth = 1.0;

            std::vector<double> hist(numBins, 0.0);
            for (int64_t i = 0; i < total; ++i) {
                double v1 = data1[i], v2 = data2[i];
                if (hasND1 && v1 == nd1) continue;
                if (hasND2 && v2 == nd2) continue;
                int bin = static_cast<int>((v1 - v2 - diffMin) / binWidth);
                bin = std::clamp(bin, 0, numBins - 1);
                hist[bin]++;
            }

            ChartResult chart;
            chart.type = ChartResult::Histogram;
            chart.title = "Frequency Spectrum Difference Distribution";
            chart.xLabel = "Magnitude Difference";
            chart.yLabel = "Frequency";

            ChartSeries series;
            series.label = "Difference";
            series.color = ChartColor(100, 100, 200);
            for (int i = 0; i < numBins; ++i) {
                series.x.push_back(diffMin + (i + 0.5) * binWidth);
                series.y.push_back(hist[i]);
            }
            chart.series.push_back(std::move(series));
            setChartResult(std::move(chart));
        }

        reportProgress(0.9, "Writing report...");

        // Write report
        QString reportPath = parameter("output_report").toString();
        QFile reportFile(reportPath);
        if (!reportFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        QTextStream out(&reportFile);
        out << "Frequency Spectrum Distance Report\n";
        out << "===================================\n\n";
        out << "Spectrum 1: " << parameter("magnitude1").toString() << "\n";
        out << "Spectrum 2: " << parameter("magnitude2").toString() << "\n\n";
        out << "Dimensions: " << cols1 << " x " << rows1 << "\n";
        out << "Valid pixels compared: " << validCount << "\n\n";
        out << "--- Distance Metrics ---\n";
        out << "Euclidean distance:    " << euclideanDist << "\n";
        out << "RMS distance:          " << rmsDist << "\n";
        out << "Pearson correlation:   " << correlation << "\n\n";
        out << "--- Summary Statistics ---\n";
        out << "Spectrum 1 mean: " << mean1 << "\n";
        out << "Spectrum 2 mean: " << mean2 << "\n";

        reportFile.close();

        reportProgress(1.0, "Report written.");
        return true;
    }
};

REGISTER_MODULE(FreqDistModule)

} // namespace aplaceholder
