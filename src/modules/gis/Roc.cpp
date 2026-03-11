#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <algorithm>

namespace aplaceholder {

class RocModule : public Module {
public:
    QString name() const override { return "ROC"; }
    QString description() const override {
        return "ROC (Receiver Operating Characteristic) curve analysis. "
               "Given a continuous prediction raster and a binary reference, "
               "computes TPR/FPR at multiple thresholds and calculates AUC.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("prediction", "Continuous prediction raster"),
            ParameterDef::file("reference", "Binary reference raster (0/1)"),
            ParameterDef::output("output_report", "Output ROC report text file"),
            ParameterDef::integer("num_thresholds", "Number of thresholds",
                100, 10, 10000,
                "Number of threshold values to evaluate along the ROC curve"),
        };
    }

    bool execute() override {
        auto rPred = GdalIO::read(parameter("prediction").toString());
        auto rRef = GdalIO::read(parameter("reference").toString());
        if (!rPred || !rRef) {
            setError("Failed to read input rasters");
            return false;
        }

        if (rPred->cols() != rRef->cols() || rPred->rows() != rRef->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = rPred->cols(), rows = rPred->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& dPred = rPred->data(0);
        const auto& dRef = rRef->data(0);
        double noDataP = rPred->noDataValue();
        double noDataR = rRef->noDataValue();
        bool hasNDP = rPred->hasNoData();
        bool hasNDR = rRef->hasNoData();
        int numThresholds = parameter("num_thresholds").toInt();
        if (numThresholds <= 0) numThresholds = 100;

        reportProgress(0.0, "Collecting valid predictions...");

        // Find min/max of prediction values and count positives/negatives
        double predMin = 1e300, predMax = -1e300;
        int64_t totalPositive = 0, totalNegative = 0;

        for (int64_t i = 0; i < total; ++i) {
            if (hasNDP && dPred[i] == noDataP) continue;
            if (hasNDR && dRef[i] == noDataR) continue;
            double p = dPred[i];
            int ref = static_cast<int>(std::round(dRef[i]));
            if (ref != 0 && ref != 1) continue; // skip non-binary reference

            if (p < predMin) predMin = p;
            if (p > predMax) predMax = p;
            if (ref == 1) ++totalPositive;
            else ++totalNegative;
        }

        if (totalPositive == 0 || totalNegative == 0) {
            setError("Reference raster must contain both positive (1) and negative (0) pixels");
            return false;
        }

        reportProgress(0.2, "Computing ROC curve...");

        // Compute ROC curve at multiple thresholds
        double step = (predMax - predMin) / numThresholds;
        struct RocPoint { double threshold; double TPR; double FPR; };
        std::vector<RocPoint> rocCurve;

        // Add point for threshold above max (all negative predictions)
        rocCurve.push_back({predMax + step, 0.0, 0.0});

        for (int t = 0; t <= numThresholds; ++t) {
            double threshold = predMax - t * step;
            int64_t TP = 0, FP = 0;

            for (int64_t i = 0; i < total; ++i) {
                if (hasNDP && dPred[i] == noDataP) continue;
                if (hasNDR && dRef[i] == noDataR) continue;
                int ref = static_cast<int>(std::round(dRef[i]));
                if (ref != 0 && ref != 1) continue;

                if (dPred[i] >= threshold) {
                    if (ref == 1) ++TP;
                    else ++FP;
                }
            }

            double TPR = static_cast<double>(TP) / totalPositive;
            double FPR = static_cast<double>(FP) / totalNegative;
            rocCurve.push_back({threshold, TPR, FPR});

            if (t % 10 == 0)
                reportProgress(0.2 + 0.6 * static_cast<double>(t) / numThresholds);
        }

        // Sort by FPR for trapezoidal integration
        std::sort(rocCurve.begin(), rocCurve.end(),
                  [](const RocPoint& a, const RocPoint& b) { return a.FPR < b.FPR; });

        // Compute AUC using trapezoidal rule
        double AUC = 0.0;
        for (size_t i = 1; i < rocCurve.size(); ++i) {
            double dFPR = rocCurve[i].FPR - rocCurve[i - 1].FPR;
            double avgTPR = (rocCurve[i].TPR + rocCurve[i - 1].TPR) / 2.0;
            AUC += dFPR * avgTPR;
        }

        reportProgress(0.9, "Writing report...");

        // Write report
        QString reportPath = parameter("output_report").toString();
        std::ofstream outFile(reportPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output report file: " + reportPath);
            return false;
        }

        outFile << "ROC Curve Analysis\n";
        outFile << "==================\n\n";
        outFile << "Prediction: " << parameter("prediction").toString().toStdString() << "\n";
        outFile << "Reference: " << parameter("reference").toString().toStdString() << "\n\n";
        outFile << "Total positive pixels: " << totalPositive << "\n";
        outFile << "Total negative pixels: " << totalNegative << "\n";
        outFile << "Number of thresholds: " << numThresholds << "\n\n";
        outFile << "Area Under the Curve (AUC): " << AUC << "\n\n";

        outFile << "Interpretation:\n";
        outFile << "  AUC = 1.0: Perfect discrimination\n";
        outFile << "  AUC = 0.5: No discrimination (random)\n";
        outFile << "  AUC < 0.5: Worse than random\n\n";

        outFile << "ROC Curve Data:\n";
        outFile << "Threshold\tFPR\tTPR\n";
        for (const auto& pt : rocCurve) {
            outFile << pt.threshold << "\t" << pt.FPR << "\t" << pt.TPR << "\n";
        }

        outFile.close();

        // Build ROC line chart
        ChartResult chart;
        chart.type = ChartResult::Line;
        chart.title = QString("ROC Curve (AUC = %1)").arg(AUC, 0, 'f', 4);
        chart.xLabel = "False Positive Rate";
        chart.yLabel = "True Positive Rate";
        ChartSeries s;
        s.label = "ROC";
        s.color = ChartColor(255, 0, 0);
        for (const auto& pt : rocCurve) {
            s.x.push_back(pt.FPR);
            s.y.push_back(pt.TPR);
        }
        chart.series.push_back(s);
        // Add diagonal reference line
        ChartSeries diag;
        diag.label = "Random";
        diag.color = ChartColor(128, 128, 128);
        diag.x = {0.0, 1.0};
        diag.y = {0.0, 1.0};
        chart.series.push_back(diag);
        setChartResult(chart);

        reportProgress(1.0,
            QString("AUC = %1, Positives = %2, Negatives = %3")
                .arg(AUC, 0, 'f', 6)
                .arg(totalPositive)
                .arg(totalNegative));

        return true;
    }
};

REGISTER_MODULE(RocModule)

} // namespace aplaceholder
