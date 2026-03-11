#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace aplaceholder {

class ErrMatModule : public Module {
public:
    QString name() const override { return "ERRMAT"; }
    QString description() const override {
        return "Error Matrix / Accuracy Assessment. Compares a classified raster "
               "against a reference (ground truth) raster and produces a confusion "
               "matrix with overall accuracy, producer's accuracy, user's accuracy, "
               "and Kappa coefficient.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("classified", "Classified image",
                "The categorical raster map to evaluate"),
            ParameterDef::file("reference", "Reference (ground truth) image",
                "Raster with true category values (0 = background/no data)"),
            ParameterDef::output("output", "Output report file",
                "Text file for the accuracy assessment report"),
        };
    }

    bool execute() override {
        auto classified = GdalIO::read(parameter("classified").toString());
        auto reference = GdalIO::read(parameter("reference").toString());
        if (!classified || !reference) {
            setError("Failed to read input rasters");
            return false;
        }

        if (classified->cols() != reference->cols() ||
            classified->rows() != reference->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = classified->cols(), rows = classified->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& dClass = classified->data(0);
        const auto& dRef = reference->data(0);
        bool hasNDClass = classified->hasNoData();
        bool hasNDRef = reference->hasNoData();
        double ndClass = classified->noDataValue();
        double ndRef = reference->noDataValue();

        // Build confusion matrix from pixels where reference != 0
        // (reference value 0 means no sample point)
        std::set<int> allClasses;
        std::map<std::pair<int,int>, int64_t> confMatrix; // (mapped, truth) -> count
        int64_t sampleCount = 0;

        reportProgress(0.0, "Building confusion matrix...");
        for (int64_t i = 0; i < total; ++i) {
            if (hasNDClass && dClass[i] == ndClass) continue;
            if (hasNDRef && dRef[i] == ndRef) continue;

            int refVal = static_cast<int>(dRef[i]);
            if (refVal == 0) continue; // background in truth image

            int classVal = static_cast<int>(dClass[i]);
            allClasses.insert(classVal);
            allClasses.insert(refVal);
            confMatrix[{classVal, refVal}]++;
            sampleCount++;

            if (i % 1000000 == 0)
                reportProgress(0.4 * static_cast<double>(i) / total);
        }

        if (sampleCount == 0) {
            setError("No valid sample pixels found (reference image has no non-zero values)");
            return false;
        }

        reportProgress(0.5, "Computing accuracy metrics...");

        // Convert classes to a sorted vector
        std::vector<int> classList(allClasses.begin(), allClasses.end());
        int numClasses = static_cast<int>(classList.size());

        // Compute row totals (mapped class totals) and column totals (truth class totals)
        std::map<int, int64_t> rowTotal; // sum across columns for each mapped class
        std::map<int, int64_t> colTotal; // sum across rows for each truth class
        int64_t diagonal = 0;

        for (int mc : classList) {
            for (int tc : classList) {
                auto it = confMatrix.find({mc, tc});
                int64_t count = (it != confMatrix.end()) ? it->second : 0;
                rowTotal[mc] += count;
                colTotal[tc] += count;
                if (mc == tc) diagonal += count;
            }
        }

        // Overall accuracy: Po = diagonal / total samples
        double Po = static_cast<double>(diagonal) / sampleCount;

        // Expected agreement: Pe = sum( (row_i_total * col_i_total) ) / N^2
        double Pe = 0.0;
        for (int c : classList) {
            Pe += static_cast<double>(rowTotal[c]) * static_cast<double>(colTotal[c]);
        }
        Pe /= (static_cast<double>(sampleCount) * sampleCount);

        // Kappa = (Po - Pe) / (1 - Pe)
        double kappa = 0.0;
        if (Pe < 1.0)
            kappa = (Po - Pe) / (1.0 - Pe);

        // Producer's accuracy per class: diagonal_i / col_total_i
        // (1 - omission error)
        std::map<int, double> producerAcc;
        for (int c : classList) {
            auto it = confMatrix.find({c, c});
            int64_t diag = (it != confMatrix.end()) ? it->second : 0;
            if (colTotal[c] > 0)
                producerAcc[c] = static_cast<double>(diag) / colTotal[c];
            else
                producerAcc[c] = 0.0;
        }

        // User's accuracy per class: diagonal_i / row_total_i
        // (1 - commission error)
        std::map<int, double> userAcc;
        for (int c : classList) {
            auto it = confMatrix.find({c, c});
            int64_t diag = (it != confMatrix.end()) ? it->second : 0;
            if (rowTotal[c] > 0)
                userAcc[c] = static_cast<double>(diag) / rowTotal[c];
            else
                userAcc[c] = 0.0;
        }

        // Per-category Kappa
        std::map<int, double> perClassKappa;
        for (int c : classList) {
            auto it = confMatrix.find({c, c});
            int64_t diag = (it != confMatrix.end()) ? it->second : 0;
            double expected = static_cast<double>(rowTotal[c]) * colTotal[c]
                              / sampleCount;
            double maxPossible = 0.5 * (rowTotal[c] + colTotal[c]);
            if (maxPossible > expected)
                perClassKappa[c] = (static_cast<double>(diag) - expected)
                                   / (maxPossible - expected);
            else
                perClassKappa[c] = 0.0;
        }

        reportProgress(0.7, "Writing report...");

        // Write the report
        QString outPath = parameter("output").toString();
        std::ofstream report(outPath.toStdString());
        if (!report.is_open()) {
            setError("Failed to open output report file: " + outPath);
            return false;
        }

        report << "ERROR MATRIX / ACCURACY ASSESSMENT REPORT\n";
        report << "==========================================\n\n";
        report << "Classified image: " << parameter("classified").toString().toStdString() << "\n";
        report << "Reference image:  " << parameter("reference").toString().toStdString() << "\n";
        report << "Total sample points: " << sampleCount << "\n";
        report << "Number of classes: " << numClasses << "\n\n";

        // Determine column width
        int colWidth = 8;

        // Print confusion matrix header
        report << "CONFUSION MATRIX (rows=mapped, columns=truth)\n";
        report << std::string(70, '-') << "\n";
        report << std::setw(colWidth) << "Class";
        for (int tc : classList)
            report << std::setw(colWidth) << tc;
        report << std::setw(colWidth + 2) << "RowTotal";
        report << std::setw(colWidth + 2) << "UserAcc";
        report << "\n";
        report << std::string(70, '-') << "\n";

        for (int mc : classList) {
            report << std::setw(colWidth) << mc;
            for (int tc : classList) {
                auto it = confMatrix.find({mc, tc});
                int64_t count = (it != confMatrix.end()) ? it->second : 0;
                report << std::setw(colWidth) << count;
            }
            report << std::setw(colWidth + 2) << rowTotal[mc];
            report << std::setw(colWidth + 2) << std::fixed << std::setprecision(4)
                   << userAcc[mc];
            report << "\n";
        }

        report << std::string(70, '-') << "\n";

        // Column totals
        report << std::setw(colWidth) << "ColTotal";
        for (int tc : classList)
            report << std::setw(colWidth) << colTotal[tc];
        report << std::setw(colWidth + 2) << sampleCount;
        report << "\n";

        // Producer's accuracy row
        report << std::setw(colWidth) << "ProdAcc";
        for (int tc : classList)
            report << std::setw(colWidth) << std::fixed << std::setprecision(4)
                   << producerAcc[tc];
        report << "\n\n";

        // Summary statistics
        report << "SUMMARY STATISTICS\n";
        report << std::string(40, '-') << "\n";
        report << "Overall Accuracy:   " << std::fixed << std::setprecision(4)
               << Po << " (" << std::setprecision(2) << (Po * 100.0) << "%)\n";
        report << "Overall Kappa (KIA): " << std::fixed << std::setprecision(4)
               << kappa << "\n\n";

        // Per-class details
        report << "PER-CLASS STATISTICS\n";
        report << std::string(70, '-') << "\n";
        report << std::setw(colWidth) << "Class"
               << std::setw(14) << "ProdAccuracy"
               << std::setw(14) << "UserAccuracy"
               << std::setw(14) << "OmissionErr"
               << std::setw(14) << "CommissionErr"
               << std::setw(12) << "KIA"
               << "\n";
        report << std::string(70, '-') << "\n";

        for (int c : classList) {
            report << std::setw(colWidth) << c
                   << std::setw(14) << std::fixed << std::setprecision(4) << producerAcc[c]
                   << std::setw(14) << std::fixed << std::setprecision(4) << userAcc[c]
                   << std::setw(14) << std::fixed << std::setprecision(4) << (1.0 - producerAcc[c])
                   << std::setw(14) << std::fixed << std::setprecision(4) << (1.0 - userAcc[c])
                   << std::setw(12) << std::fixed << std::setprecision(4) << perClassKappa[c]
                   << "\n";
        }

        report << std::string(70, '-') << "\n";
        report.close();

        QString statsMsg = QString("Accuracy assessment complete. "
                                   "Overall accuracy: %1%, Kappa: %2")
                           .arg(Po * 100.0, 0, 'f', 2)
                           .arg(kappa, 0, 'f', 4);
        reportProgress(1.0, statsMsg);
        return true;
    }
};

REGISTER_MODULE(ErrMatModule)

} // namespace aplaceholder
