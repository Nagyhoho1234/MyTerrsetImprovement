#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <cmath>

namespace aplaceholder {

class JaccardModule : public Module {
public:
    QString name() const override { return "JACCARD"; }
    QString description() const override {
        return "Jaccard similarity index between two binary rasters. Computes "
               "J = |intersection| / |union|, where intersection is the count of pixels "
               "with value=1 in both rasters and union is the count of pixels with value=1 "
               "in either raster. Produces a text report with the similarity coefficient.";
    }
    QString category() const override { return "Habitat & Biodiversity"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First binary raster",
                "Binary raster (0/1) representing species presence or habitat"),
            ParameterDef::file("input2", "Second binary raster",
                "Binary raster (0/1) representing species presence or habitat"),
            ParameterDef::output("output_report", "Output report file (.txt)",
                "Text file containing Jaccard similarity results"),
        };
    }

    bool execute() override {
        // 1. Load both binary rasters
        auto raster1 = GdalIO::read(parameter("input1").toString());
        if (!raster1) {
            setError("Failed to read first binary raster.");
            return false;
        }

        auto raster2 = GdalIO::read(parameter("input2").toString());
        if (!raster2) {
            setError("Failed to read second binary raster.");
            return false;
        }

        reportProgress(0.10, "Input rasters loaded.");

        // 2. Validate dimensions
        int cols = raster1->cols();
        int rows = raster1->rows();
        int64_t total = raster1->cellCount();

        if (raster2->cols() != cols || raster2->rows() != rows) {
            setError(QString("Dimension mismatch: input1 is %1x%2 but input2 is %3x%4")
                         .arg(cols).arg(rows)
                         .arg(raster2->cols()).arg(raster2->rows()));
            return false;
        }

        // 3. Compute intersection and union
        const auto& d1 = raster1->data(0);
        const auto& d2 = raster2->data(0);

        double noData1 = raster1->noDataValue();
        double noData2 = raster2->noDataValue();
        bool hasND1 = raster1->hasNoData();
        bool hasND2 = raster2->hasNoData();

        int64_t countBoth = 0;      // intersection (both = 1)
        int64_t countEither = 0;    // union (at least one = 1)
        int64_t count1only = 0;     // present in raster1 only
        int64_t count2only = 0;     // present in raster2 only
        int64_t totalValid = 0;

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            if (hasND1 && d1[i] == noData1) skip = true;
            if (hasND2 && d2[i] == noData2) skip = true;
            if (skip) continue;

            bool p1 = (std::round(d1[i]) != 0);
            bool p2 = (std::round(d2[i]) != 0);

            if (p1 && p2)  ++countBoth;
            if (p1 || p2)  ++countEither;
            if (p1 && !p2) ++count1only;
            if (!p1 && p2) ++count2only;
            ++totalValid;

            if (i % 500000 == 0)
                reportProgress(0.10 + 0.70 * static_cast<double>(i) / total);
        }

        // 4. Compute Jaccard index
        double jaccard = (countEither > 0)
                             ? static_cast<double>(countBoth) / static_cast<double>(countEither)
                             : 0.0;

        reportProgress(0.85, "Writing report...");

        // 5. Write report
        QString reportPath = parameter("output_report").toString();
        QFile reportFile(reportPath);
        if (!reportFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to write output report: " + reportPath);
            return false;
        }

        QTextStream out(&reportFile);
        out << "Jaccard Similarity Analysis Report\n";
        out << "==================================\n\n";
        out << "Input 1: " << parameter("input1").toString() << "\n";
        out << "Input 2: " << parameter("input2").toString() << "\n\n";
        out << "Results:\n";
        out << "  Total valid pixels:          " << totalValid << "\n";
        out << "  Intersection (both present): " << countBoth << "\n";
        out << "  Union (either present):      " << countEither << "\n";
        out << "  Input 1 only:                " << count1only << "\n";
        out << "  Input 2 only:                " << count2only << "\n\n";
        out << "  Jaccard Similarity Index:     " << QString::number(jaccard, 'f', 6) << "\n";
        out << "  (0 = no overlap, 1 = identical)\n";

        reportProgress(1.0,
            QString("Done. Jaccard index: %1 (intersection=%2, union=%3)")
                .arg(jaccard, 0, 'f', 6)
                .arg(countBoth)
                .arg(countEither));

        return true;
    }
};

REGISTER_MODULE(JaccardModule)

} // namespace aplaceholder
