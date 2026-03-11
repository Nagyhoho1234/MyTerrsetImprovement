#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <cmath>

namespace aplaceholder {

class ChangeAnalysisModule : public Module {
public:
    QString name() const override { return "CHANGE_ANALYSIS"; }
    QString description() const override {
        return "Analyzes land cover change between two time periods. "
               "Identifies gains, losses, persistence, and net change for each category.";
    }
    QString category() const override { return "Land Change Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("earlier_map", "Earlier land cover image"),
            ParameterDef::file("later_map", "Later land cover image"),
            ParameterDef::output("output", "Output change analysis"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("earlier_map").toString());
        auto r2 = GdalIO::read(parameter("later_map").toString());
        if (!r1 || !r2) {
            setError("Failed to read input rasters");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        if (cols != r2->cols() || rows != r2->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        const auto& d1 = r1->data(0);
        const auto& d2 = r2->data(0);
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Collect all unique classes
        std::set<int> classes;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (d1[i] == noData || d2[i] == noData))
                continue;
            classes.insert(static_cast<int>(d1[i]));
            classes.insert(static_cast<int>(d2[i]));
        }

        if (classes.empty()) {
            setError("No valid pixels found in the input rasters");
            return false;
        }

        // Build transition count matrix: transitions[from_class][to_class] = count
        std::map<int, std::map<int, int64_t>> transitions;
        for (int c1 : classes)
            for (int c2 : classes)
                transitions[c1][c2] = 0;

        reportProgress(0.0, "Counting transitions...");
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (d1[i] == noData || d2[i] == noData))
                continue;
            int from = static_cast<int>(d1[i]);
            int to = static_cast<int>(d2[i]);
            transitions[from][to]++;

            if (i % 1000000 == 0)
                reportProgress(0.3 * static_cast<double>(i) / total);
        }

        // Compute per-class statistics
        // area_t1[c] = total pixels of class c at time 1
        // area_t2[c] = total pixels of class c at time 2
        // gains[c] = pixels that became class c (were something else at t1)
        // losses[c] = pixels that left class c (were class c at t1 but not at t2)
        // persistence[c] = pixels that stayed as class c
        // net_change[c] = gains - losses
        std::map<int, int64_t> area_t1, area_t2, gains, losses, persistence, net_change;

        for (int c : classes) {
            int64_t a1 = 0, a2 = 0, g = 0, l = 0, p = 0;

            // area at t1: sum of row c
            for (int c2 : classes)
                a1 += transitions[c][c2];

            // area at t2: sum of column c
            for (int c1 : classes)
                a2 += transitions[c1][c];

            // persistence: diagonal
            p = transitions[c][c];

            // gains: pixels that became c from other classes
            for (int c1 : classes) {
                if (c1 != c)
                    g += transitions[c1][c];
            }

            // losses: pixels that left c to other classes
            for (int c2 : classes) {
                if (c2 != c)
                    l += transitions[c][c2];
            }

            area_t1[c] = a1;
            area_t2[c] = a2;
            gains[c] = g;
            losses[c] = l;
            persistence[c] = p;
            net_change[c] = g - l;
        }

        reportProgress(0.5, "Building transition map...");

        // Create output transition raster: pixel value = from_class * 100 + to_class
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(r1->geoTransform());
        output.setProjection(r1->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (d1[i] == noData || d2[i] == noData)) {
                out[i] = noData;
                continue;
            }
            int from = static_cast<int>(d1[i]);
            int to = static_cast<int>(d2[i]);
            out[i] = static_cast<double>(from * 100 + to);

            if (i % 1000000 == 0)
                reportProgress(0.5 + 0.3 * static_cast<double>(i) / total);
        }

        reportProgress(0.8, "Writing output raster...");
        QString outputPath = parameter("output").toString();
        if (!GdalIO::write(output, outputPath)) {
            setError("Failed to write output raster");
            return false;
        }

        // Write CSV summary file
        reportProgress(0.9, "Writing statistics CSV...");
        QString csvPath = outputPath;
        if (csvPath.endsWith(".tif", Qt::CaseInsensitive) ||
            csvPath.endsWith(".rst", Qt::CaseInsensitive)) {
            csvPath = csvPath.left(csvPath.lastIndexOf('.'));
        }
        csvPath += "_stats.csv";

        std::ofstream csv(csvPath.toStdString());
        if (!csv.is_open()) {
            setError("Failed to write CSV summary file");
            return false;
        }

        csv << "Class,Area_T1,Area_T2,Gains,Losses,Net_Change,Persistence\n";
        for (int c : classes) {
            csv << c << ","
                << area_t1[c] << ","
                << area_t2[c] << ","
                << gains[c] << ","
                << losses[c] << ","
                << net_change[c] << ","
                << persistence[c] << "\n";
        }

        // Write full cross-tabulation matrix
        csv << "\nTransition Matrix (from\\to)";
        for (int c : classes)
            csv << "," << c;
        csv << "\n";
        for (int c1 : classes) {
            csv << c1;
            for (int c2 : classes)
                csv << "," << transitions[c1][c2];
            csv << "\n";
        }

        csv.close();

        reportProgress(1.0, "Change analysis complete.");
        return true;
    }
};

REGISTER_MODULE(ChangeAnalysisModule)

} // namespace aplaceholder
