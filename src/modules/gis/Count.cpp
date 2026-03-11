#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <map>
#include <fstream>

namespace aplaceholder {

class CountModule : public Module {
public:
    QString name() const override { return "COUNT"; }
    QString description() const override {
        return "Count pixels by category in a categorical raster. "
               "Produces a summary table of category values and their pixel counts, "
               "along with proportional area for each category.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input categorical raster image"),
            ParameterDef::output("output", "Output count statistics file"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        reportProgress(0.0, "Counting pixels by category...");

        std::map<int, int64_t> counts;
        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data[i] == noData) continue;
            int cat = static_cast<int>(std::round(data[i]));
            counts[cat]++;
            validCount++;
        }

        reportProgress(0.8, "Writing results...");

        QString outputPath = parameter("output").toString();
        std::ofstream outFile(outputPath.toStdString());
        if (!outFile.is_open()) {
            setError("Failed to open output file: " + outputPath);
            return false;
        }

        outFile << "Pixel Count by Category\n";
        outFile << "=======================\n\n";
        outFile << "Input: " << parameter("input").toString().toStdString() << "\n";
        outFile << "Total valid pixels: " << validCount << "\n";
        outFile << "Number of categories: " << counts.size() << "\n\n";
        outFile << "Category\tCount\tProportion\n";
        outFile << "--------\t-----\t----------\n";

        for (const auto& kv : counts) {
            double prop = static_cast<double>(kv.second) / validCount;
            outFile << kv.first << "\t" << kv.second << "\t" << prop << "\n";
        }

        outFile.close();

        reportProgress(1.0,
            QString("Found %1 categories across %2 valid pixels")
                .arg(counts.size()).arg(validCount));
        return true;
    }
};

REGISTER_MODULE(CountModule)

} // namespace aplaceholder
