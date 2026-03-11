#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <string>

namespace aplaceholder {

class ReclassModule : public Module {
public:
    QString name() const override { return "RECLASS"; }
    QString description() const override {
        return "Reclassification of raster values by intervals. "
               "Assigns new class values based on equal intervals or "
               "user-defined breakpoints.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::output("output", "Output reclassified image"),
            ParameterDef::combo("method", "Classification method",
                {"Equal Interval", "User Defined"}, 0,
                "Method for determining class breakpoints"),
            ParameterDef::integer("num_classes", "Number of classes", 5, 2, 256,
                "Number of output classes (Equal Interval mode)"),
            ParameterDef::file("reclass_file", "Reclass definition file"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input raster"); return false; }

        int cols = input->cols(), rows = input->rows();
        int method = parameter("method").toInt();
        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        if (method == 0) {
            // Equal Interval
            int numClasses = parameter("num_classes").toInt();

            // Compute min and max of valid pixels
            double minVal = std::numeric_limits<double>::max();
            double maxVal = std::numeric_limits<double>::lowest();
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && data[i] == noData) continue;
                if (data[i] < minVal) minVal = data[i];
                if (data[i] > maxVal) maxVal = data[i];
            }

            if (minVal > maxVal) {
                setError("No valid pixels found in input raster");
                return false;
            }

            double interval = (maxVal - minVal) / numClasses;

            reportProgress(0.0, "Reclassifying (Equal Interval)...");
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && data[i] == noData) {
                    out[i] = noData;
                    continue;
                }

                // Determine which class this value falls into
                int classId;
                if (interval <= 0) {
                    classId = 1;
                } else {
                    classId = static_cast<int>((data[i] - minVal) / interval) + 1;
                    // Clamp: the maximum value should fall into the last class
                    if (classId > numClasses) classId = numClasses;
                    if (classId < 1) classId = 1;
                }
                out[i] = static_cast<double>(classId);

                if (i % 1000000 == 0)
                    reportProgress(static_cast<double>(i) / total);
            }
        } else {
            // User Defined: read reclass file
            QString reclassPath = parameter("reclass_file").toString();
            if (reclassPath.isEmpty()) {
                setError("Reclass definition file is required for User Defined method");
                return false;
            }

            // Parse reclass file: each line is "old_low old_high new_value"
            struct ReclassRule {
                double oldLow;
                double oldHigh;
                double newValue;
            };
            std::vector<ReclassRule> rules;

            std::ifstream fin(reclassPath.toStdString());
            if (!fin.is_open()) {
                setError("Failed to open reclass definition file: " + reclassPath);
                return false;
            }

            std::string line;
            while (std::getline(fin, line)) {
                // Skip empty lines and comments
                if (line.empty() || line[0] == '#') continue;
                double lo, hi, nv;
                if (sscanf(line.c_str(), "%lf %lf %lf", &lo, &hi, &nv) == 3) {
                    rules.push_back({lo, hi, nv});
                }
            }
            fin.close();

            if (rules.empty()) {
                setError("No valid rules found in reclass definition file");
                return false;
            }

            reportProgress(0.0, "Reclassifying (User Defined)...");
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && data[i] == noData) {
                    out[i] = noData;
                    continue;
                }

                double val = data[i];
                bool matched = false;
                for (const auto& rule : rules) {
                    if (val >= rule.oldLow && val <= rule.oldHigh) {
                        out[i] = rule.newValue;
                        matched = true;
                        break;
                    }
                }
                if (!matched) {
                    out[i] = noData;
                }

                if (i % 1000000 == 0)
                    reportProgress(static_cast<double>(i) / total);
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ReclassModule)

} // namespace aplaceholder
