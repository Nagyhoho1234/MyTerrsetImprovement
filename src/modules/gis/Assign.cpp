#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <QFile>
#include <QTextStream>
#include <QMap>
#include <cmath>

namespace aplaceholder {

class AssignModule : public Module {
public:
    QString name() const override { return "ASSIGN"; }
    QString description() const override {
        return "Assigns new attribute values to a raster feature definition image "
               "based on an external attribute values file (CSV). Values not specified "
               "in the file are assigned zero in the output.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Feature definition image (integer class IDs)"),
            ParameterDef::file("attribute_file", "Attribute file (CSV with 'id' and 'value' columns)",
                "CSV file mapping class IDs to new attribute values"),
            ParameterDef::output("output", "Output image"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) {
            setError("Failed to read input raster");
            return false;
        }

        // Read attribute mapping from CSV
        QMap<int, double> mapping;
        if (!readAttributeFile(parameter("attribute_file").toString(), mapping)) {
            return false;
        }

        int cols = input->cols(), rows = input->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(input->noDataValue());

        const auto& inData = input->data(0);
        auto& out = output.data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();

        reportProgress(0.0, "Assigning attribute values...");

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && inData[i] == noData) {
                out[i] = noData;
                continue;
            }

            int classId = static_cast<int>(std::round(inData[i]));
            auto it = mapping.find(classId);
            if (it != mapping.end()) {
                out[i] = it.value();
            } else {
                // Values not in the attribute file are assigned zero
                out[i] = 0.0;
            }

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }

private:
    bool readAttributeFile(const QString& path, QMap<int, double>& mapping) {
        QFile file(path);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open attribute file: " + path);
            return false;
        }

        QTextStream in(&file);
        int lineNum = 0;
        int idCol = -1, valueCol = -1;

        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty()) continue;
            ++lineNum;

            QStringList parts = line.split(",");

            if (lineNum == 1) {
                // Parse header to find 'id' and 'value' columns
                for (int c = 0; c < parts.size(); ++c) {
                    QString col = parts[c].trimmed().toLower();
                    if (col == "id") idCol = c;
                    else if (col == "value") valueCol = c;
                }
                if (idCol < 0 || valueCol < 0) {
                    setError("Attribute CSV must have 'id' and 'value' columns in header");
                    return false;
                }
                continue;
            }

            if (parts.size() <= std::max(idCol, valueCol)) {
                setError("Insufficient columns at line " + QString::number(lineNum));
                return false;
            }

            bool okId, okVal;
            int id = parts[idCol].trimmed().toInt(&okId);
            double val = parts[valueCol].trimmed().toDouble(&okVal);

            if (!okId || !okVal) {
                setError("Invalid numeric values at line " + QString::number(lineNum));
                return false;
            }

            mapping[id] = val;
        }

        if (mapping.isEmpty()) {
            setError("No valid mappings found in attribute file");
            return false;
        }

        return true;
    }
};

REGISTER_MODULE(AssignModule)

} // namespace aplaceholder
