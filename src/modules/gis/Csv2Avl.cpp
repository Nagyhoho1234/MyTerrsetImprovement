#include "Module.h"
#include "ModuleRegistry.h"
#include <QFile>
#include <QTextStream>

namespace aplaceholder {

class Csv2AvlModule : public Module {
public:
    QString name() const override { return "CSV2AVL"; }
    QString description() const override {
        return "Convert a CSV file to TerrSet attribute values (.avl) format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input CSV file"),
            ParameterDef::output("output", "Output AVL file"),
            ParameterDef::integer("id_column", "ID column index", 0, 0, 999,
                "Zero-based column index for feature IDs"),
            ParameterDef::integer("value_column", "Value column index", 1, 0, 999,
                "Zero-based column index for attribute values"),
            ParameterDef::boolean("has_header", "File has header row", true),
        };
    }

    bool execute() override {
        QFile inFile(parameter("input").toString());
        if (!inFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open input CSV file");
            return false;
        }

        QTextStream in(&inFile);
        int idCol = parameter("id_column").toInt();
        int valCol = parameter("value_column").toInt();
        bool hasHeader = parameter("has_header").toBool();

        QFile outFile(parameter("output").toString());
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file");
            return false;
        }
        QTextStream out(&outFile);

        if (hasHeader && !in.atEnd())
            in.readLine(); // skip header

        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty()) continue;
            QStringList fields = line.split(',');
            if (idCol < fields.size() && valCol < fields.size()) {
                out << fields[idCol].trimmed() << " "
                    << fields[valCol].trimmed() << "\n";
            }
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(Csv2AvlModule)

} // namespace aplaceholder
