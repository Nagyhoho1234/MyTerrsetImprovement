#include "Module.h"
#include "ModuleRegistry.h"
#include <QFile>
#include <QTextStream>
#include <QRegularExpression>

namespace aplaceholder {

class Avl2CsvModule : public Module {
public:
    QString name() const override { return "AVL2CSV"; }
    QString description() const override {
        return "Convert a TerrSet attribute values (.avl) file to CSV format.";
    }
    QString category() const override { return "Import/Export"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input AVL file"),
            ParameterDef::output("output", "Output CSV file"),
            ParameterDef::boolean("write_header", "Write header row", true),
        };
    }

    bool execute() override {
        QFile inFile(parameter("input").toString());
        if (!inFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open input AVL file");
            return false;
        }

        QTextStream in(&inFile);
        QFile outFile(parameter("output").toString());
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file");
            return false;
        }
        QTextStream out(&outFile);

        if (parameter("write_header").toBool())
            out << "ID,Value\n";

        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty()) continue;
            QStringList fields = line.split(QRegularExpression("\\s+"));
            out << fields.join(',') << "\n";
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(Avl2CsvModule)

} // namespace aplaceholder
