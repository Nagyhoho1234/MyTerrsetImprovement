#include "Module.h"
#include "ModuleRegistry.h"
#include <QFile>
#include <QTextStream>

namespace aplaceholder {

class Var2FixModule : public Module {
public:
    QString name() const override { return "VAR2FIX"; }
    QString description() const override {
        return "Convert variable-length ASCII records to fixed-length format.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input variable-length file"),
            ParameterDef::output("output", "Output fixed-length file"),
            ParameterDef::integer("record_length", "Fixed record length", 80, 1, 99999,
                "Target fixed record length in characters"),
        };
    }

    bool execute() override {
        QFile inFile(parameter("input").toString());
        if (!inFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open input file");
            return false;
        }

        QTextStream in(&inFile);
        int recLen = parameter("record_length").toInt();

        QFile outFile(parameter("output").toString());
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file");
            return false;
        }
        QTextStream out(&outFile);

        while (!in.atEnd()) {
            QString line = in.readLine();
            // Pad or truncate to fixed length
            if (line.length() < recLen)
                line = line.leftJustified(recLen, ' ');
            else
                line = line.left(recLen);
            out << line << "\n";
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(Var2FixModule)

} // namespace aplaceholder
