#include "Module.h"
#include "ModuleRegistry.h"
#include <QFile>
#include <QByteArray>

namespace aplaceholder {

class CrlfModule : public Module {
public:
    QString name() const override { return "CRLF"; }
    QString description() const override {
        return "Convert line endings between Unix (LF) and Windows (CR+LF).";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input file"),
            ParameterDef::output("output", "Output file"),
            ParameterDef::combo("mode", "Conversion mode",
                {"Unix to Windows (LF to CRLF)", "Windows to Unix (CRLF to LF)"}, 0),
        };
    }

    bool execute() override {
        QFile inFile(parameter("input").toString());
        if (!inFile.open(QIODevice::ReadOnly)) {
            setError("Failed to open input file");
            return false;
        }
        QByteArray data = inFile.readAll();
        inFile.close();

        int mode = parameter("mode").toInt();
        if (mode == 0) {
            // LF -> CRLF: first normalize to LF, then convert
            data.replace("\r\n", "\n");
            data.replace("\n", "\r\n");
        } else {
            // CRLF -> LF
            data.replace("\r\n", "\n");
        }

        reportProgress(0.5, "Writing output...");
        QFile outFile(parameter("output").toString());
        if (!outFile.open(QIODevice::WriteOnly)) {
            setError("Failed to open output file");
            return false;
        }
        outFile.write(data);
        outFile.close();
        return true;
    }
};

REGISTER_MODULE(CrlfModule)

} // namespace aplaceholder
