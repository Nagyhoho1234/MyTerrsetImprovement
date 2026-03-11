#include "Module.h"
#include "ModuleRegistry.h"
#include <QFile>
#include <QTextStream>

namespace aplaceholder {

class CreateTsfModule : public Module {
public:
    QString name() const override { return "CREATE_TSF"; }
    QString description() const override {
        return "Create Time Series File (.tsf) with metadata for temporal analysis.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::output("output_path", "Output .tsf file path"),
            ParameterDef{"start_date", "Start date (YYYY-MM-DD)", ParameterDef::String, QVariant(""), {}, 0, 0, "", true},
            ParameterDef::real("interval", "Time interval between images (days)", 1.0),
            ParameterDef::integer("num_images", "Number of images in series", 1),
        };
    }

    bool execute() override {
        QString path = parameter("output_path").toString();
        QString startDate = parameter("start_date").toString();
        double interval = parameter("interval").toDouble();
        int numImages = parameter("num_images").toInt();

        if (numImages <= 0) {
            setError("Number of images must be positive");
            return false;
        }

        QFile file(path);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to create TSF file");
            return false;
        }

        QTextStream out(&file);
        out << "time series file\n";
        out << "start_date: " << startDate << "\n";
        out << "interval: " << interval << "\n";
        out << "num_images: " << numImages << "\n";
        file.close();

        reportProgress(1.0, "TSF file created.");
        return true;
    }
};

REGISTER_MODULE(CreateTsfModule)

} // namespace aplaceholder
