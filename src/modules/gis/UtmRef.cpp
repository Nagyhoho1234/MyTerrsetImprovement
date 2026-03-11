#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <QFile>
#include <QTextStream>

namespace aplaceholder {

class UtmRefModule : public Module {
public:
    QString name() const override { return "UTMREF"; }
    QString description() const override {
        return "Compute the appropriate UTM zone from latitude/longitude bounds "
               "and create a georeferencing file.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::real("min_lon", "Minimum longitude", -180.0, -180.0, 180.0,
                "Western bound in decimal degrees"),
            ParameterDef::real("max_lon", "Maximum longitude", 180.0, -180.0, 180.0,
                "Eastern bound in decimal degrees"),
            ParameterDef::real("min_lat", "Minimum latitude", -90.0, -90.0, 90.0,
                "Southern bound in decimal degrees"),
            ParameterDef::real("max_lat", "Maximum latitude", 90.0, -90.0, 90.0,
                "Northern bound in decimal degrees"),
            ParameterDef::output("output", "Output reference file"),
        };
    }

    bool execute() override {
        double minLon = parameter("min_lon").toDouble();
        double maxLon = parameter("max_lon").toDouble();
        double minLat = parameter("min_lat").toDouble();
        double maxLat = parameter("max_lat").toDouble();

        if (minLon >= maxLon) {
            setError("Minimum longitude must be less than maximum longitude");
            return false;
        }
        if (minLat >= maxLat) {
            setError("Minimum latitude must be less than maximum latitude");
            return false;
        }
        if (minLon < -180.0 || maxLon > 180.0) {
            setError("Longitude must be between -180 and 180 degrees");
            return false;
        }
        if (minLat < -90.0 || maxLat > 90.0) {
            setError("Latitude must be between -90 and 90 degrees");
            return false;
        }

        reportProgress(0.3, "Computing UTM zone...");

        // Compute the central meridian of the bounding box
        double centralLon = (minLon + maxLon) / 2.0;
        double centralLat = (minLat + maxLat) / 2.0;

        // UTM zone number: 1-60, each 6 degrees wide starting at -180
        int zone = static_cast<int>(std::floor((centralLon + 180.0) / 6.0)) + 1;
        if (zone < 1) zone = 1;
        if (zone > 60) zone = 60;

        // Determine hemisphere
        bool northern = (centralLat >= 0.0);

        // Compute the central meridian of the UTM zone
        double zoneCentralMeridian = (zone - 1) * 6.0 - 180.0 + 3.0;

        reportProgress(0.6, "Writing reference file...");

        QString outputPath = parameter("output").toString();
        QFile file(outputPath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output file for writing: " + outputPath);
            return false;
        }

        QTextStream out(&file);
        out << "ref. system  : UTM Zone " << zone
            << (northern ? "N" : "S") << "\n";
        out << "projection   : Transverse Mercator\n";
        out << "datum        : WGS 84\n";
        out << "zone         : " << zone << "\n";
        out << "hemisphere   : " << (northern ? "North" : "South") << "\n";
        out << "central merid: " << zoneCentralMeridian << "\n";
        out << "min. lon     : " << minLon << "\n";
        out << "max. lon     : " << maxLon << "\n";
        out << "min. lat     : " << minLat << "\n";
        out << "max. lat     : " << maxLat << "\n";
        file.close();

        reportProgress(1.0, "UTM Zone " + QString::number(zone)
                       + (northern ? "N" : "S") + " reference file created.");
        return true;
    }
};

REGISTER_MODULE(UtmRefModule)

} // namespace aplaceholder
