#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class MetaUpdateModule : public Module {
public:
    QString name() const override { return "METAUPDATE"; }
    QString description() const override {
        return "Update raster metadata: projection, nodata value, or cell size.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::separator("Optional metadata changes"),
            ParameterDef({
                "nodata", "NoData value", ParameterDef::Double, -9999.0,
                {}, -1e15, 1e15, "Set nodata value (leave default to skip)", false
            }),
            ParameterDef({
                "projection", "Projection (EPSG code)", ParameterDef::String, {},
                {}, 0, 0, "e.g. EPSG:4326 or EPSG:32617", false
            }),
            ParameterDef({
                "cell_size", "Cell size", ParameterDef::Double, 0.0,
                {}, 0.0, 1e10, "Set cell size (0 to keep original)", false
            }),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        // Update nodata if provided
        QVariant noDataParam = parameter("nodata");
        if (noDataParam.isValid() && !noDataParam.isNull()) {
            raster->setNoDataValue(noDataParam.toDouble());
        }

        // Update projection if provided
        QString proj = parameter("projection").toString();
        if (!proj.isEmpty()) {
            raster->setProjection(proj);
        }

        // Update cell size if provided
        double cellSize = parameter("cell_size").toDouble();
        if (cellSize > 0.0) {
            GeoTransform gt = raster->geoTransform();
            gt.pixelWidth = cellSize;
            gt.pixelHeight = cellSize;
            raster->setGeoTransform(gt);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(*raster, parameter("output").toString());
    }
};

REGISTER_MODULE(MetaUpdateModule)

} // namespace aplaceholder
