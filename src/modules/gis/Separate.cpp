#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class SeparateModule : public Module {
public:
    QString name() const override { return "SEPARATE"; }
    QString description() const override {
        return "Extract individual bands from a multi-band raster into separate files.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input multi-band raster"),
            ParameterDef({
                "output_prefix", "Output prefix",
                ParameterDef::String, {}, {}, 0, 0,
                "Creates prefix_band1.tif, prefix_band2.tif, etc.", true
            }),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int bands = raster->bands();
        QString prefix = parameter("output_prefix").toString();

        for (int b = 0; b < bands; ++b) {
            Raster single(raster->cols(), raster->rows(), 1, raster->dataType());
            single.setGeoTransform(raster->geoTransform());
            single.setProjection(raster->projection());
            if (raster->hasNoData())
                single.setNoDataValue(raster->noDataValue());

            single.data(0) = raster->data(b);

            QString outPath = QString("%1_band%2.tif").arg(prefix).arg(b + 1);
            if (!GdalIO::write(single, outPath)) {
                setError("Failed to write band " + QString::number(b + 1));
                return false;
            }

            reportProgress(static_cast<double>(b + 1) / bands);
        }

        reportProgress(1.0, "Done");
        return true;
    }
};

REGISTER_MODULE(SeparateModule)

} // namespace aplaceholder
