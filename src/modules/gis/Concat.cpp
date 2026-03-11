#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

namespace aplaceholder {

class ConcatModule : public Module {
public:
    QString name() const override { return "CONCAT"; }
    QString description() const override {
        return "Concatenate/stack multiple single-band rasters into one multi-band raster.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef({
                "inputs", "Input rasters (comma-separated paths)",
                ParameterDef::String, {}, {}, 0, 0,
                "Comma-separated list of input raster file paths", true
            }),
            ParameterDef::output("output", "Output multi-band raster"),
        };
    }

    bool execute() override {
        QString inputsStr = parameter("inputs").toString();
        QStringList paths = inputsStr.split(',', Qt::SkipEmptyParts);

        if (paths.isEmpty()) {
            setError("No input rasters specified");
            return false;
        }

        // Read first raster to get dimensions
        auto first = GdalIO::read(paths[0].trimmed());
        if (!first) {
            setError("Failed to read first input raster: " + paths[0].trimmed());
            return false;
        }

        int cols = first->cols(), rows = first->rows();
        int totalBands = first->bands();

        // Count total bands across all inputs
        std::vector<std::unique_ptr<Raster>> rasters;
        rasters.push_back(std::move(first));

        for (int i = 1; i < paths.size(); ++i) {
            auto r = GdalIO::read(paths[i].trimmed());
            if (!r) {
                setError("Failed to read input raster: " + paths[i].trimmed());
                return false;
            }
            if (r->cols() != cols || r->rows() != rows) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
            totalBands += r->bands();
            rasters.push_back(std::move(r));
        }

        Raster output(cols, rows, totalBands, DataType::Float64);
        output.setGeoTransform(rasters[0]->geoTransform());
        output.setProjection(rasters[0]->projection());
        if (rasters[0]->hasNoData())
            output.setNoDataValue(rasters[0]->noDataValue());

        int bandIdx = 0;
        for (int r = 0; r < static_cast<int>(rasters.size()); ++r) {
            for (int b = 0; b < rasters[r]->bands(); ++b) {
                output.data(bandIdx) = rasters[r]->data(b);
                ++bandIdx;
            }
            reportProgress(static_cast<double>(r + 1) / rasters.size());
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ConcatModule)

} // namespace aplaceholder
