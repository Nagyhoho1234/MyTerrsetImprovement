#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class UpdateModule : public Module {
public:
    QString name() const override { return "UPDATE"; }
    QString description() const override {
        return "Update raster values where a condition is met, replacing with "
               "values from an update raster or a constant.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::file("update", "Update raster (or leave empty for constant)"),
            ParameterDef::real("update_value", "Update constant value", 0.0, -999999.0, 999999.0,
                "Constant value used when no update raster is provided"),
            ParameterDef::file("condition", "Condition raster (optional)"),
            ParameterDef::output("output", "Output raster"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) {
            setError("Failed to read input raster");
            return false;
        }

        int cols = input->cols(), rows = input->rows(), bands = input->bands();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Read update raster if provided
        QString updatePath = parameter("update").toString();
        std::unique_ptr<Raster> updateRaster;
        bool useConstant = updatePath.isEmpty();
        double updateConst = parameter("update_value").toDouble();

        if (!useConstant) {
            updateRaster = GdalIO::read(updatePath);
            if (!updateRaster) {
                setError("Failed to read update raster");
                return false;
            }
            if (updateRaster->cols() != cols || updateRaster->rows() != rows) {
                setError("Update raster dimensions do not match input raster");
                return false;
            }
        }

        // Read condition raster if provided
        QString condPath = parameter("condition").toString();
        std::unique_ptr<Raster> condRaster;
        bool hasCondition = !condPath.isEmpty();

        if (hasCondition) {
            condRaster = GdalIO::read(condPath);
            if (!condRaster) {
                setError("Failed to read condition raster");
                return false;
            }
            if (condRaster->cols() != cols || condRaster->rows() != rows) {
                setError("Condition raster dimensions do not match input raster");
                return false;
            }
        }

        Raster output(cols, rows, bands, input->dataType());
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        if (input->hasNoData())
            output.setNoDataValue(input->noDataValue());

        for (int b = 0; b < bands; ++b) {
            const auto& src = input->data(b);
            auto& dst = output.data(b);

            for (int64_t i = 0; i < total; ++i) {
                bool condMet = true;
                if (hasCondition) {
                    double cv = condRaster->data(0)[i];
                    condMet = (cv != 0.0 && !std::isnan(cv));
                }

                if (condMet) {
                    dst[i] = useConstant ? updateConst
                                         : updateRaster->data(b < updateRaster->bands() ? b : 0)[i];
                } else {
                    dst[i] = src[i];
                }
            }

            reportProgress(static_cast<double>(b + 1) / bands);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(UpdateModule)

} // namespace aplaceholder
