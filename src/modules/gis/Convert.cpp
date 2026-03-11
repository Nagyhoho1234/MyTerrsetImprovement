#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class ConvertModule : public Module {
public:
    QString name() const override { return "CONVERT"; }
    QString description() const override {
        return "Convert raster data type (byte, int16, int32, float32, float64).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::combo("data_type", "Output data type",
                {"Byte", "Int16", "Int32", "Float32", "Float64"}, 3,
                "Target data type for the output raster"),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        int typeIdx = parameter("data_type").toInt();
        DataType targetType;
        switch (typeIdx) {
            case 0: targetType = DataType::Byte;    break;
            case 1: targetType = DataType::Int16;   break;
            case 2: targetType = DataType::Int32;   break;
            case 3: targetType = DataType::Float32; break;
            case 4: targetType = DataType::Float64; break;
            default: targetType = DataType::Float64; break;
        }

        int cols = raster->cols(), rows = raster->rows(), bands = raster->bands();
        Raster output(cols, rows, bands, targetType);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        if (raster->hasNoData())
            output.setNoDataValue(raster->noDataValue());

        for (int b = 0; b < bands; ++b) {
            const auto& src = raster->data(b);
            auto& dst = output.data(b);
            int64_t total = static_cast<int64_t>(cols) * rows;

            for (int64_t i = 0; i < total; ++i) {
                double val = src[i];
                // Clamp for integer types
                switch (typeIdx) {
                    case 0: val = std::round(std::max(0.0, std::min(255.0, val))); break;
                    case 1: val = std::round(std::max(-32768.0, std::min(32767.0, val))); break;
                    case 2: val = std::round(std::max(-2147483648.0, std::min(2147483647.0, val))); break;
                    default: break;
                }
                dst[i] = val;
            }

            reportProgress(static_cast<double>(b + 1) / bands);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ConvertModule)

} // namespace aplaceholder
