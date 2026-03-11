#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class LogicModule : public Module {
public:
    QString name() const override { return "LOGIC"; }
    QString description() const override {
        return "Combined logical operations on two rasters. "
               "Supports AND (min), OR (max), XOR (a+b-2*min(a,b)), "
               "and NOT (1-a, uses first input only).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First input raster image"),
            ParameterDef::file("input2", "Second input raster image (ignored for NOT)"),
            ParameterDef::combo("operation", "Logical operation", {"AND", "OR", "XOR", "NOT"}, 0),
            ParameterDef::output("output", "Output raster image"),
        };
    }

    bool execute() override {
        QString op = parameter("operation").toString().toUpper();

        auto r1 = GdalIO::read(parameter("input1").toString());
        if (!r1) {
            setError("Failed to read first input raster");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& d1 = r1->data(0);
        double noData1 = r1->noDataValue();
        bool hasND1 = r1->hasNoData();

        // For NOT, we only need the first raster
        std::unique_ptr<Raster> r2;
        double noData2 = 0.0;
        bool hasND2 = false;

        if (op != "NOT") {
            r2 = GdalIO::read(parameter("input2").toString());
            if (!r2) {
                setError("Failed to read second input raster");
                return false;
            }
            if (r2->cols() != cols || r2->rows() != rows) {
                setError("Input rasters must have the same dimensions");
                return false;
            }
            noData2 = r2->noDataValue();
            hasND2 = r2->hasNoData();
        }

        reportProgress(0.1, "Computing logical operation: " + op + "...");

        Raster result(cols, rows, 1, DataType::Float64);
        result.setGeoTransform(r1->geoTransform());
        result.setProjection(r1->projection());
        result.setNoDataValue(r1->noDataValue());
        auto& out = result.data(0);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND1 && d1[i] == noData1) {
                out[i] = noData1;
                continue;
            }

            if (op == "NOT") {
                out[i] = 1.0 - d1[i];
            } else {
                const auto& d2 = r2->data(0);
                if (hasND2 && d2[i] == noData2) {
                    out[i] = noData1;
                    continue;
                }

                double a = d1[i];
                double b = d2[i];

                if (op == "AND") {
                    out[i] = std::min(a, b);
                } else if (op == "OR") {
                    out[i] = std::max(a, b);
                } else if (op == "XOR") {
                    out[i] = a + b - 2.0 * std::min(a, b);
                } else {
                    setError("Unknown operation: " + op);
                    return false;
                }
            }

            if (i % 100000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        reportProgress(0.95, "Writing output raster...");

        if (!GdalIO::write(result, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0, "Logical operation " + op + " complete.");
        return true;
    }
};

REGISTER_MODULE(LogicModule)

} // namespace aplaceholder
