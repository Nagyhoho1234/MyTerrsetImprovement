#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class NdviCompModule : public Module {
public:
    QString name() const override { return "NDVICOMP"; }
    QString description() const override {
        return "NDVI maximum value composite. Selects the highest NDVI value per pixel "
               "across multiple date images, producing a cloud-free composite "
               "that represents peak vegetation conditions.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("inputs", "Input NDVI rasters (comma-separated)",
                "Comma-separated list of NDVI raster file paths from multiple dates"),
            ParameterDef::output("output", "Output maximum value composite"),
        };
    }

    bool execute() override {
        QString inputsParam = parameter("inputs").toString();
        QStringList inputPaths = inputsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : inputPaths)
            p = p.trimmed();

        QString outputPath = parameter("output").toString();

        int numInputs = inputPaths.size();
        if (numInputs < 2) {
            setError("At least 2 NDVI images are required for compositing.");
            return false;
        }

        // Read all inputs
        reportProgress(0.0, "Reading input NDVI images...");
        std::vector<std::unique_ptr<Raster>> rasters(numInputs);
        int cols = 0, rows = 0;

        for (int n = 0; n < numInputs; ++n) {
            rasters[n] = GdalIO::read(inputPaths[n]);
            if (!rasters[n]) {
                setError("Failed to read input: " + inputPaths[n]);
                return false;
            }
            if (n == 0) {
                cols = rasters[n]->cols();
                rows = rasters[n]->rows();
            } else if (rasters[n]->cols() != cols || rasters[n]->rows() != rows) {
                setError("All input images must have the same dimensions. "
                         "Mismatch in: " + inputPaths[n]);
                return false;
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = rasters[0]->hasNoData();
        double noData = rasters[0]->noDataValue();
        double outNoData = -9999.0;

        Raster output(cols, rows, 1, DataType::Float32);
        output.setGeoTransform(rasters[0]->geoTransform());
        output.setProjection(rasters[0]->projection());
        output.setNoDataValue(outNoData);
        auto& outData = output.data(0);

        reportProgress(0.2, "Computing maximum value composite...");

        // Collect data pointers
        std::vector<const std::vector<double>*> data(numInputs);
        for (int n = 0; n < numInputs; ++n)
            data[n] = &rasters[n]->data(0);

        for (int64_t i = 0; i < total; ++i) {
            double maxVal = -1e30;
            bool anyValid = false;

            for (int n = 0; n < numInputs; ++n) {
                double val = (*data[n])[i];
                if (hasND && val == noData)
                    continue;
                if (!anyValid || val > maxVal) {
                    maxVal = val;
                    anyValid = true;
                }
            }

            outData[i] = anyValid ? maxVal : outNoData;

            if (i % 1000000 == 0)
                reportProgress(0.2 + 0.7 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
    }
};

REGISTER_MODULE(NdviCompModule)

} // namespace aplaceholder
