#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class MaxSetModule : public Module {
public:
    QString name() const override { return "MAXSET"; }
    QString description() const override {
        return "Select Maximum from Multiple Class Maps. Takes multiple soft "
               "classification outputs (e.g., from BAYCLASS, FUZCLASS, or UNMIX) "
               "and for each pixel selects the input with the highest value. "
               "Outputs a single raster where each pixel value is the index "
               "(1-based) of the winning input, plus a raster of the maximum values.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("inputs", "Input images (comma-separated paths)",
                "Comma-separated list of soft classification raster paths"),
            ParameterDef::output("output", "Output maximum-class image"),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Parse input file paths
        // ------------------------------------------------------------------
        QString inputsParam = parameter("inputs").toString();
        QStringList inputPaths = inputsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : inputPaths)
            p = p.trimmed();

        if (inputPaths.size() < 2) {
            setError("At least two input images are required");
            return false;
        }

        int numInputs = inputPaths.size();

        // ------------------------------------------------------------------
        // 2. Read all input rasters
        // ------------------------------------------------------------------
        reportProgress(0.0, "Reading input images...");

        std::vector<std::unique_ptr<Raster>> rasters(numInputs);
        int cols = 0, rows = 0;

        for (int r = 0; r < numInputs; ++r) {
            rasters[r] = GdalIO::read(inputPaths[r]);
            if (!rasters[r]) {
                setError("Failed to read input image: " + inputPaths[r]);
                return false;
            }
            if (r == 0) {
                cols = rasters[r]->cols();
                rows = rasters[r]->rows();
            } else if (rasters[r]->cols() != cols || rasters[r]->rows() != rows) {
                setError("Input image dimensions do not match: " + inputPaths[r]);
                return false;
            }
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = rasters[0]->hasNoData();
        double noData = rasters[0]->noDataValue();

        std::vector<const std::vector<double>*> data(numInputs);
        for (int r = 0; r < numInputs; ++r)
            data[r] = &rasters[r]->data(0);

        // ------------------------------------------------------------------
        // 3. For each pixel, find the input with maximum value
        // ------------------------------------------------------------------
        reportProgress(0.1, "Computing maximum class...");

        // Output 1: class index image (1-based index of winning input)
        Raster classOutput(cols, rows, 1, DataType::Int32);
        classOutput.setGeoTransform(rasters[0]->geoTransform());
        classOutput.setProjection(rasters[0]->projection());
        classOutput.setNoDataValue(0);

        // Output 2: maximum value image
        Raster maxValOutput(cols, rows, 1, DataType::Float32);
        maxValOutput.setGeoTransform(rasters[0]->geoTransform());
        maxValOutput.setProjection(rasters[0]->projection());
        if (hasND) maxValOutput.setNoDataValue(noData);

        auto& classOut = classOutput.data(0);
        auto& maxValOut = maxValOutput.data(0);

        for (int64_t i = 0; i < total; ++i) {
            // Check for NoData in the first input as indicator
            if (hasND && (*data[0])[i] == noData) {
                classOut[i] = 0;
                maxValOut[i] = noData;
                continue;
            }

            double maxVal = -std::numeric_limits<double>::max();
            int maxIdx = 0;

            for (int r = 0; r < numInputs; ++r) {
                double val = (*data[r])[i];
                if (val > maxVal) {
                    maxVal = val;
                    maxIdx = r + 1; // 1-based
                }
            }

            classOut[i] = maxIdx;
            maxValOut[i] = maxVal;

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        // ------------------------------------------------------------------
        // 4. Write outputs
        // ------------------------------------------------------------------
        reportProgress(0.95, "Writing output images...");

        QString outPath = parameter("output").toString();
        if (!GdalIO::write(classOutput, outPath)) {
            setError("Failed to write class output image: " + outPath);
            return false;
        }

        // Write max-value companion image
        QString maxValPath = outPath;
        int dotPos = maxValPath.lastIndexOf('.');
        if (dotPos > 0)
            maxValPath = maxValPath.left(dotPos) + "_maxval" + maxValPath.mid(dotPos);
        else
            maxValPath += "_maxval.tif";

        if (!GdalIO::write(maxValOutput, maxValPath)) {
            setError("Failed to write max-value image: " + maxValPath);
            return false;
        }

        reportProgress(1.0, QString("MaxSet complete: %1 inputs evaluated per pixel.")
                       .arg(numInputs));
        return true;
    }
};

REGISTER_MODULE(MaxSetModule)

} // namespace aplaceholder
