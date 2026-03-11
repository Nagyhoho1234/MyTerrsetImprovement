#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class HardenModule : public Module {
public:
    QString name() const override { return "HARDEN"; }
    QString description() const override {
        return "Converts soft classification (probability/membership surfaces) to hard "
               "classification by assigning each pixel to the class with highest probability.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_rasters", "Input probability rasters (comma-separated)",
                "Comma-separated list of probability/membership surface file paths, "
                "one per class. Class ID corresponds to the order (1-based)."),
            ParameterDef::output("output", "Output hard classified image"),
            ParameterDef::real("threshold", "Minimum probability threshold", 0.0, 0.0, 1.0,
                "Minimum probability value for assignment. Pixels below this in all "
                "classes are left unclassified (0). Set to 0 to disable."),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Read input probability rasters
        // ------------------------------------------------------------------
        QString rastersParam = parameter("input_rasters").toString();
        QStringList rasterPaths = rastersParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : rasterPaths) p = p.trimmed();

        if (rasterPaths.size() < 2) {
            setError("At least 2 probability rasters are required");
            return false;
        }

        int numClasses = rasterPaths.size();
        std::vector<std::unique_ptr<Raster>> probRasters;
        for (const auto& path : rasterPaths) {
            auto r = GdalIO::read(path);
            if (!r) {
                setError("Failed to read probability raster: " + path);
                return false;
            }
            probRasters.push_back(std::move(r));
        }

        int cols = probRasters[0]->cols();
        int rows = probRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = probRasters[0]->hasNoData();
        double noData = probRasters[0]->noDataValue();

        for (int c = 1; c < numClasses; ++c) {
            if (probRasters[c]->cols() != cols || probRasters[c]->rows() != rows) {
                setError("All input probability rasters must have the same dimensions");
                return false;
            }
        }

        double threshold = parameter("threshold").toDouble();
        bool useThreshold = (threshold > 0.0);

        reportProgress(0.1, "Hardening classification...");

        // ------------------------------------------------------------------
        // 2. Collect data pointers
        // ------------------------------------------------------------------
        std::vector<const std::vector<double>*> probData(numClasses);
        for (int c = 0; c < numClasses; ++c)
            probData[c] = &probRasters[c]->data(0);

        // ------------------------------------------------------------------
        // 3. Create output classified image
        // ------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(probRasters[0]->geoTransform());
        output.setProjection(probRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);

        // ------------------------------------------------------------------
        // 4. For each pixel, find class with maximum probability
        // ------------------------------------------------------------------
        for (int64_t i = 0; i < total; ++i) {
            // Check NoData (use first raster as reference)
            if (hasND && (*probData[0])[i] == noData) {
                out[i] = 0;
                continue;
            }

            double bestProb = -std::numeric_limits<double>::max();
            int bestClass = 0; // 0 = unclassified

            for (int c = 0; c < numClasses; ++c) {
                double prob = (*probData[c])[i];
                if (prob > bestProb) {
                    bestProb = prob;
                    bestClass = c + 1; // 1-based class ID
                }
            }

            // Apply threshold
            if (useThreshold && bestProb < threshold) {
                bestClass = 0; // unclassified
            }

            out[i] = bestClass;

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        // ------------------------------------------------------------------
        // 5. Write output
        // ------------------------------------------------------------------
        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(HardenModule)

} // namespace aplaceholder
