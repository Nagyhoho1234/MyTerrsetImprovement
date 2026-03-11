#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <unordered_map>

namespace aplaceholder {

class SegClassModule : public Module {
public:
    QString name() const override { return "SEGCLASS"; }
    QString description() const override {
        return "Majority rule segment classifier. Assigns each segment the class that "
               "dominates within it, based on a prior pixel-based classification.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("segment_raster", "Segment image (region IDs)",
                "Segmented raster from SEGMENT or SEGMENTATION"),
            ParameterDef::file("training_raster", "Training/classification raster",
                "Pixel-based classification or training site raster"),
            ParameterDef::output("output", "Output classified image"),
        };
    }

    bool execute() override {
        auto segRaster = GdalIO::read(parameter("segment_raster").toString());
        if (!segRaster) {
            setError("Failed to read segment raster");
            return false;
        }

        auto trainRaster = GdalIO::read(parameter("training_raster").toString());
        if (!trainRaster) {
            setError("Failed to read training raster");
            return false;
        }

        int cols = segRaster->cols();
        int rows = segRaster->rows();

        if (trainRaster->cols() != cols || trainRaster->rows() != rows) {
            setError("Segment and training rasters must have the same dimensions");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        double segND = segRaster->noDataValue();
        bool segHasND = segRaster->hasNoData();
        double trainND = trainRaster->noDataValue();
        bool trainHasND = trainRaster->hasNoData();

        const auto& segData = segRaster->data(0);
        const auto& trainData = trainRaster->data(0);

        reportProgress(0.1, "Counting class frequencies per segment...");

        // For each segment ID, count occurrences of each class
        // segClassCounts[segID][classID] = count
        std::unordered_map<int, std::unordered_map<int, int64_t>> segClassCounts;

        for (int64_t i = 0; i < total; ++i) {
            double sv = segData[i];
            double tv = trainData[i];

            if (segHasND && sv == segND) continue;
            if (trainHasND && tv == trainND) continue;

            int segID = static_cast<int>(sv);
            int classID = static_cast<int>(tv);

            segClassCounts[segID][classID]++;

            if (i % (total / 20 + 1) == 0)
                reportProgress(0.1 + 0.5 * static_cast<double>(i) / total);
        }

        reportProgress(0.6, "Determining majority class per segment...");

        // Determine majority class for each segment
        std::unordered_map<int, int> segMajorityClass;
        for (auto& [segID, classCounts] : segClassCounts) {
            int bestClass = 0;
            int64_t bestCount = -1;
            for (auto& [classID, count] : classCounts) {
                if (count > bestCount) {
                    bestCount = count;
                    bestClass = classID;
                }
            }
            segMajorityClass[segID] = bestClass;
        }

        reportProgress(0.8, "Writing output...");

        // Create output raster
        double outND = -9999.0;
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(segRaster->geoTransform());
        output.setProjection(segRaster->projection());
        output.setNoDataValue(outND);

        auto& outData = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            double sv = segData[i];
            if (segHasND && sv == segND) {
                outData[i] = outND;
                continue;
            }
            int segID = static_cast<int>(sv);
            auto it = segMajorityClass.find(segID);
            if (it != segMajorityClass.end()) {
                outData[i] = static_cast<double>(it->second);
            } else {
                outData[i] = outND;
            }
        }

        reportProgress(1.0, "Complete.");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(SegClassModule)

} // namespace aplaceholder
