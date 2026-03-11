#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <memory>

namespace aplaceholder {

class OverlapUseModule : public Module {
public:
    QString name() const override { return "OVERLAP_USE"; }
    QString description() const override {
        return "Overlapping land use analysis. Counts how many input binary rasters "
               "have value=1 at each pixel, producing an overlap count raster. Useful "
               "for identifying areas of competing or complementary land uses, species "
               "richness from individual presence maps, or conservation priority zones.";
    }
    QString category() const override { return "Habitat & Biodiversity"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("inputs", "Input binary rasters (comma-separated)",
                "Comma-separated paths to binary (0/1) rasters to overlay"),
            ParameterDef::output("output", "Output overlap count raster",
                "Integer raster where each pixel = count of inputs with value 1"),
        };
    }

    bool execute() override {
        // 1. Parse comma-separated input paths
        QString inputsStr = parameter("inputs").toString();
        QStringList inputPaths = inputsStr.split(',', Qt::SkipEmptyParts);
        for (int i = 0; i < inputPaths.size(); ++i)
            inputPaths[i] = inputPaths[i].trimmed();

        int numInputs = inputPaths.size();
        if (numInputs < 1) {
            setError("At least one input raster is required.");
            return false;
        }

        // 2. Load all input rasters
        std::vector<std::unique_ptr<Raster>> rasters;
        for (int i = 0; i < numInputs; ++i) {
            auto r = GdalIO::read(inputPaths[i]);
            if (!r) {
                setError("Failed to read input raster: " + inputPaths[i]);
                return false;
            }
            rasters.push_back(std::move(r));
        }

        reportProgress(0.10, QString("%1 input rasters loaded.").arg(numInputs));

        // 3. Validate dimensions (all must match first raster)
        int cols = rasters[0]->cols();
        int rows = rasters[0]->rows();
        int64_t total = rasters[0]->cellCount();

        for (int i = 1; i < numInputs; ++i) {
            if (rasters[i]->cols() != cols || rasters[i]->rows() != rows) {
                setError(QString("Dimension mismatch: input 0 is %1x%2 but input %3 is %4x%5")
                             .arg(cols).arg(rows).arg(i)
                             .arg(rasters[i]->cols()).arg(rasters[i]->rows()));
                return false;
            }
        }

        // 4. Count overlaps at each pixel
        double noData = rasters[0]->noDataValue();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(rasters[0]->geoTransform());
        output.setProjection(rasters[0]->projection());
        output.setNoDataValue(noData);

        auto& outData = output.data(0);

        // Collect data pointers and nodata info
        std::vector<const std::vector<double>*> dataPtrs;
        std::vector<bool> hasNDVec;
        std::vector<double> ndVec;
        for (int i = 0; i < numInputs; ++i) {
            dataPtrs.push_back(&rasters[i]->data(0));
            hasNDVec.push_back(rasters[i]->hasNoData());
            ndVec.push_back(rasters[i]->noDataValue());
        }

        int64_t validCount = 0;

        for (int64_t i = 0; i < total; ++i) {
            // If any input is nodata at this pixel, output nodata
            bool anyNoData = false;
            for (int r = 0; r < numInputs; ++r) {
                if (hasNDVec[r] && (*dataPtrs[r])[i] == ndVec[r]) {
                    anyNoData = true;
                    break;
                }
            }

            if (anyNoData) {
                outData[i] = noData;
                continue;
            }

            int count = 0;
            for (int r = 0; r < numInputs; ++r) {
                if (std::round((*dataPtrs[r])[i]) != 0)
                    ++count;
            }

            outData[i] = static_cast<double>(count);
            ++validCount;

            if (i % 500000 == 0)
                reportProgress(0.10 + 0.80 * static_cast<double>(i) / total);
        }

        reportProgress(0.90, "Writing output...");

        // 5. Write output
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        auto stats = output.computeStats(0);
        reportProgress(1.0,
            QString("Done. Overlap count — min: %1, max: %2, mean: %3, valid pixels: %4")
                .arg(stats.min, 0, 'f', 0)
                .arg(stats.max, 0, 'f', 0)
                .arg(stats.mean, 0, 'f', 2)
                .arg(validCount));

        return true;
    }
};

REGISTER_MODULE(OverlapUseModule)

} // namespace aplaceholder
