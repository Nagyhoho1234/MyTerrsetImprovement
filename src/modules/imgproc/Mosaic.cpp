#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <QStringList>
#include <algorithm>
#include <limits>

namespace aplaceholder {

class MosaicModule : public Module {
public:
    QString name() const override { return "MOSAIC"; }
    QString description() const override {
        return "Mosaic multiple rasters into a single output image. Computes the combined "
               "output extent from all inputs, handles NoData regions, and supports "
               "multiple blending methods for overlap zones (first, last, average, max).";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_rasters", "Input rasters (comma-separated file paths)"),
            ParameterDef::output("output", "Output mosaic image"),
            ParameterDef::combo("method", "Overlap blending method",
                                {"first", "last", "average", "max"}, 0,
                                "How to resolve pixel values where rasters overlap."),
        };
    }

    bool execute() override {
        QString inputStr = parameter("input_rasters").toString();
        QString outputPath = parameter("output").toString();
        int methodIdx = parameter("method").toInt();

        QStringList rasterPaths = inputStr.split(',', Qt::SkipEmptyParts);
        for (auto& p : rasterPaths)
            p = p.trimmed();

        if (rasterPaths.size() < 2) {
            setError("At least two input rasters are required for mosaicing");
            return false;
        }

        // Read all input rasters
        std::vector<std::unique_ptr<Raster>> rasters;
        for (const auto& path : rasterPaths) {
            auto raster = GdalIO::read(path);
            if (!raster) {
                setError("Failed to read input raster: " + path);
                return false;
            }
            rasters.push_back(std::move(raster));
        }

        reportProgress(0.1, "Computing output extent...");

        // Compute the combined bounding box in geographic coordinates
        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();

        // Use pixel size from the first raster as reference
        double pixelW = rasters[0]->geoTransform().pixelWidth;
        double pixelH = rasters[0]->geoTransform().pixelHeight;

        for (const auto& r : rasters) {
            const auto& gt = r->geoTransform();
            double x0 = gt.originX;
            double y0 = gt.originY;
            double x1 = x0 + r->cols() * gt.pixelWidth;
            double y1 = y0 + r->rows() * gt.pixelHeight;

            minX = std::min(minX, std::min(x0, x1));
            maxX = std::max(maxX, std::max(x0, x1));
            minY = std::min(minY, std::min(y0, y1));
            maxY = std::max(maxY, std::max(y0, y1));
        }

        // Compute output dimensions
        int outCols = static_cast<int>(std::ceil((maxX - minX) / std::abs(pixelW)));
        int outRows = static_cast<int>(std::ceil((maxY - minY) / std::abs(pixelH)));

        if (outCols <= 0 || outRows <= 0) {
            setError("Computed output dimensions are invalid");
            return false;
        }

        double noData = -9999.0;

        // Create output raster
        Raster output(outCols, outRows, 1, DataType::Float64);
        GeoTransform outGT;
        outGT.originX = minX;
        // If pixelHeight is negative (common for north-up images), origin Y is maxY
        outGT.originY = (pixelH < 0) ? maxY : minY;
        outGT.pixelWidth = pixelW;
        outGT.pixelHeight = pixelH;
        output.setGeoTransform(outGT);
        output.setProjection(rasters[0]->projection());
        output.setNoDataValue(noData);

        auto& outData = output.data(0);
        int64_t outTotal = static_cast<int64_t>(outCols) * outRows;

        // Initialize output to NoData
        for (int64_t i = 0; i < outTotal; ++i)
            outData[i] = noData;

        // For average method, track count of contributions
        std::vector<int> countMap;
        if (methodIdx == 2) { // average
            countMap.resize(outTotal, 0);
        }

        reportProgress(0.2, "Mosaicing rasters...");

        // Process each input raster
        int numRasters = static_cast<int>(rasters.size());
        for (int r = 0; r < numRasters; ++r) {
            const auto& raster = rasters[r];
            const auto& gt = raster->geoTransform();
            const auto& inData = raster->data(0);
            bool hasND = raster->hasNoData();
            double rasterND = raster->noDataValue();
            int rCols = raster->cols();
            int rRows = raster->rows();

            for (int row = 0; row < rRows; ++row) {
                for (int col = 0; col < rCols; ++col) {
                    double val = inData[static_cast<int64_t>(row) * rCols + col];
                    if (hasND && val == rasterND)
                        continue;

                    // Convert input pixel to geographic coordinates
                    double geoX = gt.originX + col * gt.pixelWidth;
                    double geoY = gt.originY + row * gt.pixelHeight;

                    // Convert to output pixel coordinates
                    int outCol = static_cast<int>((geoX - outGT.originX) / outGT.pixelWidth);
                    int outRow = static_cast<int>((geoY - outGT.originY) / outGT.pixelHeight);

                    if (outCol < 0 || outCol >= outCols || outRow < 0 || outRow >= outRows)
                        continue;

                    int64_t outIdx = static_cast<int64_t>(outRow) * outCols + outCol;

                    switch (methodIdx) {
                    case 0: // first: keep first valid value
                        if (outData[outIdx] == noData)
                            outData[outIdx] = val;
                        break;

                    case 1: // last: always overwrite
                        outData[outIdx] = val;
                        break;

                    case 2: // average: accumulate sum and count
                        if (outData[outIdx] == noData) {
                            outData[outIdx] = val;
                            countMap[outIdx] = 1;
                        } else {
                            outData[outIdx] += val;
                            countMap[outIdx]++;
                        }
                        break;

                    case 3: // max: keep maximum value
                        if (outData[outIdx] == noData || val > outData[outIdx])
                            outData[outIdx] = val;
                        break;
                    }
                }
            }

            reportProgress(0.2 + 0.7 * static_cast<double>(r + 1) / numRasters);
        }

        // Finalize average method: divide accumulated sums by counts
        if (methodIdx == 2) {
            for (int64_t i = 0; i < outTotal; ++i) {
                if (countMap[i] > 1)
                    outData[i] /= countMap[i];
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
    }
};

REGISTER_MODULE(MosaicModule)

} // namespace aplaceholder
