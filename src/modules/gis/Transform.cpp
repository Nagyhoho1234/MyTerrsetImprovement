#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class TransformModule : public Module {
public:
    QString name() const override { return "TRANSFORM"; }
    QString description() const override {
        return "Affine geometric transformation: rotation, scaling, and translation.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::real("rotation_degrees", "Rotation (degrees)", 0.0, -360.0, 360.0),
            ParameterDef::real("scale_x", "Scale X", 1.0, 0.001, 1000.0),
            ParameterDef::real("scale_y", "Scale Y", 1.0, 0.001, 1000.0),
            ParameterDef::real("translate_x", "Translate X", 0.0, -1e15, 1e15),
            ParameterDef::real("translate_y", "Translate Y", 0.0, -1e15, 1e15),
            ParameterDef::combo("method", "Resampling method",
                {"Nearest", "Bilinear"}, 0),
        };
    }

    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        double rotDeg = parameter("rotation_degrees").toDouble();
        double scaleX = parameter("scale_x").toDouble();
        double scaleY = parameter("scale_y").toDouble();
        double transX = parameter("translate_x").toDouble();
        double transY = parameter("translate_y").toDouble();
        int method = parameter("method").toInt();

        double rotRad = rotDeg * M_PI / 180.0;
        double cosR = std::cos(rotRad);
        double sinR = std::sin(rotRad);

        int srcCols = raster->cols(), srcRows = raster->rows();
        GeoTransform srcGt = raster->geoTransform();

        // Compute output dimensions to encompass transformed raster
        double corners[4][2] = {
            {0, 0},
            {static_cast<double>(srcCols), 0},
            {static_cast<double>(srcCols), static_cast<double>(srcRows)},
            {0, static_cast<double>(srcRows)}
        };

        double minX = 1e30, minY = 1e30, maxX = -1e30, maxY = -1e30;
        for (auto& c : corners) {
            double tx = (c[0] * cosR * scaleX - c[1] * sinR * scaleY) + transX / srcGt.pixelWidth;
            double ty = (c[0] * sinR * scaleX + c[1] * cosR * scaleY) + transY / srcGt.pixelHeight;
            minX = std::min(minX, tx);
            minY = std::min(minY, ty);
            maxX = std::max(maxX, tx);
            maxY = std::max(maxY, ty);
        }

        int outCols = static_cast<int>(std::ceil(maxX - minX));
        int outRows = static_cast<int>(std::ceil(maxY - minY));
        if (outCols <= 0 || outRows <= 0) {
            setError("Transformed raster has zero extent");
            return false;
        }

        Raster output(outCols, outRows, raster->bands(), DataType::Float64);
        double noData = raster->hasNoData() ? raster->noDataValue() : -9999.0;
        output.setNoDataValue(noData);
        output.setProjection(raster->projection());

        GeoTransform outGt = srcGt;
        outGt.originX += minX * srcGt.pixelWidth;
        outGt.originY += minY * srcGt.pixelHeight;
        output.setGeoTransform(outGt);

        // Inverse transform: from output pixel to source pixel
        double invScaleX = 1.0 / scaleX;
        double invScaleY = 1.0 / scaleY;

        for (int b = 0; b < raster->bands(); ++b) {
            const auto& srcData = raster->data(b);
            auto& dstData = output.data(b);

            for (int row = 0; row < outRows; ++row) {
                for (int col = 0; col < outCols; ++col) {
                    // Output pixel in transformed space
                    double px = col + minX - transX / srcGt.pixelWidth;
                    double py = row + minY - transY / srcGt.pixelHeight;

                    // Inverse rotation and scale
                    double srcX = (px * cosR + py * sinR) * invScaleX;
                    double srcY = (-px * sinR + py * cosR) * invScaleY;

                    int64_t outIdx = static_cast<int64_t>(row) * outCols + col;

                    if (method == 0) {
                        // Nearest neighbor
                        int sc = static_cast<int>(std::round(srcX));
                        int sr = static_cast<int>(std::round(srcY));
                        if (sc >= 0 && sc < srcCols && sr >= 0 && sr < srcRows) {
                            dstData[outIdx] = srcData[static_cast<int64_t>(sr) * srcCols + sc];
                        } else {
                            dstData[outIdx] = noData;
                        }
                    } else {
                        // Bilinear interpolation
                        int x0 = static_cast<int>(std::floor(srcX));
                        int y0 = static_cast<int>(std::floor(srcY));
                        int x1 = x0 + 1, y1 = y0 + 1;
                        double fx = srcX - x0, fy = srcY - y0;

                        if (x0 >= 0 && x1 < srcCols && y0 >= 0 && y1 < srcRows) {
                            double v00 = srcData[static_cast<int64_t>(y0) * srcCols + x0];
                            double v10 = srcData[static_cast<int64_t>(y0) * srcCols + x1];
                            double v01 = srcData[static_cast<int64_t>(y1) * srcCols + x0];
                            double v11 = srcData[static_cast<int64_t>(y1) * srcCols + x1];

                            bool hasNd = raster->hasNoData();
                            if (hasNd && (v00 == noData || v10 == noData ||
                                          v01 == noData || v11 == noData)) {
                                dstData[outIdx] = noData;
                            } else {
                                dstData[outIdx] = v00 * (1 - fx) * (1 - fy) +
                                                  v10 * fx * (1 - fy) +
                                                  v01 * (1 - fx) * fy +
                                                  v11 * fx * fy;
                            }
                        } else {
                            dstData[outIdx] = noData;
                        }
                    }
                }

                if (row % 100 == 0)
                    reportProgress(static_cast<double>(b * outRows + row) /
                                   (raster->bands() * outRows));
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(TransformModule)

} // namespace aplaceholder
