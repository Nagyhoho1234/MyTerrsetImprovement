#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class ResampleModule : public Module {
public:
    QString name() const override { return "RESAMPLE"; }
    QString description() const override {
        return "Resample raster to a different cell size using nearest-neighbor, "
               "bilinear, or cubic interpolation.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input image"),
            ParameterDef::output("output", "Output image"),
            ParameterDef::real("cell_size", "Output cell size", 30.0, 0.0001, 999999,
                               "Target cell size in map units"),
            ParameterDef::combo("method", "Interpolation method",
                {"nearest", "bilinear", "cubic"}, 0,
                "Resampling interpolation method"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input").toString());
        if (!r1) {
            setError("Failed to read input raster");
            return false;
        }

        double cellSize = parameter("cell_size").toDouble();
        int method = parameter("method").toInt();

        const auto& gt = r1->geoTransform();
        int srcCols = r1->cols(), srcRows = r1->rows();
        double noData = r1->noDataValue();
        bool hasND = r1->hasNoData();

        // Compute extent of input raster
        double extentW = std::abs(gt.pixelWidth) * srcCols;
        double extentH = std::abs(gt.pixelHeight) * srcRows;

        // Compute output dimensions
        int outCols = std::max(1, static_cast<int>(std::round(extentW / cellSize)));
        int outRows = std::max(1, static_cast<int>(std::round(extentH / cellSize)));

        Raster output(outCols, outRows, 1, DataType::Float64);
        GeoTransform outGt = gt;
        outGt.pixelWidth = (gt.pixelWidth > 0) ? cellSize : -cellSize;
        outGt.pixelHeight = (gt.pixelHeight > 0) ? cellSize : -cellSize;
        output.setGeoTransform(outGt);
        output.setProjection(r1->projection());
        if (hasND) output.setNoDataValue(noData);

        auto& out = output.data(0);

        for (int row = 0; row < outRows; ++row) {
            for (int col = 0; col < outCols; ++col) {
                // Map output pixel center to source coordinates
                double outX, outY;
                output.colRowToXY(col, row, outX, outY);

                // Convert to source fractional col/row
                double srcColF = (outX - gt.originX) / gt.pixelWidth - 0.5;
                double srcRowF = (outY - gt.originY) / gt.pixelHeight - 0.5;

                double val;
                if (method == 0) {
                    val = sampleNearest(*r1, srcColF, srcRowF);
                } else if (method == 1) {
                    val = sampleBilinear(*r1, srcColF, srcRowF);
                } else {
                    val = sampleCubic(*r1, srcColF, srcRowF);
                }

                output.setValue(col, row, val);
            }

            if (row % 100 == 0)
                reportProgress(static_cast<double>(row) / outRows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }

private:
    double safeValue(const Raster& r, int col, int row) const {
        if (col < 0 || col >= r.cols() || row < 0 || row >= r.rows())
            return r.hasNoData() ? r.noDataValue() : 0.0;
        return r.value(col, row);
    }

    double sampleNearest(const Raster& r, double col, double row) const {
        int c = static_cast<int>(std::round(col));
        int rr = static_cast<int>(std::round(row));
        return safeValue(r, c, rr);
    }

    double sampleBilinear(const Raster& r, double col, double row) const {
        int c0 = static_cast<int>(std::floor(col));
        int r0 = static_cast<int>(std::floor(row));
        double dc = col - c0;
        double dr = row - r0;

        double v00 = safeValue(r, c0, r0);
        double v10 = safeValue(r, c0 + 1, r0);
        double v01 = safeValue(r, c0, r0 + 1);
        double v11 = safeValue(r, c0 + 1, r0 + 1);

        double noData = r.noDataValue();
        bool hasND = r.hasNoData();
        if (hasND && (v00 == noData || v10 == noData || v01 == noData || v11 == noData))
            return noData;

        return v00 * (1 - dc) * (1 - dr) +
               v10 * dc * (1 - dr) +
               v01 * (1 - dc) * dr +
               v11 * dc * dr;
    }

    double sampleCubic(const Raster& r, double col, double row) const {
        int c0 = static_cast<int>(std::floor(col));
        int r0 = static_cast<int>(std::floor(row));
        double dc = col - c0;
        double dr = row - r0;

        double noData = r.noDataValue();
        bool hasND = r.hasNoData();

        // Get 4x4 neighborhood
        double rows4[4];
        for (int j = -1; j <= 2; ++j) {
            double cols4[4];
            for (int i = -1; i <= 2; ++i) {
                cols4[i + 1] = safeValue(r, c0 + i, r0 + j);
                if (hasND && cols4[i + 1] == noData)
                    return noData;
            }
            rows4[j + 1] = cubicInterp(cols4, dc);
        }
        return cubicInterp(rows4, dr);
    }

    static double cubicInterp(const double v[4], double t) {
        double a = -0.5 * v[0] + 1.5 * v[1] - 1.5 * v[2] + 0.5 * v[3];
        double b = v[0] - 2.5 * v[1] + 2.0 * v[2] - 0.5 * v[3];
        double c = -0.5 * v[0] + 0.5 * v[2];
        double d = v[1];
        return a * t * t * t + b * t * t + c * t + d;
    }
};

REGISTER_MODULE(ResampleModule)

} // namespace aplaceholder
