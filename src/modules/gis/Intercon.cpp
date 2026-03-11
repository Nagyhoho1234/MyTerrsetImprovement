#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class InterconModule : public Module {
public:
    QString name() const override { return "INTERCON"; }
    QString description() const override {
        return "Interpolates a raster surface from rasterized isoline (contour) data "
               "using linear interpolation along four directional lines. For each unknown "
               "pixel, extends lines in 8 directions to find known values, then interpolates "
               "using the line with the greatest slope.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input contour raster",
                "Raster with contour cells holding elevation values, background = NoData"),
            ParameterDef::output("output", "Output interpolated surface"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) {
            setError("Failed to read input contour raster");
            return false;
        }

        int cols = input->cols();
        int rows = input->rows();
        double noData = input->hasNoData() ? input->noDataValue() : -9999.0;

        const auto& srcData = input->data(0);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);

        double dx = std::abs(input->geoTransform().pixelWidth);
        double dy = std::abs(input->geoTransform().pixelHeight);
        double dd = std::sqrt(dx * dx + dy * dy);

        // 4 line directions: horizontal, vertical, diagonal-1, diagonal-2
        // Each has two opposing search directions
        // dc, dr pairs for the 8 directions (4 lines x 2 ends)
        static const int dirs[4][2][2] = {
            {{ 1, 0}, {-1,  0}},   // horizontal
            {{ 0, 1}, { 0, -1}},   // vertical
            {{ 1, 1}, {-1, -1}},   // diagonal /
            {{ 1,-1}, {-1,  1}},   // diagonal backslash
        };

        auto isContour = [&](int c, int r) -> bool {
            if (c < 0 || c >= cols || r < 0 || r >= rows) return false;
            double v = srcData[static_cast<size_t>(r) * cols + c];
            return (v != noData && !std::isnan(v));
        };

        auto getValue = [&](int c, int r) -> double {
            return srcData[static_cast<size_t>(r) * cols + c];
        };

        reportProgress(0.0, "Interpolating from contours...");

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;

                // If this cell is already a contour, keep its value
                if (isContour(c, r)) {
                    outData[idx] = getValue(c, r);
                    continue;
                }

                double bestSlope = -1.0;
                double bestValue = noData;

                // Search along 4 lines
                for (int d = 0; d < 4; ++d) {
                    double stepDist = (d < 2) ? ((d == 0) ? dx : dy) : dd;

                    // Search in positive direction
                    double val1 = noData;
                    double dist1 = 0.0;
                    {
                        int sc = c, sr = r;
                        while (true) {
                            sc += dirs[d][0][0];
                            sr += dirs[d][0][1];
                            dist1 += stepDist;
                            if (sc < 0 || sc >= cols || sr < 0 || sr >= rows) {
                                dist1 = 0.0;
                                break;
                            }
                            if (isContour(sc, sr)) {
                                val1 = getValue(sc, sr);
                                break;
                            }
                        }
                    }

                    // Search in negative direction
                    double val2 = noData;
                    double dist2 = 0.0;
                    {
                        int sc = c, sr = r;
                        while (true) {
                            sc += dirs[d][1][0];
                            sr += dirs[d][1][1];
                            dist2 += stepDist;
                            if (sc < 0 || sc >= cols || sr < 0 || sr >= rows) {
                                dist2 = 0.0;
                                break;
                            }
                            if (isContour(sc, sr)) {
                                val2 = getValue(sc, sr);
                                break;
                            }
                        }
                    }

                    // Both endpoints must be found
                    if (dist1 <= 0.0 || dist2 <= 0.0) continue;
                    if (val1 == noData || val2 == noData) continue;

                    double totalDist = dist1 + dist2;
                    double slope = std::abs(val1 - val2) / totalDist;

                    // Use the line with greatest slope (per spec)
                    if (slope > bestSlope) {
                        bestSlope = slope;
                        // Linear interpolation weighted by distance
                        bestValue = val2 + (val1 - val2) * (dist2 / totalDist);
                    }
                }

                outData[idx] = bestValue;
            }

            if (r % 50 == 0)
                reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(InterconModule)

} // namespace aplaceholder
