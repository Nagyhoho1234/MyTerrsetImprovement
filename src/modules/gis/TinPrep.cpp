#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>

namespace aplaceholder {

class TinPrepModule : public Module {
public:
    QString name() const override { return "TINPREP"; }
    QString description() const override {
        return "Prepares data for TIN generation by thinning vertices along contour lines. "
               "Removes points where the angle between consecutive line segments is nearly "
               "straight (within tolerance), reducing point count while preserving shape.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_raster", "Input contour raster",
                "Raster with polyline/contour cells holding attribute values, background = NoData"),
            ParameterDef::output("output_points", "Output thinned point file (CSV)",
                "CSV file with x, y, value columns"),
            ParameterDef::real("angle_tolerance", "Angle tolerance (degrees)", 5.0, 0.1, 90.0,
                "Points where segment angle deviates less than this from straight are removed"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input_raster").toString());
        if (!input) {
            setError("Failed to read input raster");
            return false;
        }

        QString outPath = parameter("output_points").toString();
        double angleTol = parameter("angle_tolerance").toDouble();
        double cosThreshold = std::cos((180.0 - angleTol) * 3.14159265358979323846 / 180.0);

        int cols = input->cols();
        int rows = input->rows();
        double noData = input->hasNoData() ? input->noDataValue() : -9999.0;
        const auto& srcData = input->data(0);

        reportProgress(0.0, "Extracting contour points...");

        // Extract contour cell coordinates grouped by value
        struct PixelPt {
            int col, row;
            double value;
        };

        // Collect all contour pixels
        std::vector<PixelPt> allContourPixels;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                double val = srcData[static_cast<size_t>(r) * cols + c];
                if (val != noData && !std::isnan(val)) {
                    allContourPixels.push_back({c, r, val});
                }
            }
        }

        if (allContourPixels.empty()) {
            setError("No contour pixels found in input raster");
            return false;
        }

        reportProgress(0.3, "Tracing contour lines...");

        // Trace connected contour segments by following 8-connected neighbors
        // with same value. Build ordered chains.
        std::vector<std::vector<PixelPt>> chains;
        std::vector<std::vector<bool>> visited(rows, std::vector<bool>(cols, false));

        auto isContour = [&](int c, int r) -> bool {
            if (c < 0 || c >= cols || r < 0 || r >= rows) return false;
            double v = srcData[static_cast<size_t>(r) * cols + c];
            return (v != noData && !std::isnan(v));
        };

        // 8-connected neighbor offsets
        static const int dx8[] = {1, 1, 0, -1, -1, -1, 0, 1};
        static const int dy8[] = {0, 1, 1, 1, 0, -1, -1, -1};

        for (const auto& startPx : allContourPixels) {
            if (visited[startPx.row][startPx.col]) continue;

            // BFS/chain trace for connected contour pixels with same value
            std::vector<PixelPt> chain;
            std::queue<std::pair<int, int>> q;
            q.push({startPx.col, startPx.row});
            visited[startPx.row][startPx.col] = true;

            while (!q.empty()) {
                auto [cc, cr] = q.front();
                q.pop();
                double val = srcData[static_cast<size_t>(cr) * cols + cc];
                chain.push_back({cc, cr, val});

                for (int d = 0; d < 8; ++d) {
                    int nc = cc + dx8[d];
                    int nr = cr + dy8[d];
                    if (nc < 0 || nc >= cols || nr < 0 || nr >= rows) continue;
                    if (visited[nr][nc]) continue;
                    double nv = srcData[static_cast<size_t>(nr) * cols + nc];
                    if (nv == val) {
                        visited[nr][nc] = true;
                        q.push({nc, nr});
                    }
                }
            }

            if (!chain.empty())
                chains.push_back(std::move(chain));
        }

        reportProgress(0.6, "Thinning points by angle...");

        // For each chain, thin points where the angle is nearly straight
        struct OutPoint { double x, y, value; };
        std::vector<OutPoint> outputPoints;

        for (const auto& chain : chains) {
            if (chain.size() <= 2) {
                // Keep all points in very short chains
                for (const auto& px : chain) {
                    double x, y;
                    input->colRowToXY(px.col, px.row, x, y);
                    outputPoints.push_back({x, y, px.value});
                }
                continue;
            }

            // Always keep first point
            {
                double x, y;
                input->colRowToXY(chain[0].col, chain[0].row, x, y);
                outputPoints.push_back({x, y, chain[0].value});
            }

            // Check angle at each interior point
            for (size_t i = 1; i < chain.size() - 1; ++i) {
                double x0, y0, x1, y1, x2, y2;
                input->colRowToXY(chain[i-1].col, chain[i-1].row, x0, y0);
                input->colRowToXY(chain[i].col, chain[i].row, x1, y1);
                input->colRowToXY(chain[i+1].col, chain[i+1].row, x2, y2);

                double ax = x0 - x1, ay = y0 - y1;
                double bx = x2 - x1, by = y2 - y1;

                double lenA = std::sqrt(ax * ax + ay * ay);
                double lenB = std::sqrt(bx * bx + by * by);

                if (lenA < 1e-15 || lenB < 1e-15) continue; // skip coincident points

                double cosAngle = (ax * bx + ay * by) / (lenA * lenB);
                cosAngle = std::max(-1.0, std::min(1.0, cosAngle));

                // If angle is nearly straight (cos close to -1), skip this point
                if (cosAngle <= cosThreshold) {
                    continue; // nearly straight, remove
                }

                // Keep this point - significant direction change
                outputPoints.push_back({x1, y1, chain[i].value});
            }

            // Always keep last point
            {
                const auto& last = chain.back();
                double x, y;
                input->colRowToXY(last.col, last.row, x, y);
                outputPoints.push_back({x, y, last.value});
            }
        }

        reportProgress(0.9, "Writing output...");

        // Write CSV output
        {
            std::ofstream ofs(outPath.toStdString());
            if (!ofs.is_open()) {
                setError("Failed to write output file: " + outPath);
                return false;
            }
            ofs.precision(12);
            ofs << "x,y,value\n";
            for (const auto& pt : outputPoints) {
                ofs << pt.x << "," << pt.y << "," << pt.value << "\n";
            }
        }

        reportProgress(1.0, "Done");
        return true;
    }
};

REGISTER_MODULE(TinPrepModule)

} // namespace aplaceholder
