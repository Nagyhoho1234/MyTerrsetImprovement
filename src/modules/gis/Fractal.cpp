#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>

namespace aplaceholder {

class FractalModule : public Module {
public:
    QString name() const override { return "FRACTAL"; }
    QString description() const override {
        return "Compute local fractal dimension from a DEM using the triangular prism "
               "surface area method. The fractal dimension ranges from 2.0 (smooth) "
               "to 3.0 (extremely rough).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input DEM raster"),
            ParameterDef::output("output", "Output fractal dimension image"),
            ParameterDef::integer("max_window", "Maximum window size", 9, 3, 25,
                "Maximum window size for multi-scale analysis (must be odd)"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input DEM"); return false; }

        int cols = input->cols(), rows = input->rows();
        int maxWin = parameter("max_window").toInt();
        if (maxWin % 2 == 0) maxWin++;

        const auto& data = input->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        double cs = std::abs(input->geoTransform().pixelWidth);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(input->geoTransform());
        output.setProjection(input->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing fractal dimension...");

        // Multi-scale surface area approach
        // For each scale (window size), compute surface area in the neighborhood
        // Then regress log(area) vs log(scale) to get fractal dimension

        // Collect window sizes: 3, 5, 7, ... up to maxWin
        std::vector<int> winSizes;
        for (int w = 3; w <= maxWin; w += 2) winSizes.push_back(w);

        if (winSizes.size() < 2) {
            setError("Need at least 2 scales; increase maximum window size");
            return false;
        }

        int nScales = static_cast<int>(winSizes.size());

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t idx = static_cast<size_t>(r) * cols + c;
                if (hasND && data[idx] == noData) {
                    out[idx] = noData;
                    continue;
                }

                // For each scale, compute surface area using triangular prism method
                double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

                for (int s = 0; s < nScales; ++s) {
                    int half = winSizes[s] / 2;
                    double area = 0.0;
                    int triangles = 0;

                    for (int kr = -half; kr < half; ++kr) {
                        for (int kc = -half; kc < half; ++kc) {
                            int r0 = r + kr, c0 = c + kc;
                            int r1 = r + kr + 1, c1 = c + kc + 1;

                            auto safeVal = [&](int rr, int cc) -> double {
                                rr = std::max(0, std::min(rows - 1, rr));
                                cc = std::max(0, std::min(cols - 1, cc));
                                size_t i = static_cast<size_t>(rr) * cols + cc;
                                if (hasND && data[i] == noData) return data[idx];
                                return data[i];
                            };

                            double z00 = safeVal(r0, c0);
                            double z01 = safeVal(r0, c1);
                            double z10 = safeVal(r1, c0);
                            double z11 = safeVal(r1, c1);

                            // Two triangles per cell
                            // Triangle 1: (0,0), (1,0), (0,1)
                            double a1 = cs, b1 = cs, c1h = std::sqrt(cs*cs + (z10-z00)*(z10-z00));
                            double s1 = (a1 + b1 + c1h) / 2.0;
                            double t1 = s1 * (s1-a1) * (s1-b1) * (s1-c1h);
                            if (t1 > 0) area += std::sqrt(t1);

                            // Triangle 2: (1,1), (1,0), (0,1)
                            double a2 = cs, b2 = cs, c2 = std::sqrt(cs*cs + (z11-z01)*(z11-z01));
                            double s2 = (a2 + b2 + c2) / 2.0;
                            double t2 = s2 * (s2-a2) * (s2-b2) * (s2-c2);
                            if (t2 > 0) area += std::sqrt(t2);

                            triangles += 2;
                        }
                    }

                    if (area <= 0.0 || triangles == 0) continue;

                    double logScale = std::log(winSizes[s] * cs);
                    double logArea = std::log(area);
                    sumX += logScale;
                    sumY += logArea;
                    sumXY += logScale * logArea;
                    sumXX += logScale * logScale;
                }

                // Linear regression: logArea = slope * logScale + intercept
                // Fractal dimension D = 2.0 + slope (from area ~ scale^D)
                // More precisely D = 2 - slope since area grows with scale
                double denom = nScales * sumXX - sumX * sumX;
                if (std::abs(denom) < 1e-15) {
                    out[idx] = 2.0;
                } else {
                    double slope = (nScales * sumXY - sumX * sumY) / denom;
                    // Fractal dimension: area ~ scale^(D-2), so slope = D-2
                    double D = 2.0 + std::abs(slope - 2.0);
                    // Clamp to valid range [2, 3]
                    D = std::max(2.0, std::min(3.0, D));
                    out[idx] = D;
                }
            }
            if (r % 50 == 0) reportProgress(static_cast<double>(r) / rows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(FractalModule)

} // namespace aplaceholder
