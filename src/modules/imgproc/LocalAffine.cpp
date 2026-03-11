#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <array>

namespace aplaceholder {

class LocalAffineModule : public Module {
public:
    QString name() const override { return "LOCALAFFINE"; }
    QString description() const override {
        return "Local affine transformation for geometric correction and image registration. "
               "Applies a piecewise affine transform using ground control points (GCPs) to "
               "warp the input raster to the correct geometry.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster image"),
            ParameterDef::file("gcp_file", "Ground control points file"),
            ParameterDef::output("output", "Output corrected image"),
        };
    }

private:
    struct GCP {
        double srcX, srcY;  // source (image) coordinates
        double dstX, dstY;  // destination (map) coordinates
    };

    struct Triangle {
        int i0, i1, i2;  // indices into GCP array
    };

    // Affine coefficients: dstX = a0 + a1*srcX + a2*srcY (and similarly for Y)
    struct AffineCoeffs {
        double a0, a1, a2;  // for X transform
        double b0, b1, b2;  // for Y transform
    };

    // Read GCPs from file: each line is "srcX srcY dstX dstY"
    bool readGCPs(const QString& path, std::vector<GCP>& gcps) {
        std::ifstream ifs(path.toStdString());
        if (!ifs.is_open()) return false;

        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            GCP g;
            if (ss >> g.srcX >> g.srcY >> g.dstX >> g.dstY) {
                gcps.push_back(g);
            }
        }
        return gcps.size() >= 3;
    }

    // Build triangulation using a simple fan approach from sorted GCPs
    void buildTriangles(const std::vector<GCP>& gcps, std::vector<Triangle>& tris) {
        int n = static_cast<int>(gcps.size());
        if (n < 3) return;

        // Compute centroid of destination points
        double cx = 0, cy = 0;
        for (const auto& g : gcps) { cx += g.dstX; cy += g.dstY; }
        cx /= n; cy /= n;

        // Sort GCP indices by angle from centroid
        std::vector<int> idx(n);
        for (int i = 0; i < n; ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            double angA = std::atan2(gcps[a].dstY - cy, gcps[a].dstX - cx);
            double angB = std::atan2(gcps[b].dstY - cy, gcps[b].dstX - cx);
            return angA < angB;
        });

        // Create fan triangulation
        for (int i = 0; i < n - 2; ++i) {
            tris.push_back({idx[0], idx[i + 1], idx[i + 2]});
        }
    }

    // Solve affine coefficients for a triangle (3 GCPs)
    bool solveAffine(const GCP& g0, const GCP& g1, const GCP& g2, AffineCoeffs& c) {
        double d0x = g0.dstX, d1x = g1.dstX, d2x = g2.dstX;
        double d0y = g0.dstY, d1y = g1.dstY, d2y = g2.dstY;
        double s0x = g0.srcX, s1x = g1.srcX, s2x = g2.srcX;
        double s0y = g0.srcY, s1y = g1.srcY, s2y = g2.srcY;

        // Determinant of the source coordinate matrix
        double det = (s1x - s0x) * (s2y - s0y) - (s2x - s0x) * (s1y - s0y);
        if (std::abs(det) < 1e-12) return false;

        double invDet = 1.0 / det;

        // X coefficients
        c.a1 = ((d1x - d0x) * (s2y - s0y) - (d2x - d0x) * (s1y - s0y)) * invDet;
        c.a2 = ((d2x - d0x) * (s1x - s0x) - (d1x - d0x) * (s2x - s0x)) * invDet;
        c.a0 = d0x - c.a1 * s0x - c.a2 * s0y;

        // Y coefficients
        c.b1 = ((d1y - d0y) * (s2y - s0y) - (d2y - d0y) * (s1y - s0y)) * invDet;
        c.b2 = ((d2y - d0y) * (s1x - s0x) - (d1y - d0y) * (s2x - s0x)) * invDet;
        c.b0 = d0y - c.b1 * s0x - c.b2 * s0y;

        return true;
    }

    // Check if point (px, py) is inside triangle defined by three destination GCPs
    bool pointInTriangle(double px, double py,
                         const GCP& g0, const GCP& g1, const GCP& g2) {
        double d1 = (px - g1.dstX) * (g0.dstY - g1.dstY) -
                     (g0.dstX - g1.dstX) * (py - g1.dstY);
        double d2 = (px - g2.dstX) * (g1.dstY - g2.dstY) -
                     (g1.dstX - g2.dstX) * (py - g2.dstY);
        double d3 = (px - g0.dstX) * (g2.dstY - g0.dstY) -
                     (g2.dstX - g0.dstX) * (py - g0.dstY);

        bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
        bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);
        return !(hasNeg && hasPos);
    }

public:
    bool execute() override {
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input raster");
            return false;
        }

        std::vector<GCP> gcps;
        if (!readGCPs(parameter("gcp_file").toString(), gcps)) {
            setError("Failed to read GCP file or fewer than 3 GCPs provided");
            return false;
        }

        int srcCols = raster->cols(), srcRows = raster->rows();
        const auto& srcGT = raster->geoTransform();

        // Build triangulation and compute affine coefficients per triangle
        std::vector<Triangle> triangles;
        buildTriangles(gcps, triangles);

        std::vector<AffineCoeffs> coeffs(triangles.size());
        for (size_t t = 0; t < triangles.size(); ++t) {
            const auto& tri = triangles[t];
            if (!solveAffine(gcps[tri.i0], gcps[tri.i1], gcps[tri.i2], coeffs[t])) {
                setError("Degenerate triangle in GCP triangulation");
                return false;
            }
        }

        // Use a global affine from first 3 GCPs as fallback
        AffineCoeffs globalCoeffs;
        solveAffine(gcps[0], gcps[1], gcps[2], globalCoeffs);

        // Output has same dimensions as input
        int dstCols = srcCols, dstRows = srcRows;
        Raster output(dstCols, dstRows, 1, DataType::Float64);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(raster->noDataValue());

        const auto& src = raster->data(0);
        auto& dst = output.data(0);
        double noData = raster->noDataValue();
        bool hasND = raster->hasNoData();

        reportProgress(0.0, "Applying local affine transformation...");

        for (int dstR = 0; dstR < dstRows; ++dstR) {
            for (int dstC = 0; dstC < dstCols; ++dstC) {
                int64_t dstIdx = static_cast<int64_t>(dstR) * dstCols + dstC;

                // Convert pixel coordinates to map coordinates
                double mapX = srcGT.originX + dstC * srcGT.pixelWidth + dstR * srcGT.rotationX;
                double mapY = srcGT.originY + dstC * srcGT.rotationY + dstR * srcGT.pixelHeight;

                // Find which triangle contains this point and apply its affine
                const AffineCoeffs* useCoeffs = &globalCoeffs;
                for (size_t t = 0; t < triangles.size(); ++t) {
                    const auto& tri = triangles[t];
                    if (pointInTriangle(mapX, mapY,
                                        gcps[tri.i0], gcps[tri.i1], gcps[tri.i2])) {
                        useCoeffs = &coeffs[t];
                        break;
                    }
                }

                // Inverse transform: map destination coords back to source pixel
                double srcMapX = useCoeffs->a0 + useCoeffs->a1 * mapX + useCoeffs->a2 * mapY;
                double srcMapY = useCoeffs->b0 + useCoeffs->b1 * mapX + useCoeffs->b2 * mapY;

                // Convert map coordinates back to pixel coordinates
                double srcPixX = (srcMapX - srcGT.originX) / srcGT.pixelWidth;
                double srcPixY = (srcMapY - srcGT.originY) / srcGT.pixelHeight;

                int sx = static_cast<int>(std::round(srcPixX));
                int sy = static_cast<int>(std::round(srcPixY));

                if (sx < 0 || sx >= srcCols || sy < 0 || sy >= srcRows) {
                    dst[dstIdx] = noData;
                    continue;
                }

                double val = src[static_cast<int64_t>(sy) * srcCols + sx];
                if (hasND && val == noData) {
                    dst[dstIdx] = noData;
                } else {
                    dst[dstIdx] = val;
                }
            }

            if (dstR % 100 == 0)
                reportProgress(static_cast<double>(dstR) / dstRows);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(LocalAffineModule)

} // namespace aplaceholder
