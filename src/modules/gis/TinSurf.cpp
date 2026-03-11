#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace aplaceholder {

class TinSurfModule : public Module {
public:
    QString name() const override { return "TINSURF"; }
    QString description() const override {
        return "Creates a raster surface from an existing TIN model. For each triangular "
               "facet, solves a planar equation and assigns interpolated values to pixels "
               "falling within that facet.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("tin_file", "Input TIN file",
                "TIN file produced by the TIN module"),
            ParameterDef::file("reference_raster", "Reference raster (extent/dimensions)",
                "Raster defining the output extent and resolution"),
            ParameterDef::output("output", "Output raster surface"),
        };
    }

    bool execute() override {
        QString tinPath = parameter("tin_file").toString();
        QString refPath = parameter("reference_raster").toString();
        QString outPath = parameter("output").toString();

        // Read TIN file
        struct Vertex { double x, y, z; };
        struct Triangle { int v[3]; };
        std::vector<Vertex> vertices;
        std::vector<Triangle> triangles;

        {
            std::ifstream ifs(tinPath.toStdString());
            if (!ifs.is_open()) {
                setError("Failed to open TIN file: " + tinPath);
                return false;
            }

            std::string line;
            int numVerts = 0, numTris = 0;
            enum Section { NONE, VERTS, TRIS } section = NONE;

            while (std::getline(ifs, line)) {
                if (line.empty()) continue;
                std::istringstream ss(line);
                std::string token;
                ss >> token;

                if (token == "NUM_VERTICES") {
                    ss >> numVerts;
                    vertices.reserve(numVerts);
                } else if (token == "NUM_TRIANGLES") {
                    ss >> numTris;
                    triangles.reserve(numTris);
                } else if (token == "VERTICES") {
                    section = VERTS;
                } else if (token == "TRIANGLES") {
                    section = TRIS;
                } else if (token == "END") {
                    break;
                } else if (section == VERTS) {
                    // index x y z
                    Vertex v;
                    ss >> v.x >> v.y >> v.z;
                    vertices.push_back(v);
                } else if (section == TRIS) {
                    // index v0 v1 v2
                    Triangle t;
                    ss >> t.v[0] >> t.v[1] >> t.v[2];
                    triangles.push_back(t);
                }
            }
        }

        if (vertices.empty() || triangles.empty()) {
            setError("TIN file contains no valid data");
            return false;
        }

        // Read reference raster metadata
        auto ref = GdalIO::readMetadata(refPath);
        if (!ref) {
            setError("Failed to read reference raster: " + refPath);
            return false;
        }

        int cols = ref->cols();
        int rows = ref->rows();
        double noData = -9999.0;

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(ref->geoTransform());
        output.setProjection(ref->projection());
        output.setNoDataValue(noData);

        auto& outData = output.data(0);
        // Initialize to NoData
        for (auto& v : outData) v = noData;

        reportProgress(0.0, "Rasterizing TIN...");

        // For each triangle, solve planar equation and fill pixels
        // H = A*x + B*y + C
        for (size_t ti = 0; ti < triangles.size(); ++ti) {
            const auto& tri = triangles[ti];
            const Vertex& v0 = vertices[tri.v[0]];
            const Vertex& v1 = vertices[tri.v[1]];
            const Vertex& v2 = vertices[tri.v[2]];

            // Solve for A, B, C: H = Ax + By + C
            // Using Cramer's rule on the 3x3 system
            double det = v0.x * (v1.y - v2.y) + v1.x * (v2.y - v0.y) + v2.x * (v0.y - v1.y);
            if (std::abs(det) < 1e-20) continue; // degenerate triangle

            double A = (v0.z * (v1.y - v2.y) + v1.z * (v2.y - v0.y) + v2.z * (v0.y - v1.y)) / det;
            double B = (v0.x * (v1.z - v2.z) + v1.x * (v2.z - v0.z) + v2.x * (v0.z - v1.z)) / det;
            double C = (v0.x * (v1.y * v2.z - v2.y * v1.z)
                      + v1.x * (v2.y * v0.z - v0.y * v2.z)
                      + v2.x * (v0.y * v1.z - v1.y * v0.z)) / det;

            // Bounding box of triangle in pixel coords
            double txMin = std::min({v0.x, v1.x, v2.x});
            double txMax = std::max({v0.x, v1.x, v2.x});
            double tyMin = std::min({v0.y, v1.y, v2.y});
            double tyMax = std::max({v0.y, v1.y, v2.y});

            int cMin, rMin, cMax, rMax;
            output.xyToColRow(txMin, tyMax, cMin, rMin); // note: y is inverted
            output.xyToColRow(txMax, tyMin, cMax, rMax);

            cMin = std::max(0, cMin - 1);
            rMin = std::max(0, rMin - 1);
            cMax = std::min(cols - 1, cMax + 1);
            rMax = std::min(rows - 1, rMax + 1);

            for (int r = rMin; r <= rMax; ++r) {
                for (int c = cMin; c <= cMax; ++c) {
                    double px, py;
                    output.colRowToXY(c, r, px, py);

                    // Barycentric coordinates to check containment
                    double d00 = (v1.x - v0.x) * (v1.x - v0.x) + (v1.y - v0.y) * (v1.y - v0.y);
                    double d01 = (v1.x - v0.x) * (v2.x - v0.x) + (v1.y - v0.y) * (v2.y - v0.y);
                    double d11 = (v2.x - v0.x) * (v2.x - v0.x) + (v2.y - v0.y) * (v2.y - v0.y);
                    double d20 = (px - v0.x) * (v1.x - v0.x) + (py - v0.y) * (v1.y - v0.y);
                    double d21 = (px - v0.x) * (v2.x - v0.x) + (py - v0.y) * (v2.y - v0.y);

                    double denom = d00 * d11 - d01 * d01;
                    if (std::abs(denom) < 1e-20) continue;

                    double u = (d11 * d20 - d01 * d21) / denom;
                    double w = (d00 * d21 - d01 * d20) / denom;

                    // Point is inside triangle if u >= 0, w >= 0, u + w <= 1
                    if (u >= -1e-10 && w >= -1e-10 && (u + w) <= 1.0 + 1e-10) {
                        double val = A * px + B * py + C;
                        outData[static_cast<size_t>(r) * cols + c] = val;
                    }
                }
            }

            if (ti % 100 == 0)
                reportProgress(0.9 * ti / triangles.size());
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outPath);
    }
};

REGISTER_MODULE(TinSurfModule)

} // namespace aplaceholder
