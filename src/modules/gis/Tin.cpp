#include "Module.h"
#include "ModuleRegistry.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <array>

namespace aplaceholder {

struct TINVertex {
    double x, y, z;
};

struct TINTriangle {
    int v[3]; // vertex indices
};

class TinModule : public Module {
public:
    QString name() const override { return "TIN"; }
    QString description() const override {
        return "Creates a Delaunay triangulation from point data using incremental "
               "insertion. Outputs a TIN file with triangle vertex indices and coordinates.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("points_file", "Input point data (CSV: x,y,z)",
                "CSV file with columns x, y, z"),
            ParameterDef::output("output_tin_file", "Output TIN file",
                "Text file containing triangulation data"),
        };
    }

    bool execute() override {
        QString pointsPath = parameter("points_file").toString();
        QString outPath = parameter("output_tin_file").toString();

        // Read points
        std::vector<TINVertex> vertices;
        {
            std::ifstream ifs(pointsPath.toStdString());
            if (!ifs.is_open()) {
                setError("Failed to open points file: " + pointsPath);
                return false;
            }
            std::string line;
            std::getline(ifs, line); // skip header
            while (std::getline(ifs, line)) {
                if (line.empty()) continue;
                std::istringstream ss(line);
                TINVertex v;
                char sep;
                if (ss >> v.x >> sep >> v.y >> sep >> v.z) {
                    vertices.push_back(v);
                }
            }
        }

        if (vertices.size() < 3) {
            setError("Need at least 3 points for triangulation");
            return false;
        }

        reportProgress(0.0, "Building Delaunay triangulation...");

        // Compute bounding box
        double minX = vertices[0].x, maxX = vertices[0].x;
        double minY = vertices[0].y, maxY = vertices[0].y;
        for (const auto& v : vertices) {
            minX = std::min(minX, v.x); maxX = std::max(maxX, v.x);
            minY = std::min(minY, v.y); maxY = std::max(maxY, v.y);
        }

        double rangeX = maxX - minX;
        double rangeY = maxY - minY;
        double margin = std::max(rangeX, rangeY) * 10.0;

        // Add super-triangle vertices (3 points enclosing all data)
        int superBase = static_cast<int>(vertices.size());
        vertices.push_back({minX - margin, minY - margin, 0.0});
        vertices.push_back({maxX + margin * 2, minY - margin, 0.0});
        vertices.push_back({minX + rangeX / 2, maxY + margin * 2, 0.0});

        // Start with super-triangle
        std::vector<TINTriangle> triangles;
        triangles.push_back({{superBase, superBase + 1, superBase + 2}});

        // Incremental insertion
        int totalPts = superBase; // original point count
        for (int pi = 0; pi < totalPts; ++pi) {
            double px = vertices[pi].x;
            double py = vertices[pi].y;

            // Find triangles whose circumcircle contains the point
            std::vector<TINTriangle> badTriangles;
            std::vector<TINTriangle> goodTriangles;

            for (const auto& tri : triangles) {
                if (inCircumcircle(vertices, tri, px, py)) {
                    badTriangles.push_back(tri);
                } else {
                    goodTriangles.push_back(tri);
                }
            }

            // Find boundary polygon (edges of bad triangles not shared)
            struct Edge {
                int a, b;
                bool operator==(const Edge& o) const {
                    return (a == o.a && b == o.b) || (a == o.b && b == o.a);
                }
            };
            std::vector<Edge> polygon;

            for (const auto& tri : badTriangles) {
                Edge edges[3] = {
                    {tri.v[0], tri.v[1]},
                    {tri.v[1], tri.v[2]},
                    {tri.v[2], tri.v[0]},
                };
                for (const auto& e : edges) {
                    bool shared = false;
                    for (const auto& other : badTriangles) {
                        if (&other == &tri) continue;
                        // Check if this edge is in 'other'
                        Edge otherEdges[3] = {
                            {other.v[0], other.v[1]},
                            {other.v[1], other.v[2]},
                            {other.v[2], other.v[0]},
                        };
                        for (const auto& oe : otherEdges) {
                            if (e == oe) { shared = true; break; }
                        }
                        if (shared) break;
                    }
                    if (!shared) {
                        polygon.push_back(e);
                    }
                }
            }

            // Create new triangles from polygon edges to new point
            triangles = std::move(goodTriangles);
            for (const auto& e : polygon) {
                triangles.push_back({{pi, e.a, e.b}});
            }

            if (pi % 100 == 0)
                reportProgress(0.8 * pi / totalPts);
        }

        // Remove triangles that reference super-triangle vertices
        std::vector<TINTriangle> finalTriangles;
        for (const auto& tri : triangles) {
            bool usesSuperVertex = false;
            for (int i = 0; i < 3; ++i) {
                if (tri.v[i] >= superBase) {
                    usesSuperVertex = true;
                    break;
                }
            }
            if (!usesSuperVertex)
                finalTriangles.push_back(tri);
        }

        // Remove super-triangle vertices
        vertices.resize(superBase);

        reportProgress(0.9, "Writing TIN file...");

        // Write TIN file
        {
            std::ofstream ofs(outPath.toStdString());
            if (!ofs.is_open()) {
                setError("Failed to write output TIN file: " + outPath);
                return false;
            }

            ofs << "TIN_FILE\n";
            ofs << "NUM_VERTICES " << vertices.size() << "\n";
            ofs << "NUM_TRIANGLES " << finalTriangles.size() << "\n";
            ofs << "VERTICES\n";
            ofs.precision(12);
            for (size_t i = 0; i < vertices.size(); ++i) {
                ofs << i << " " << vertices[i].x << " " << vertices[i].y
                    << " " << vertices[i].z << "\n";
            }
            ofs << "TRIANGLES\n";
            for (size_t i = 0; i < finalTriangles.size(); ++i) {
                ofs << i << " " << finalTriangles[i].v[0] << " "
                    << finalTriangles[i].v[1] << " " << finalTriangles[i].v[2] << "\n";
            }
            ofs << "END\n";
        }

        reportProgress(1.0, "Done");
        return true;
    }

private:
    // Check if point (px,py) is inside the circumcircle of triangle tri
    bool inCircumcircle(const std::vector<TINVertex>& verts,
                        const TINTriangle& tri,
                        double px, double py) const {
        double ax = verts[tri.v[0]].x;
        double ay = verts[tri.v[0]].y;
        double bx = verts[tri.v[1]].x;
        double by = verts[tri.v[1]].y;
        double cx = verts[tri.v[2]].x;
        double cy = verts[tri.v[2]].y;

        double dax = ax - px;
        double day = ay - py;
        double dbx = bx - px;
        double dby = by - py;
        double dcx = cx - px;
        double dcy = cy - py;

        double det = dax * (dby * (dcx * dcx + dcy * dcy) - dcy * (dbx * dbx + dby * dby))
                   - day * (dbx * (dcx * dcx + dcy * dcy) - dcx * (dbx * dbx + dby * dby))
                   + (dax * dax + day * day) * (dbx * dcy - dby * dcx);

        // Ensure consistent orientation
        double cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
        if (cross < 0) det = -det;

        return det > 0;
    }
};

REGISTER_MODULE(TinModule)

} // namespace aplaceholder
