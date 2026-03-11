#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class PipedModule : public Module {
public:
    QString name() const override { return "PIPED"; }
    QString description() const override {
        return "Parallelepiped classifier. Assigns each pixel to a class if its values "
               "in all bands fall within the min-max range defined by the training signature.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::file("signature_file", "Signature file (CSV)",
                "CSV with rows: class_id, min_b1, max_b1, min_b2, max_b2, ..."),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::combo("overlap_mode", "Overlap handling",
                {"Unclassified", "Nearest centroid"}, 0,
                "How to handle pixels falling in overlapping class boxes"),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Read input bands
        // ------------------------------------------------------------------
        QString bandsParam = parameter("input_bands").toString();
        QStringList bandPaths = bandsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : bandPaths) p = p.trimmed();

        if (bandPaths.isEmpty()) {
            setError("No input bands specified");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> bandRasters;
        for (const auto& path : bandPaths) {
            auto r = GdalIO::read(path);
            if (!r) {
                setError("Failed to read band image: " + path);
                return false;
            }
            bandRasters.push_back(std::move(r));
        }

        int numBands = static_cast<int>(bandRasters.size());
        int cols = bandRasters[0]->cols();
        int rows = bandRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = bandRasters[0]->hasNoData();
        double noData = bandRasters[0]->noDataValue();

        // Verify all bands have same dimensions
        for (int b = 1; b < numBands; ++b) {
            if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                setError("All input bands must have the same dimensions");
                return false;
            }
        }

        int overlapMode = parameter("overlap_mode").toInt();

        // ------------------------------------------------------------------
        // 2. Parse signature file
        // ------------------------------------------------------------------
        // Format: class_id, min_b1, max_b1, min_b2, max_b2, ...
        QString sigPath = parameter("signature_file").toString();

        struct ClassBox {
            int classId;
            std::vector<double> minVal;  // per band
            std::vector<double> maxVal;  // per band
            std::vector<double> centroid; // midpoint per band (for overlap resolution)
        };

        std::vector<ClassBox> classes;

        {
            std::ifstream file(sigPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open signature file: " + sigPath);
                return false;
            }

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;

                std::istringstream iss(line);
                std::vector<double> values;
                std::string token;
                while (std::getline(iss, token, ',')) {
                    size_t start = token.find_first_not_of(" \t");
                    size_t end = token.find_last_not_of(" \t\r\n");
                    if (start != std::string::npos && end != std::string::npos)
                        values.push_back(std::stod(token.substr(start, end - start + 1)));
                }

                // Expected: 1 (class_id) + 2*numBands (min/max pairs)
                int expected = 1 + 2 * numBands;
                if (static_cast<int>(values.size()) < expected) continue;

                ClassBox box;
                box.classId = static_cast<int>(values[0]);
                box.minVal.resize(numBands);
                box.maxVal.resize(numBands);
                box.centroid.resize(numBands);

                for (int b = 0; b < numBands; ++b) {
                    box.minVal[b] = values[1 + 2 * b];
                    box.maxVal[b] = values[1 + 2 * b + 1];
                    box.centroid[b] = (box.minVal[b] + box.maxVal[b]) / 2.0;
                }

                classes.push_back(std::move(box));
            }
        }

        if (classes.empty()) {
            setError("No valid class signatures found in signature file");
            return false;
        }

        reportProgress(0.1, "Classifying pixels...");

        // ------------------------------------------------------------------
        // 3. Collect band data pointers
        // ------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // ------------------------------------------------------------------
        // 4. Create output classified image
        // ------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(bandRasters[0]->geoTransform());
        output.setProjection(bandRasters[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);

        // ------------------------------------------------------------------
        // 5. Classify each pixel
        // ------------------------------------------------------------------
        int numClasses = static_cast<int>(classes.size());

        for (int64_t i = 0; i < total; ++i) {
            // Check NoData
            if (hasND && (*bands[0])[i] == noData) {
                out[i] = 0;
                continue;
            }

            // Find which class boxes contain this pixel
            std::vector<int> matchingClasses;

            for (int c = 0; c < numClasses; ++c) {
                const auto& box = classes[c];
                bool inside = true;

                for (int b = 0; b < numBands; ++b) {
                    double val = (*bands[b])[i];
                    if (val < box.minVal[b] || val > box.maxVal[b]) {
                        inside = false;
                        break;
                    }
                }

                if (inside)
                    matchingClasses.push_back(c);
            }

            if (matchingClasses.empty()) {
                // Outside all boxes: unclassified
                out[i] = 0;
            } else if (matchingClasses.size() == 1) {
                // Unique assignment
                out[i] = classes[matchingClasses[0]].classId;
            } else {
                // Overlap
                if (overlapMode == 0) {
                    // Mark as unclassified
                    out[i] = 0;
                } else {
                    // Assign to nearest centroid (Euclidean distance)
                    double bestDist = std::numeric_limits<double>::max();
                    int bestClass = 0;

                    for (int ci : matchingClasses) {
                        const auto& box = classes[ci];
                        double distSq = 0.0;
                        for (int b = 0; b < numBands; ++b) {
                            double diff = (*bands[b])[i] - box.centroid[b];
                            distSq += diff * diff;
                        }
                        if (distSq < bestDist) {
                            bestDist = distSq;
                            bestClass = box.classId;
                        }
                    }

                    out[i] = bestClass;
                }
            }

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        // ------------------------------------------------------------------
        // 6. Write output
        // ------------------------------------------------------------------
        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PipedModule)

} // namespace aplaceholder
