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

class MinDistModule : public Module {
public:
    QString name() const override { return "MINDIST"; }
    QString description() const override {
        return "Minimum Distance to Means classifier. Assigns each pixel to the class "
               "whose mean signature vector is closest in spectral space.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input multi-band image"),
            ParameterDef::file("signature_file", "Signature file (CSV)"),
            ParameterDef::output("output", "Output classified image"),
            ParameterDef::combo("distance_type", "Distance type",
                {"Euclidean", "Standardized"}, 0,
                "Euclidean uses raw band-space distance; Standardized divides by "
                "class standard deviation to account for differing class variability"),
            ParameterDef::real("max_distance", "Maximum distance threshold",
                0.0, 0.0, 999999.0,
                "Pixels farther than this from all class means are left unclassified. "
                "Set to 0 to disable threshold (classify all pixels)."),
        };
    }

    bool execute() override {
        // ------------------------------------------------------------------
        // 1. Read input raster
        // ------------------------------------------------------------------
        auto raster = GdalIO::read(parameter("input").toString());
        if (!raster) {
            setError("Failed to read input image");
            return false;
        }

        int cols = raster->cols(), rows = raster->rows();
        int numBands = raster->bands();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = raster->hasNoData();
        double noData = raster->noDataValue();

        if (numBands < 1) {
            setError("Input image must have at least one band");
            return false;
        }

        // ------------------------------------------------------------------
        // 2. Read parameters
        // ------------------------------------------------------------------
        int distanceType = parameter("distance_type").toInt();  // 0 = Euclidean, 1 = Standardized
        double maxDist = parameter("max_distance").toDouble();
        bool useThreshold = (maxDist > 0.0);

        // ------------------------------------------------------------------
        // 3. Parse signature file
        // ------------------------------------------------------------------
        // Expected CSV format:
        //   class_id, mean_b1, mean_b2, ..., mean_bN, stddev_b1, stddev_b2, ..., stddev_bN
        // One row per class. The stddev columns are required for Standardized
        // distance, but also read for Euclidean (ignored if present).
        // Minimal format (Euclidean only):
        //   class_id, mean_b1, mean_b2, ..., mean_bN

        QString sigPath = parameter("signature_file").toString();

        struct ClassSig {
            int classId;
            std::vector<double> mean;    // length = numBands
            std::vector<double> stddev;  // length = numBands (for standardized distance)
        };

        std::vector<ClassSig> classes;

        {
            std::ifstream file(sigPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open signature file: " + sigPath);
                return false;
            }

            std::string line;
            while (std::getline(file, line)) {
                // Skip empty lines and comments
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

                // Minimum: 1 (class_id) + numBands (means)
                int expectedMin = 1 + numBands;
                if (static_cast<int>(values.size()) < expectedMin) continue;

                ClassSig sig;
                sig.classId = static_cast<int>(values[0]);

                sig.mean.resize(numBands);
                for (int b = 0; b < numBands; ++b)
                    sig.mean[b] = values[1 + b];

                // Read standard deviations if present
                sig.stddev.resize(numBands, 1.0);  // default to 1.0 (no scaling)
                int expectedWithStddev = 1 + numBands + numBands;
                if (static_cast<int>(values.size()) >= expectedWithStddev) {
                    for (int b = 0; b < numBands; ++b) {
                        double sd = values[1 + numBands + b];
                        // Guard against zero or negative stddev
                        sig.stddev[b] = (sd > 1e-12) ? sd : 1e-12;
                    }
                }

                classes.push_back(std::move(sig));
            }
        }

        if (classes.empty()) {
            setError("No valid class signatures found in signature file");
            return false;
        }

        // If standardized distance was requested but no stddevs were provided,
        // warn but continue (stddevs default to 1.0, which reduces to Euclidean).
        // This is intentional — it degrades gracefully.

        reportProgress(0.1, "Classifying pixels...");

        // ------------------------------------------------------------------
        // 4. Collect band data pointers
        // ------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &raster->data(b);

        // ------------------------------------------------------------------
        // 5. Create output classified image
        // ------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Int32);
        output.setGeoTransform(raster->geoTransform());
        output.setProjection(raster->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);

        // ------------------------------------------------------------------
        // 6. Classify each pixel
        // ------------------------------------------------------------------
        // For Euclidean distance:
        //   d(x, c) = sqrt( sum_b (x_b - mean_c_b)^2 )
        //
        // For Standardized distance:
        //   d(x, c) = sqrt( sum_b ((x_b - mean_c_b) / stddev_c_b)^2 )
        //
        // We compare squared distances to avoid the sqrt in the inner loop
        // (sqrt is monotonic, so argmin is preserved). We only compute the
        // actual distance when checking the threshold.

        int numClasses = static_cast<int>(classes.size());

        for (int64_t i = 0; i < total; ++i) {
            // Check for NoData — use first band as indicator
            if (hasND && (*bands[0])[i] == noData) {
                out[i] = 0;  // unclassified
                continue;
            }

            double bestDistSq = std::numeric_limits<double>::max();
            int bestClass = 0;  // 0 = unclassified

            for (int c = 0; c < numClasses; ++c) {
                const auto& cls = classes[c];
                double distSq = 0.0;

                if (distanceType == 0) {
                    // Euclidean distance (squared)
                    for (int b = 0; b < numBands; ++b) {
                        double diff = (*bands[b])[i] - cls.mean[b];
                        distSq += diff * diff;
                    }
                } else {
                    // Standardized distance (squared)
                    for (int b = 0; b < numBands; ++b) {
                        double diff = ((*bands[b])[i] - cls.mean[b]) / cls.stddev[b];
                        distSq += diff * diff;
                    }
                }

                if (distSq < bestDistSq) {
                    bestDistSq = distSq;
                    bestClass = cls.classId;
                }
            }

            // Apply threshold — compare actual distance (need sqrt now)
            if (useThreshold) {
                double bestDist = std::sqrt(bestDistSq);
                if (bestDist > maxDist) {
                    bestClass = 0;  // unclassified
                }
            }

            out[i] = bestClass;

            if (i % 1000000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(i) / total);
        }

        // ------------------------------------------------------------------
        // 7. Write output
        // ------------------------------------------------------------------
        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(MinDistModule)

} // namespace aplaceholder
