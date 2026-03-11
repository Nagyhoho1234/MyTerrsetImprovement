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

class HyperMinModule : public Module {
public:
    QString name() const override { return "HYPERMIN"; }
    QString description() const override {
        return "Minimum-distance classifier for hyperspectral data. For each pixel, finds the "
               "endmember with the minimum Euclidean spectral distance. Outputs class assignment "
               "and distance map.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("bands", "Input band images (comma-separated)",
                "Comma-separated list of input band image file paths"),
            ParameterDef::file("endmember_file", "Endmember file (CSV)",
                "CSV: each row is an endmember, columns are band values"),
            ParameterDef::output("output_class", "Output class image",
                "Output raster with class assignments"),
            ParameterDef::output("output_distance", "Output distance image",
                "Output raster with minimum spectral distance per pixel"),
        };
    }

    bool execute() override {
        // --------------------------------------------------------------------
        // 1. Read input bands
        // --------------------------------------------------------------------
        QString bandsParam = parameter("bands").toString();
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

        for (int b = 1; b < numBands; ++b) {
            if (bandRasters[b]->cols() != cols || bandRasters[b]->rows() != rows) {
                setError("All input bands must have the same dimensions");
                return false;
            }
        }

        // --------------------------------------------------------------------
        // 2. Read endmember file (CSV)
        // --------------------------------------------------------------------
        QString emPath = parameter("endmember_file").toString();
        std::vector<std::vector<double>> endmembers;

        {
            std::ifstream file(emPath.toStdString());
            if (!file.is_open()) {
                setError("Failed to open endmember file: " + emPath);
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

                if (static_cast<int>(values.size()) >= numBands) {
                    values.resize(numBands);
                    endmembers.push_back(std::move(values));
                }
            }
        }

        int numEM = static_cast<int>(endmembers.size());
        if (numEM == 0) {
            setError("No valid endmembers found in endmember file");
            return false;
        }

        reportProgress(0.05, "Classifying pixels by minimum spectral distance...");

        // --------------------------------------------------------------------
        // 3. Collect band data pointers
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> bands(numBands);
        for (int b = 0; b < numBands; ++b)
            bands[b] = &bandRasters[b]->data(0);

        // --------------------------------------------------------------------
        // 4. Create output rasters
        // --------------------------------------------------------------------
        Raster outputClass(cols, rows, 1, DataType::Float32);
        outputClass.setGeoTransform(bandRasters[0]->geoTransform());
        outputClass.setProjection(bandRasters[0]->projection());
        outputClass.setNoDataValue(noData);
        auto& classData = outputClass.data(0);

        Raster outputDist(cols, rows, 1, DataType::Float64);
        outputDist.setGeoTransform(bandRasters[0]->geoTransform());
        outputDist.setProjection(bandRasters[0]->projection());
        outputDist.setNoDataValue(noData);
        auto& distData = outputDist.data(0);

        // --------------------------------------------------------------------
        // 5. For each pixel, find endmember with minimum Euclidean distance
        // --------------------------------------------------------------------
        std::vector<double> pixelVals(numBands);

        for (int64_t idx = 0; idx < total; ++idx) {
            if (hasND && (*bands[0])[idx] == noData) {
                classData[idx] = noData;
                distData[idx] = noData;
                continue;
            }

            for (int b = 0; b < numBands; ++b)
                pixelVals[b] = (*bands[b])[idx];

            double minDist = std::numeric_limits<double>::max();
            int bestClass = 0;

            for (int e = 0; e < numEM; ++e) {
                double dist = 0.0;
                for (int b = 0; b < numBands; ++b) {
                    double diff = pixelVals[b] - endmembers[e][b];
                    dist += diff * diff;
                }
                dist = std::sqrt(dist);

                if (dist < minDist) {
                    minDist = dist;
                    bestClass = e + 1; // 1-based
                }
            }

            classData[idx] = static_cast<double>(bestClass);
            distData[idx] = minDist;

            if (idx % 1000000 == 0)
                reportProgress(0.05 + 0.85 * static_cast<double>(idx) / total);
        }

        // --------------------------------------------------------------------
        // 6. Write outputs
        // --------------------------------------------------------------------
        reportProgress(0.9, "Writing output rasters...");

        QString classPath = parameter("output_class").toString();
        if (!GdalIO::write(outputClass, classPath)) {
            setError("Failed to write class output: " + classPath);
            return false;
        }

        QString distPath = parameter("output_distance").toString();
        if (!GdalIO::write(outputDist, distPath)) {
            setError("Failed to write distance output: " + distPath);
            return false;
        }

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(HyperMinModule)

} // namespace aplaceholder
