#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <memory>

namespace aplaceholder {

class SpeciesDistributionModule : public Module {
public:
    QString name() const override { return "SPECIES_DISTRIBUTION"; }
    QString description() const override {
        return "Species distribution modeling using a bioclimatic envelope (BIOCLIM) "
               "approach. Extracts environmental variable values at species presence "
               "locations, computes a percentile-based envelope, and scores every cell "
               "by the proportion of variables that fall within the envelope.";
    }
    QString category() const override { return "Habitat & Biodiversity"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("presence_points", "Presence points file (CSV with x,y)",
                "CSV file containing species occurrence locations. Must have at least "
                "two columns: x (longitude/easting) and y (latitude/northing). "
                "First row is treated as a header if non-numeric."),
            ParameterDef::file("env_rasters", "Environmental rasters (comma-separated)",
                "Comma-separated list of environmental variable rasters "
                "(e.g. bioclimatic variables, elevation, precipitation)."),
            ParameterDef::output("output", "Output suitability raster",
                "Output raster where each cell contains the proportion of "
                "environmental variables within the species envelope (0-1)."),
            ParameterDef::real("percentile", "Envelope percentile threshold",
                5.0, 0.0, 50.0,
                "Lower percentile for the bioclimatic envelope. The envelope is "
                "defined as [P, 100-P] where P is this value. E.g., 5 gives the "
                "5th-95th percentile range. Lower values produce wider envelopes."),
        };
    }

    bool execute() override {
        // --- Parse parameters --------------------------------------------------
        const QString presencePath   = parameter("presence_points").toString();
        const QString envRastersStr  = parameter("env_rasters").toString();
        const QString outputPath     = parameter("output").toString();
        const double  percentile     = parameter("percentile").toDouble();

        if (presencePath.isEmpty() || envRastersStr.isEmpty() || outputPath.isEmpty()) {
            setError("Presence points file, environmental rasters, and output "
                     "path are all required.");
            return false;
        }

        if (percentile < 0.0 || percentile >= 50.0) {
            setError("Percentile threshold must be in the range [0, 50).");
            return false;
        }

        QStringList envPaths = envRastersStr.split(",", Qt::SkipEmptyParts);
        for (auto& p : envPaths) p = p.trimmed();

        const int numVars = envPaths.size();
        if (numVars == 0) {
            setError("At least one environmental raster is required.");
            return false;
        }

        // --- Read presence points from CSV ------------------------------------
        reportProgress(0.0, "Reading presence points...");

        std::vector<double> ptX, ptY;
        if (!readPresenceCSV(presencePath, ptX, ptY)) {
            return false; // error already set
        }

        if (ptX.empty()) {
            setError("No valid presence points found in " + presencePath);
            return false;
        }

        const int numPoints = static_cast<int>(ptX.size());

        // --- Load environmental rasters ---------------------------------------
        reportProgress(0.05, "Loading environmental rasters...");

        std::vector<std::unique_ptr<Raster>> envRasters(numVars);
        envRasters[0] = GdalIO::read(envPaths[0]);
        if (!envRasters[0]) {
            setError("Failed to read environmental raster: " + envPaths[0]);
            return false;
        }

        const int cols = envRasters[0]->cols();
        const int rows = envRasters[0]->rows();

        for (int v = 1; v < numVars; ++v) {
            envRasters[v] = GdalIO::read(envPaths[v]);
            if (!envRasters[v]) {
                setError("Failed to read environmental raster: " + envPaths[v]);
                return false;
            }
            if (envRasters[v]->cols() != cols || envRasters[v]->rows() != rows) {
                setError(QString("Raster dimension mismatch: '%1' is %2x%3 but "
                                 "expected %4x%5.")
                             .arg(envPaths[v])
                             .arg(envRasters[v]->cols())
                             .arg(envRasters[v]->rows())
                             .arg(cols)
                             .arg(rows));
                return false;
            }
        }

        // --- Extract environmental values at presence points -------------------
        reportProgress(0.15, "Extracting environmental values at presence locations...");

        // envValues[var][point]
        std::vector<std::vector<double>> envValues(numVars);
        for (int v = 0; v < numVars; ++v) {
            envValues[v].reserve(numPoints);
        }

        int validPoints = 0;
        for (int p = 0; p < numPoints; ++p) {
            int col, row;
            envRasters[0]->xyToColRow(ptX[p], ptY[p], col, row);

            // Check bounds
            if (col < 0 || col >= cols || row < 0 || row >= rows) {
                continue; // point outside raster extent
            }

            // Check for nodata at this location in any variable
            bool anyNoData = false;
            std::vector<double> vals(numVars);
            for (int v = 0; v < numVars; ++v) {
                vals[v] = envRasters[v]->value(col, row, 0);
                if (envRasters[v]->hasNoData() &&
                    vals[v] == envRasters[v]->noDataValue()) {
                    anyNoData = true;
                    break;
                }
            }

            if (anyNoData) continue;

            for (int v = 0; v < numVars; ++v) {
                envValues[v].push_back(vals[v]);
            }
            validPoints++;
        }

        if (validPoints < 2) {
            setError(QString("Insufficient valid presence points (%1). At least 2 "
                             "points with valid environmental data are required.")
                         .arg(validPoints));
            return false;
        }

        // --- Compute percentile envelopes per variable -------------------------
        reportProgress(0.25, "Computing bioclimatic envelopes...");

        struct Envelope {
            double lower;
            double upper;
        };
        std::vector<Envelope> envelopes(numVars);

        for (int v = 0; v < numVars; ++v) {
            std::vector<double> sorted = envValues[v];
            std::sort(sorted.begin(), sorted.end());

            const int n = static_cast<int>(sorted.size());

            // Compute lower and upper percentile values using linear interpolation
            envelopes[v].lower = interpolatePercentile(sorted, percentile);
            envelopes[v].upper = interpolatePercentile(sorted, 100.0 - percentile);
        }

        // --- Score every cell --------------------------------------------------
        reportProgress(0.3, "Scoring cells against envelopes...");

        const double noData = envRasters[0]->noDataValue();

        auto output = std::make_unique<Raster>(cols, rows, 1, DataType::Float32);
        output->setGeoTransform(envRasters[0]->geoTransform());
        output->setProjection(envRasters[0]->projection());
        output->setNoDataValue(noData);
        output->allocate();

        auto& outData = output->data(0);

        // Collect data vector pointers for fast access
        std::vector<const std::vector<double>*> envData(numVars);
        for (int v = 0; v < numVars; ++v) {
            envData[v] = &envRasters[v]->data(0);
        }

        for (int r = 0; r < rows; ++r) {
            if (r % 100 == 0) {
                reportProgress(0.3 + 0.65 * static_cast<double>(r) / rows,
                               QString("Scoring row %1 / %2").arg(r + 1).arg(rows));
            }

            for (int c = 0; c < cols; ++c) {
                const int64_t idx = static_cast<int64_t>(r) * cols + c;

                // Check for nodata in any layer
                bool isNoData = false;
                for (int v = 0; v < numVars; ++v) {
                    double val = (*envData[v])[idx];
                    if (envRasters[v]->hasNoData() &&
                        val == envRasters[v]->noDataValue()) {
                        isNoData = true;
                        break;
                    }
                }

                if (isNoData) {
                    outData[idx] = noData;
                    continue;
                }

                // Count how many variables fall within their envelope
                int withinCount = 0;
                for (int v = 0; v < numVars; ++v) {
                    double val = (*envData[v])[idx];
                    if (val >= envelopes[v].lower && val <= envelopes[v].upper) {
                        withinCount++;
                    }
                }

                // Score = proportion of variables within envelope
                outData[idx] = static_cast<double>(withinCount) / numVars;
            }
        }

        // --- Write output ------------------------------------------------------
        reportProgress(0.95, "Writing output raster...");
        if (!GdalIO::write(*output, outputPath)) {
            setError("Failed to write output raster: " + outputPath);
            return false;
        }

        reportProgress(1.0, "Complete");
        return true;
    }

private:
    // Read a CSV file with x,y columns. Supports header rows (auto-detected).
    // Returns false on error (with m_lastError set).
    bool readPresenceCSV(const QString& path, std::vector<double>& xs,
                         std::vector<double>& ys) {
        std::ifstream file(path.toStdString());
        if (!file.is_open()) {
            setError("Cannot open presence points file: " + path);
            return false;
        }

        std::string line;
        bool firstLine = true;

        while (std::getline(file, line)) {
            // Skip empty lines
            if (line.empty()) continue;

            // Trim carriage return if present
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }

            // Detect delimiter: comma, tab, or space
            char delim = ',';
            if (line.find(',') == std::string::npos) {
                if (line.find('\t') != std::string::npos) {
                    delim = '\t';
                } else {
                    delim = ' ';
                }
            }

            std::stringstream ss(line);
            std::string tokenX, tokenY;

            if (delim == ',') {
                std::getline(ss, tokenX, ',');
                std::getline(ss, tokenY, ',');
            } else if (delim == '\t') {
                std::getline(ss, tokenX, '\t');
                std::getline(ss, tokenY, '\t');
            } else {
                ss >> tokenX >> tokenY;
            }

            // Try to parse as doubles
            try {
                double x = std::stod(tokenX);
                double y = std::stod(tokenY);
                xs.push_back(x);
                ys.push_back(y);
            } catch (...) {
                // If first line fails to parse, treat as header and skip
                if (firstLine) {
                    firstLine = false;
                    continue;
                }
                // Otherwise skip malformed lines silently
            }

            firstLine = false;
        }

        return true;
    }

    // Linear interpolation of a percentile from a sorted vector.
    // percentile is in [0, 100].
    static double interpolatePercentile(const std::vector<double>& sorted,
                                        double percentile) {
        if (sorted.empty()) return 0.0;
        if (sorted.size() == 1) return sorted[0];

        const int n = static_cast<int>(sorted.size());

        // Use the C = 1 variant (Excel PERCENTILE.INC equivalent)
        double rank = (percentile / 100.0) * (n - 1);
        int lower = static_cast<int>(std::floor(rank));
        int upper = static_cast<int>(std::ceil(rank));

        if (lower < 0) lower = 0;
        if (upper >= n) upper = n - 1;

        if (lower == upper) return sorted[lower];

        double frac = rank - lower;
        return sorted[lower] + frac * (sorted[upper] - sorted[lower]);
    }
};

REGISTER_MODULE(SpeciesDistributionModule)

} // namespace aplaceholder
