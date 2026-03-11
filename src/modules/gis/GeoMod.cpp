#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

namespace aplaceholder {

class GeoModModule : public Module {
public:
    QString name() const override { return "GEOMOD"; }
    QString description() const override {
        return "Land use change simulation. Compares two land cover maps to "
               "model change patterns and predicts a future map using simple "
               "proximity-weighted allocation for the specified quantity of change.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("time1", "Land cover map at time 1"),
            ParameterDef::file("time2", "Land cover map at time 2"),
            ParameterDef::integer("quantity", "Quantity of change (number of cells)", 100, 1, 1000000000,
                "Number of cells to convert in the predicted output"),
            ParameterDef::output("output", "Output predicted land cover map"),
        };
    }

    bool execute() override {
        auto map1 = GdalIO::read(parameter("time1").toString());
        if (!map1) { setError("Failed to read time 1 land cover map"); return false; }

        auto map2 = GdalIO::read(parameter("time2").toString());
        if (!map2) { setError("Failed to read time 2 land cover map"); return false; }

        int cols = map1->cols(), rows = map1->rows();

        if (map2->cols() != cols || map2->rows() != rows) {
            setError("Time 1 and time 2 maps must have the same dimensions");
            return false;
        }

        double noData = map1->noDataValue();
        bool hasND = map1->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int quantity = parameter("quantity").toInt();

        const auto& data1 = map1->data(0);
        const auto& data2 = map2->data(0);

        // Identify changed cells between time1 and time2 (locations of past change)
        reportProgress(0.0, "Analyzing change patterns...");
        std::vector<bool> changedCell(total, false);
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (data1[i] == noData || data2[i] == noData)) continue;
            if (static_cast<int>(data1[i]) != static_cast<int>(data2[i])) {
                changedCell[i] = true;
            }
        }

        // Compute proximity weight for each cell: inverse distance to nearest
        // changed cell. Cells closer to past change are more likely to change.
        reportProgress(0.2, "Computing proximity weights...");
        std::vector<double> proximity(total, 0.0);
        // Use a simple scan: for each candidate cell, find minimum distance
        // to any changed cell. To keep it tractable, collect changed cell coords.
        struct Coord { int r; int c; };
        std::vector<Coord> changedCoords;
        changedCoords.reserve(total / 10);
        for (int64_t i = 0; i < total; ++i) {
            if (changedCell[i]) {
                int r = static_cast<int>(i / cols);
                int c = static_cast<int>(i % cols);
                changedCoords.push_back({r, c});
            }
        }

        if (changedCoords.empty()) {
            setError("No change detected between time 1 and time 2 maps");
            return false;
        }

        // For efficiency with large rasters, subsample change coords if too many
        std::vector<Coord> sampleCoords;
        const int maxSample = 5000;
        if (static_cast<int>(changedCoords.size()) > maxSample) {
            int step = static_cast<int>(changedCoords.size()) / maxSample;
            for (int i = 0; i < static_cast<int>(changedCoords.size()); i += step) {
                sampleCoords.push_back(changedCoords[i]);
            }
        } else {
            sampleCoords = changedCoords;
        }

        // Identify candidate cells: cells in time2 that could potentially change
        // (i.e., cells that have NOT already changed, are valid, and match the
        //  dominant "from" class)
        // For simplicity, candidates are all non-nodata, non-changed cells
        struct Candidate {
            int64_t index;
            double weight;
        };
        std::vector<Candidate> candidates;

        reportProgress(0.4, "Evaluating candidate cells...");
        for (int64_t i = 0; i < total; ++i) {
            if (hasND && data2[i] == noData) continue;
            if (changedCell[i]) continue; // already changed

            int r = static_cast<int>(i / cols);
            int c = static_cast<int>(i % cols);

            // Compute minimum distance to a changed cell
            double minDist = std::numeric_limits<double>::max();
            for (const auto& cc : sampleCoords) {
                double dr = r - cc.r;
                double dc = c - cc.c;
                double dist = std::sqrt(dr * dr + dc * dc);
                if (dist < minDist) minDist = dist;
            }

            // Proximity weight: inverse distance (add 1 to avoid division by zero)
            double weight = 1.0 / (minDist + 1.0);
            candidates.push_back({i, weight});

            if (i % 1000000 == 0)
                reportProgress(0.4 + 0.3 * (static_cast<double>(i) / total));
        }

        if (candidates.empty()) {
            setError("No candidate cells available for allocation");
            return false;
        }

        // Sort candidates by weight descending (highest proximity first)
        reportProgress(0.7, "Ranking candidates by proximity...");
        std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& a, const Candidate& b) {
                return a.weight > b.weight;
            });

        // Determine "to" class: the most common destination class in observed changes
        std::map<int, int> toClassCounts;
        for (int64_t i = 0; i < total; ++i) {
            if (changedCell[i]) {
                int toVal = static_cast<int>(data2[i]);
                toClassCounts[toVal]++;
            }
        }
        int toClass = toClassCounts.begin()->first;
        int maxCnt = 0;
        for (const auto& pair : toClassCounts) {
            if (pair.second > maxCnt) {
                maxCnt = pair.second;
                toClass = pair.first;
            }
        }

        // Allocate: copy time2 as base, then convert top-weighted candidates
        reportProgress(0.8, "Allocating change...");
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(map1->geoTransform());
        output.setProjection(map1->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        for (int64_t i = 0; i < total; ++i) {
            out[i] = data2[i];
        }

        int allocated = 0;
        int maxAllocate = std::min(quantity, static_cast<int>(candidates.size()));
        for (int i = 0; i < maxAllocate; ++i) {
            out[candidates[i].index] = static_cast<double>(toClass);
            ++allocated;
        }

        reportProgress(1.0, "Writing output (" + QString::number(allocated) + " cells allocated)...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(GeoModModule)

} // namespace aplaceholder
