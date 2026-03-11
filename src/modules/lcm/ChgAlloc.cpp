#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <QFile>
#include <QTextStream>
#include <QStringList>

namespace aplaceholder {

class ChgAllocModule : public Module {
public:
    QString name() const override { return "CHGALLOC"; }
    QString description() const override {
        return "Change allocation for the Land Change Modeler. Allocates predicted "
               "change quantities to specific cells based on transition potential "
               "surfaces and a change demand matrix.";
    }
    QString category() const override { return "LCM"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("base_landcover", "Base land cover raster",
                "Current land cover classification raster"),
            ParameterDef::file("transition_potentials", "Transition potential rasters (comma-separated)",
                "Comma-separated list of transition potential surface raster files, one per transition"),
            ParameterDef::file("change_demand_file", "Change demand file (CSV)",
                "CSV file with columns: from_class, to_class, num_cells_to_change"),
            ParameterDef::output("output", "Output predicted land cover raster"),
        };
    }

    bool execute() override {
        // Read base land cover
        auto baseLc = GdalIO::read(parameter("base_landcover").toString());
        if (!baseLc) {
            setError("Failed to read base land cover raster");
            return false;
        }

        int cols = baseLc->cols();
        int rows = baseLc->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double lcND = baseLc->noDataValue();
        bool lcHasND = baseLc->hasNoData();
        const auto& lcData = baseLc->data(0);

        // Read transition potential rasters
        QStringList tpPaths = parameter("transition_potentials").toString()
                                  .split(",", Qt::SkipEmptyParts);
        int numTP = tpPaths.size();
        if (numTP == 0) {
            setError("No transition potential rasters specified");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> tpRasters;
        for (int t = 0; t < numTP; ++t) {
            auto tp = GdalIO::read(tpPaths[t].trimmed());
            if (!tp) {
                setError("Failed to read transition potential raster: " +
                         tpPaths[t].trimmed());
                return false;
            }
            if (tp->cols() != cols || tp->rows() != rows) {
                setError("Transition potential dimensions do not match base land cover: " +
                         tpPaths[t].trimmed());
                return false;
            }
            tpRasters.push_back(std::move(tp));
        }

        reportProgress(0.1, "Reading change demand matrix...");

        // Parse change demand CSV
        // Expected format: from_class, to_class, num_cells
        struct TransitionDemand {
            int fromClass;
            int toClass;
            int64_t numCells;
            int tpIndex; // index into tpRasters
        };
        std::vector<TransitionDemand> demands;

        QString demandPath = parameter("change_demand_file").toString();
        QFile demandFile(demandPath);
        if (!demandFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
            setError("Failed to open change demand file: " + demandPath);
            return false;
        }

        QTextStream in(&demandFile);
        int tpIdx = 0;
        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty() || line.startsWith("#") || line.startsWith("from"))
                continue; // skip comments and header

            QStringList parts = line.split(",");
            if (parts.size() < 3) continue;

            TransitionDemand td;
            td.fromClass = parts[0].trimmed().toInt();
            td.toClass = parts[1].trimmed().toInt();
            td.numCells = parts[2].trimmed().toLongLong();
            td.tpIndex = tpIdx;
            if (tpIdx < numTP)
                tpIdx++;

            if (td.numCells > 0)
                demands.push_back(td);
        }
        demandFile.close();

        if (demands.empty()) {
            setError("No valid change demands found in demand file");
            return false;
        }

        // Validate tpIndex mapping
        for (auto& td : demands) {
            if (td.tpIndex >= numTP) {
                td.tpIndex = numTP - 1; // clamp to last available
            }
        }

        reportProgress(0.2, "Allocating changes...");

        // Copy base land cover to output
        std::vector<double> outputData(total);
        for (int64_t i = 0; i < total; ++i)
            outputData[i] = lcData[i];

        // Track which pixels have been changed (can only change once)
        std::vector<bool> changed(total, false);

        // Process each transition demand
        int demandIdx = 0;
        for (auto& td : demands) {
            reportProgress(0.2 + 0.7 * static_cast<double>(demandIdx) / demands.size(),
                           "Allocating transition " +
                               QString::number(td.fromClass) + " -> " +
                               QString::number(td.toClass));

            // Collect candidate pixels: must be of fromClass and not yet changed
            struct Candidate {
                int64_t pixelIdx;
                double potential;
            };
            std::vector<Candidate> candidates;

            const auto& tpData = tpRasters[td.tpIndex]->data(0);
            bool tpHasND = tpRasters[td.tpIndex]->hasNoData();
            double tpND = tpRasters[td.tpIndex]->noDataValue();

            for (int64_t i = 0; i < total; ++i) {
                if (changed[i]) continue;
                if (lcHasND && lcData[i] == lcND) continue;

                int currentClass = static_cast<int>(outputData[i]);
                if (currentClass != td.fromClass) continue;

                double pot = tpData[i];
                if (tpHasND && pot == tpND) continue;
                if (pot <= 0.0) continue;

                candidates.push_back({i, pot});
            }

            // Sort candidates by transition potential descending (highest potential first)
            std::sort(candidates.begin(), candidates.end(),
                      [](const Candidate& a, const Candidate& b) {
                          return a.potential > b.potential;
                      });

            // Allocate the demanded number of cells
            int64_t allocated = 0;
            int64_t toAllocate = std::min(td.numCells,
                                           static_cast<int64_t>(candidates.size()));

            for (int64_t j = 0; j < toAllocate; ++j) {
                int64_t pIdx = candidates[j].pixelIdx;
                outputData[pIdx] = static_cast<double>(td.toClass);
                changed[pIdx] = true;
                allocated++;
            }

            demandIdx++;
        }

        reportProgress(0.95, "Writing output...");

        // Create output raster
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(baseLc->geoTransform());
        output.setProjection(baseLc->projection());
        if (lcHasND)
            output.setNoDataValue(lcND);

        auto& outRasterData = output.data(0);
        for (int64_t i = 0; i < total; ++i)
            outRasterData[i] = outputData[i];

        reportProgress(1.0, "Change allocation complete.");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(ChgAllocModule)

} // namespace aplaceholder
