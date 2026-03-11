#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <map>
#include <vector>

namespace aplaceholder {

class ConsensusModule : public Module {
public:
    QString name() const override { return "CONSENSUS"; }
    QString description() const override {
        return "Consensus analysis across multiple classification maps. "
               "Finds pixels where all maps agree (unanimous), a majority agree, "
               "or where the most common class wins (plurality).";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("inputs", "Input classification rasters (comma-separated paths)"),
            ParameterDef::output("output", "Output consensus raster"),
            ParameterDef::combo("method", "Consensus method",
                {"unanimous", "majority", "plurality"}, 0,
                "unanimous: all must agree; majority: >50% agree; plurality: most common class"),
        };
    }

    bool execute() override {
        QStringList inputPaths = parameter("inputs").toString().split(",", Qt::SkipEmptyParts);
        int numInputs = inputPaths.size();
        if (numInputs < 2) {
            setError("At least two input classification rasters are required");
            return false;
        }

        int methodIdx = parameter("method").toInt();
        // 0 = unanimous, 1 = majority, 2 = plurality

        reportProgress(0.0, "Loading input rasters...");

        std::vector<std::unique_ptr<Raster>> rasters(numInputs);
        for (int i = 0; i < numInputs; ++i) {
            rasters[i] = GdalIO::read(inputPaths[i].trimmed());
            if (!rasters[i]) {
                setError("Failed to read raster: " + inputPaths[i].trimmed());
                return false;
            }
        }

        int cols = rasters[0]->cols(), rows = rasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        for (int i = 1; i < numInputs; ++i) {
            if (rasters[i]->cols() != cols || rasters[i]->rows() != rows) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
        }

        reportProgress(0.1, "Computing consensus...");

        double outNoData = -9999.0;
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(rasters[0]->geoTransform());
        output.setProjection(rasters[0]->projection());
        output.setNoDataValue(outNoData);
        auto& out = output.data(0);

        int64_t consensusCount = 0;
        int64_t noConsensusCount = 0;

        for (int64_t p = 0; p < total; ++p) {
            // Collect class values from all inputs
            std::map<int, int> classCounts;
            bool valid = true;
            for (int i = 0; i < numInputs; ++i) {
                const auto& data = rasters[i]->data(0);
                if (rasters[i]->hasNoData() && data[p] == rasters[i]->noDataValue()) {
                    valid = false;
                    break;
                }
                int cls = static_cast<int>(data[p]);
                classCounts[cls]++;
            }

            if (!valid) {
                out[p] = outNoData;
                continue;
            }

            // Find the most common class
            int bestClass = 0;
            int bestCount = 0;
            for (const auto& kv : classCounts) {
                if (kv.second > bestCount) {
                    bestCount = kv.second;
                    bestClass = kv.first;
                }
            }

            bool hasConsensus = false;
            switch (methodIdx) {
                case 0: // unanimous
                    hasConsensus = (bestCount == numInputs);
                    break;
                case 1: // majority
                    hasConsensus = (bestCount > numInputs / 2);
                    break;
                case 2: // plurality
                    hasConsensus = true; // plurality always picks the most common
                    break;
            }

            if (hasConsensus) {
                out[p] = static_cast<double>(bestClass);
                ++consensusCount;
            } else {
                out[p] = outNoData;
                ++noConsensusCount;
            }

            if (p % 1000000 == 0)
                reportProgress(0.1 + 0.7 * static_cast<double>(p) / total);
        }

        reportProgress(0.85, "Writing output raster...");

        if (!GdalIO::write(output, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        static const char* methodNames[] = {"unanimous", "majority", "plurality"};
        reportProgress(1.0,
            QString("Consensus (%1): %2 consensus pixels, %3 no-consensus pixels")
                .arg(methodNames[methodIdx])
                .arg(consensusCount)
                .arg(noConsensusCount));

        return true;
    }
};

REGISTER_MODULE(ConsensusModule)

} // namespace aplaceholder
