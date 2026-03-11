#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cstdint>

namespace aplaceholder {

class MolaModule : public Module {
public:
    QString name() const override { return "MOLA"; }
    QString description() const override {
        return "Multi-objective land allocation. "
               "Allocates land among competing objectives based on "
               "suitability surfaces and target area requirements using "
               "iterative conflict resolution with rank-based assignment.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            {"objectives", "Objective suitability images (comma-separated paths)",
             ParameterDef::String, {}, {}, 0, 0,
             "Comma-separated paths to suitability images (0-255) for each objective", true},
            {"areas", "Target areas (comma-separated pixel counts)", ParameterDef::String,
             {}, {}, 0, 0, "Target number of pixels to allocate for each objective", true},
            ParameterDef::output("output", "Output land allocation image"),
        };
    }

    bool execute() override {
        QString outputPath = parameter("output").toString();

        // Parse objective file paths
        std::vector<QString> objPaths = parsePaths(parameter("objectives").toString());
        std::vector<int64_t> targetAreas = parseInts(parameter("areas").toString());

        size_t nObj = objPaths.size();
        if (nObj == 0) {
            setError("At least one objective suitability image is required");
            return false;
        }
        if (targetAreas.size() != nObj) {
            setError("Number of target areas (" + QString::number(targetAreas.size()) +
                     ") must match number of objectives (" + QString::number(nObj) + ")");
            return false;
        }

        // Load suitability rasters
        std::vector<std::unique_ptr<Raster>> suitMaps;
        for (const auto& path : objPaths) {
            auto r = GdalIO::read(path.trimmed());
            if (!r) {
                setError("Failed to read suitability image: " + path);
                return false;
            }
            suitMaps.push_back(std::move(r));
        }

        int cols = suitMaps[0]->cols(), rows = suitMaps[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Validate dimensions match
        for (size_t i = 1; i < nObj; ++i) {
            if (suitMaps[i]->cols() != cols || suitMaps[i]->rows() != rows) {
                setError("Objective image " + objPaths[i] + " has mismatched dimensions");
                return false;
            }
        }

        // Validate target areas
        int64_t totalTarget = 0;
        for (size_t i = 0; i < nObj; ++i) {
            if (targetAreas[i] <= 0) {
                setError("Target area for objective " + QString::number(i + 1) + " must be positive");
                return false;
            }
            totalTarget += targetAreas[i];
        }
        if (totalTarget > total) {
            setError("Sum of target areas exceeds total number of pixels");
            return false;
        }

        reportProgress(0.0, "Computing suitability ranks...");

        // For each objective, create sorted index arrays (descending by suitability)
        // rank[obj][pixel] = rank of that pixel for that objective (0 = best)
        std::vector<std::vector<int64_t>> sortedIndices(nObj);
        std::vector<std::vector<int64_t>> ranks(nObj);

        for (size_t obj = 0; obj < nObj; ++obj) {
            const auto& data = suitMaps[obj]->data(0);

            // Create index array and sort by suitability descending
            sortedIndices[obj].resize(total);
            std::iota(sortedIndices[obj].begin(), sortedIndices[obj].end(), 0);
            std::sort(sortedIndices[obj].begin(), sortedIndices[obj].end(),
                [&data](int64_t a, int64_t b) {
                    return data[a] > data[b]; // descending
                });

            // Compute rank for each pixel
            ranks[obj].resize(total);
            for (int64_t r = 0; r < total; ++r) {
                ranks[obj][sortedIndices[obj][r]] = r;
            }

            reportProgress(0.05 + 0.15 * static_cast<double>(obj + 1) / nObj,
                "Ranked objective " + QString::number(obj + 1));
        }

        reportProgress(0.2, "Starting iterative allocation...");

        // allocation[i] = 0 means unallocated, obj+1 means allocated to objective obj
        std::vector<int> allocation(total, 0);
        std::vector<int64_t> allocated(nObj, 0); // current allocation count per objective

        // Track the current decision line (cutoff rank) for each objective
        std::vector<int64_t> cutoff(nObj);
        for (size_t obj = 0; obj < nObj; ++obj) {
            cutoff[obj] = targetAreas[obj];
        }

        const int MAX_ITERATIONS = 200;
        int iteration = 0;

        while (iteration < MAX_ITERATIONS) {
            ++iteration;

            // Reset allocation
            std::fill(allocation.begin(), allocation.end(), 0);
            std::fill(allocated.begin(), allocated.end(), 0);

            // Step 1: Each objective claims its top-ranked pixels up to the cutoff
            // Use a claim count to detect conflicts
            std::vector<std::vector<int>> claims(total); // which objectives claim each pixel

            for (size_t obj = 0; obj < nObj; ++obj) {
                for (int64_t r = 0; r < cutoff[obj] && r < total; ++r) {
                    int64_t pixIdx = sortedIndices[obj][r];
                    claims[pixIdx].push_back(static_cast<int>(obj));
                }
            }

            // Step 2: Assign non-conflict pixels; identify conflicts
            int64_t conflictCount = 0;
            for (int64_t i = 0; i < total; ++i) {
                if (claims[i].size() == 1) {
                    int obj = claims[i][0];
                    allocation[i] = obj + 1;
                    allocated[obj]++;
                } else if (claims[i].size() > 1) {
                    ++conflictCount;

                    // Step 3: Resolve conflict - assign to the objective with the
                    // minimum (best) rank at this pixel.
                    // This implements the minimum-distance-to-ideal-point rule:
                    // the objective that would lose the most by not getting this pixel wins.
                    int bestObj = claims[i][0];
                    int64_t bestRank = ranks[bestObj][i];

                    for (size_t c = 1; c < claims[i].size(); ++c) {
                        int obj = claims[i][c];
                        int64_t objRank = ranks[obj][i];
                        if (objRank < bestRank) {
                            bestRank = objRank;
                            bestObj = obj;
                        }
                    }

                    allocation[i] = bestObj + 1;
                    allocated[bestObj]++;
                }
            }

            // Step 4: Check if all objectives have met their targets
            bool allMet = true;
            for (size_t obj = 0; obj < nObj; ++obj) {
                if (allocated[obj] < targetAreas[obj]) {
                    allMet = false;
                    break;
                }
            }

            if (allMet && conflictCount == 0) break;

            // Step 5: Adjust cutoff lines for objectives that are short on area
            // Expand the cutoff for objectives that haven't met their target
            bool anyAdjusted = false;
            for (size_t obj = 0; obj < nObj; ++obj) {
                int64_t deficit = targetAreas[obj] - allocated[obj];
                if (deficit > 0) {
                    // Expand cutoff to try to claim more pixels
                    cutoff[obj] += deficit;
                    if (cutoff[obj] > total) cutoff[obj] = total;
                    anyAdjusted = true;
                }
            }

            if (!anyAdjusted) break;

            // If all met but there were conflicts, we may have over-allocated
            // some objectives. Trim from objectives that exceeded their target.
            if (allMet) {
                // Trim excess allocations by removing worst-ranked pixels
                for (size_t obj = 0; obj < nObj; ++obj) {
                    if (allocated[obj] > targetAreas[obj]) {
                        int64_t excess = allocated[obj] - targetAreas[obj];
                        // Walk from worst rank up and deallocate
                        for (int64_t r = cutoff[obj] - 1; r >= 0 && excess > 0; --r) {
                            int64_t pixIdx = sortedIndices[obj][r];
                            if (allocation[pixIdx] == static_cast<int>(obj) + 1) {
                                allocation[pixIdx] = 0;
                                allocated[obj]--;
                                --excess;
                            }
                        }
                    }
                }
                break;
            }

            if (iteration % 10 == 0) {
                reportProgress(0.2 + 0.7 * std::min(1.0, static_cast<double>(iteration) / MAX_ITERATIONS),
                    QString("Iteration %1, conflicts: %2").arg(iteration).arg(conflictCount));
            }
        }

        reportProgress(0.95, "Writing output...");

        // Create output raster: 0 = unallocated, 1..nObj = objective assignment
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(suitMaps[0]->geoTransform());
        output.setProjection(suitMaps[0]->projection());
        output.setNoDataValue(0);

        auto& out = output.data(0);
        for (int64_t i = 0; i < total; ++i) {
            out[i] = static_cast<double>(allocation[i]);
        }

        bool writeOk = GdalIO::write(output, outputPath);

        reportProgress(1.0,
            QString("Allocation complete after %1 iterations, %2 objectives")
                .arg(iteration).arg(nObj));

        return writeOk;
    }

private:
    std::vector<QString> parsePaths(const QString& str) const {
        std::vector<QString> result;
        if (str.isEmpty()) return result;
        QStringList parts = str.split(",", Qt::SkipEmptyParts);
        for (const auto& p : parts) {
            QString trimmed = p.trimmed();
            if (!trimmed.isEmpty())
                result.push_back(trimmed);
        }
        return result;
    }

    std::vector<int64_t> parseInts(const QString& str) const {
        std::vector<int64_t> result;
        if (str.isEmpty()) return result;
        QStringList parts = str.split(",", Qt::SkipEmptyParts);
        for (const auto& p : parts) {
            bool ok;
            int64_t val = p.trimmed().toLongLong(&ok);
            if (ok) result.push_back(val);
        }
        return result;
    }
};

REGISTER_MODULE(MolaModule)

} // namespace aplaceholder
