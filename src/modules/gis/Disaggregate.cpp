#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <string>

namespace aplaceholder {

class DisaggregateModule : public Module {
public:
    QString name() const override { return "DISAGGREGATE"; }
    QString description() const override {
        return "Spatial disaggregation of coarse-resolution raster data. "
               "Takes a coarse-resolution raster and a fine-resolution template, "
               "disaggregates values proportionally across finer cells.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input coarse-resolution raster"),
            ParameterDef::file("template_raster", "Fine-resolution template raster"),
            ParameterDef::output("output", "Output disaggregated raster"),
            ParameterDef::combo("method", "Disaggregation method",
                {"Proportional", "Equal Distribution"}, 0,
                "Method for distributing coarse values to fine cells"),
        };
    }

    bool execute() override {
        auto input = GdalIO::read(parameter("input").toString());
        if (!input) { setError("Failed to read input coarse-resolution raster"); return false; }

        auto templateRst = GdalIO::read(parameter("template_raster").toString());
        if (!templateRst) { setError("Failed to read template raster"); return false; }

        int fineCols = templateRst->cols(), fineRows = templateRst->rows();
        int coarseCols = input->cols(), coarseRows = input->rows();
        int method = parameter("method").toInt();

        const auto& coarseData = input->data(0);
        const auto& templateData = templateRst->data(0);
        double noData = input->noDataValue();
        bool hasND = input->hasNoData();
        double templateND = templateRst->noDataValue();
        bool templateHasND = templateRst->hasNoData();

        // Determine the ratio between coarse and fine resolution
        auto coarseGT = input->geoTransform();
        auto fineGT = templateRst->geoTransform();
        double xRatio = std::abs(coarseGT.pixelWidth / fineGT.pixelWidth);
        double yRatio = std::abs(coarseGT.pixelHeight / fineGT.pixelHeight);
        int cellsPerCoarseX = static_cast<int>(std::round(xRatio));
        int cellsPerCoarseY = static_cast<int>(std::round(yRatio));

        if (cellsPerCoarseX < 1 || cellsPerCoarseY < 1) {
            setError("Template resolution must be finer than input resolution");
            return false;
        }

        int64_t totalFine = static_cast<int64_t>(fineCols) * fineRows;

        Raster output(fineCols, fineRows, 1, DataType::Float64);
        output.setGeoTransform(fineGT);
        output.setProjection(templateRst->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        reportProgress(0.0, "Disaggregating...");

        if (method == 0) {
            // Proportional: distribute coarse value proportionally based on
            // template weights within each coarse cell
            for (int cr = 0; cr < coarseRows; ++cr) {
                for (int cc = 0; cc < coarseCols; ++cc) {
                    int64_t coarseIdx = static_cast<int64_t>(cr) * coarseCols + cc;
                    double coarseVal = coarseData[coarseIdx];

                    if (hasND && coarseVal == noData) {
                        // Fill corresponding fine cells with noData
                        for (int fy = 0; fy < cellsPerCoarseY; ++fy) {
                            for (int fx = 0; fx < cellsPerCoarseX; ++fx) {
                                int fineR = cr * cellsPerCoarseY + fy;
                                int fineC = cc * cellsPerCoarseX + fx;
                                if (fineR < fineRows && fineC < fineCols) {
                                    out[static_cast<int64_t>(fineR) * fineCols + fineC] = noData;
                                }
                            }
                        }
                        continue;
                    }

                    // Sum template weights within this coarse cell
                    double weightSum = 0.0;
                    std::vector<std::pair<int64_t, double>> fineCells;
                    for (int fy = 0; fy < cellsPerCoarseY; ++fy) {
                        for (int fx = 0; fx < cellsPerCoarseX; ++fx) {
                            int fineR = cr * cellsPerCoarseY + fy;
                            int fineC = cc * cellsPerCoarseX + fx;
                            if (fineR >= fineRows || fineC >= fineCols) continue;
                            int64_t fineIdx = static_cast<int64_t>(fineR) * fineCols + fineC;
                            double tw = templateData[fineIdx];
                            if (templateHasND && tw == templateND) {
                                fineCells.push_back({fineIdx, 0.0});
                            } else {
                                fineCells.push_back({fineIdx, tw});
                                weightSum += tw;
                            }
                        }
                    }

                    // Distribute proportionally
                    for (const auto& cell : fineCells) {
                        if (weightSum > 0.0 && cell.second > 0.0) {
                            out[cell.first] = coarseVal * (cell.second / weightSum);
                        } else if (weightSum <= 0.0) {
                            // Equal distribution fallback when all weights are zero
                            int validCount = static_cast<int>(fineCells.size());
                            out[cell.first] = (validCount > 0) ? coarseVal / validCount : noData;
                        } else {
                            out[cell.first] = noData;
                        }
                    }
                }

                if (cr % 100 == 0)
                    reportProgress(static_cast<double>(cr) / coarseRows);
            }
        } else {
            // Equal Distribution: spread coarse value equally across fine cells
            for (int cr = 0; cr < coarseRows; ++cr) {
                for (int cc = 0; cc < coarseCols; ++cc) {
                    int64_t coarseIdx = static_cast<int64_t>(cr) * coarseCols + cc;
                    double coarseVal = coarseData[coarseIdx];

                    int validCount = 0;
                    std::vector<int64_t> fineIndices;
                    for (int fy = 0; fy < cellsPerCoarseY; ++fy) {
                        for (int fx = 0; fx < cellsPerCoarseX; ++fx) {
                            int fineR = cr * cellsPerCoarseY + fy;
                            int fineC = cc * cellsPerCoarseX + fx;
                            if (fineR < fineRows && fineC < fineCols) {
                                int64_t fineIdx = static_cast<int64_t>(fineR) * fineCols + fineC;
                                fineIndices.push_back(fineIdx);
                                if (!(templateHasND && templateData[fineIdx] == templateND)) {
                                    validCount++;
                                }
                            }
                        }
                    }

                    if (hasND && coarseVal == noData) {
                        for (int64_t idx : fineIndices)
                            out[idx] = noData;
                        continue;
                    }

                    double share = (validCount > 0) ? coarseVal / validCount : noData;
                    for (int64_t idx : fineIndices) {
                        if (templateHasND && templateData[idx] == templateND) {
                            out[idx] = noData;
                        } else {
                            out[idx] = share;
                        }
                    }
                }

                if (cr % 100 == 0)
                    reportProgress(static_cast<double>(cr) / coarseRows);
            }
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(DisaggregateModule)

} // namespace aplaceholder
