#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace aplaceholder {

class BelCalcModule : public Module {
public:
    QString name() const override { return "BELCALC"; }
    QString description() const override {
        return "Dempster-Shafer belief calculation utility. Combines multiple evidence source "
               "rasters (mass functions) using Dempster's rule of combination to produce a "
               "combined belief image.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("inputs", "Input mass function rasters (comma-separated)",
                "Comma-separated list of raster file paths, each representing a mass function / "
                "evidence source (values 0-1)"),
            ParameterDef::output("output", "Output combined belief image",
                "Output raster with combined Dempster-Shafer belief values"),
        };
    }

    bool execute() override {
        // --------------------------------------------------------------------
        // 1. Read input rasters (mass function sources)
        // --------------------------------------------------------------------
        QString inputsParam = parameter("inputs").toString();
        QStringList inputPaths = inputsParam.split(",", Qt::SkipEmptyParts);
        for (auto& p : inputPaths) p = p.trimmed();

        if (inputPaths.size() < 2) {
            setError("At least 2 input mass function rasters required for combination");
            return false;
        }

        std::vector<std::unique_ptr<Raster>> inputRasters;
        for (const auto& path : inputPaths) {
            auto r = GdalIO::read(path);
            if (!r) {
                setError("Failed to read input raster: " + path);
                return false;
            }
            inputRasters.push_back(std::move(r));
        }

        int numSources = static_cast<int>(inputRasters.size());
        int cols = inputRasters[0]->cols();
        int rows = inputRasters[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool hasND = inputRasters[0]->hasNoData();
        double noData = inputRasters[0]->noDataValue();

        for (int s = 1; s < numSources; ++s) {
            if (inputRasters[s]->cols() != cols || inputRasters[s]->rows() != rows) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
        }

        reportProgress(0.05, "Combining evidence using Dempster's rule...");

        // --------------------------------------------------------------------
        // 2. Collect data pointers
        // --------------------------------------------------------------------
        std::vector<const std::vector<double>*> sources(numSources);
        for (int s = 0; s < numSources; ++s)
            sources[s] = &inputRasters[s]->data(0);

        // --------------------------------------------------------------------
        // 3. Create output raster
        // --------------------------------------------------------------------
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(inputRasters[0]->geoTransform());
        output.setProjection(inputRasters[0]->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);

        // --------------------------------------------------------------------
        // 4. Apply Dempster's rule of combination pixel-by-pixel
        //
        // For the binary hypothesis case (H vs. not-H):
        //   Each source i provides:
        //     m_i(H) = belief in hypothesis (input value)
        //     m_i(not-H) = 1 - m_i(H) - m_i(theta)
        //     m_i(theta) = uncertainty = assigned as complement
        //
        // Simplified Dempster's rule for multiple independent sources:
        //   For two sources m1, m2 with mass on focal element A:
        //     Combined mass = (sum of products of consistent masses) / (1 - K)
        //     where K = sum of products of conflicting masses
        //
        // For simple support functions (mass on {H} and {Omega}):
        //   m1(H)=a, m1(Omega)=1-a
        //   m2(H)=b, m2(Omega)=1-b
        //   Combined: m(H) = (a + b - a*b) after normalization
        //   Belief(H) = m(H) (for simple support, belief = mass)
        //
        // Sequential combination for multiple sources:
        // --------------------------------------------------------------------

        for (int64_t idx = 0; idx < total; ++idx) {
            if (hasND && (*sources[0])[idx] == noData) {
                outData[idx] = noData;
                continue;
            }

            // Start with first source as initial combined belief
            double combinedBel = (*sources[0])[idx];
            // Clamp to [0,1]
            combinedBel = std::max(0.0, std::min(1.0, combinedBel));

            // Sequentially combine with each additional source using Dempster's rule
            // For simple support functions:
            //   m_combined(H) = 1 - (1-m1(H)) * (1-m2(H))  [when no conflict on not-H]
            //
            // More generally (allowing conflict):
            //   m1(H)=a, m1(Omega)=1-a
            //   m2(H)=b, m2(Omega)=1-b
            //   Conflict K = 0 (since both support same hypothesis)
            //   m_combined(H) = a*b + a*(1-b) + (1-a)*b = a + b - a*b
            //   m_combined(Omega) = (1-a)*(1-b)
            //
            // This is equivalent to: combined = 1 - product(1 - m_i) for all i

            for (int s = 1; s < numSources; ++s) {
                double m2 = (*sources[s])[idx];
                m2 = std::max(0.0, std::min(1.0, m2));

                // Dempster's combination for simple support functions
                double m1 = combinedBel;

                // Conflict mass
                // For the binary case with simple support:
                // No conflict between support for H and uncertainty
                // K = 0 when both mass functions are simple support for same hypothesis
                double K = 0.0;

                // Combined mass on H
                // m(H) = (m1_H * m2_H + m1_H * m2_Omega + m1_Omega * m2_H) / (1 - K)
                double m1_H = m1;
                double m1_Omega = 1.0 - m1;
                double m2_H = m2;
                double m2_Omega = 1.0 - m2;

                double combined_H = m1_H * m2_H + m1_H * m2_Omega + m1_Omega * m2_H;
                double combined_Omega = m1_Omega * m2_Omega;

                // Normalize (1 - K = combined_H + combined_Omega when K=0)
                double normFactor = combined_H + combined_Omega;
                if (normFactor > 1e-12) {
                    combinedBel = combined_H / normFactor;
                } else {
                    combinedBel = 0.0;
                }
            }

            outData[idx] = combinedBel;

            if (idx % 1000000 == 0)
                reportProgress(0.05 + 0.85 * static_cast<double>(idx) / total);
        }

        // --------------------------------------------------------------------
        // 5. Write output
        // --------------------------------------------------------------------
        reportProgress(0.9, "Writing output...");
        QString outputPath = parameter("output").toString();
        if (!GdalIO::write(output, outputPath)) {
            setError("Failed to write output: " + outputPath);
            return false;
        }

        reportProgress(1.0, "Done.");
        return true;
    }
};

REGISTER_MODULE(BelCalcModule)

} // namespace aplaceholder
