#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class KnnRegressModule : public Module {
public:
    QString name() const override { return "KNNREGRESS"; }
    QString description() const override {
        return "K-nearest neighbor regression. Takes a dependent raster, "
               "independent raster(s), K value, and sample points. For each "
               "prediction pixel, finds K nearest training samples and computes "
               "a distance-weighted average.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("dependent", "Dependent variable raster"),
            ParameterDef::file("independent", "Independent variable raster(s)"),
            ParameterDef::file("sample_points", "Sample points raster (nonzero = training site)"),
            ParameterDef::integer("k", "Number of nearest neighbors", 5),
            ParameterDef::output("output", "Output prediction raster"),
        };
    }

    bool execute() override {
        auto depRaster = GdalIO::read(parameter("dependent").toString());
        auto sampleRaster = GdalIO::read(parameter("sample_points").toString());
        if (!depRaster || !sampleRaster) {
            setError("Failed to read dependent or sample points raster");
            return false;
        }

        int cols = depRaster->cols(), rows = depRaster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        if (sampleRaster->cols() != cols || sampleRaster->rows() != rows) {
            setError("Sample points raster must have the same dimensions as dependent raster");
            return false;
        }

        // Read independent rasters
        QStringList indepPaths = parameter("independent").toStringList();
        std::vector<std::unique_ptr<Raster>> indepRasters;
        for (const auto& path : indepPaths) {
            auto r = GdalIO::read(path);
            if (!r) {
                setError("Failed to read independent raster: " + path);
                return false;
            }
            if (r->cols() != cols || r->rows() != rows) {
                setError("Independent raster dimensions mismatch: " + path);
                return false;
            }
            indepRasters.push_back(std::move(r));
        }

        int numIndep = static_cast<int>(indepRasters.size());
        if (numIndep == 0) {
            setError("At least one independent raster is required");
            return false;
        }

        int K = parameter("k").toInt();
        if (K < 1) {
            setError("K must be at least 1");
            return false;
        }

        double noDataDep = depRaster->noDataValue();
        bool hasNDDep = depRaster->hasNoData();
        const auto& depData = depRaster->data(0);
        const auto& sampleData = sampleRaster->data(0);
        double noDataSample = sampleRaster->noDataValue();
        bool hasNDSample = sampleRaster->hasNoData();

        reportProgress(0.05, "Collecting training samples...");

        // Collect training samples: positions where sample raster is nonzero
        struct TrainingSample {
            std::vector<double> features;  // independent variable values
            double dependent;               // dependent variable value
        };
        std::vector<TrainingSample> samples;

        for (int64_t i = 0; i < total; ++i) {
            if (hasNDSample && sampleData[i] == noDataSample) continue;
            if (sampleData[i] == 0.0) continue;
            if (hasNDDep && depData[i] == noDataDep) continue;

            TrainingSample s;
            s.dependent = depData[i];
            bool valid = true;
            for (int f = 0; f < numIndep; ++f) {
                const auto& fdata = indepRasters[f]->data(0);
                double noDataF = indepRasters[f]->noDataValue();
                bool hasNDF = indepRasters[f]->hasNoData();
                if (hasNDF && fdata[i] == noDataF) {
                    valid = false;
                    break;
                }
                s.features.push_back(fdata[i]);
            }
            if (valid) {
                samples.push_back(std::move(s));
            }
        }

        if (static_cast<int>(samples.size()) < K) {
            setError(QString("Insufficient training samples (%1) for K=%2")
                         .arg(samples.size()).arg(K));
            return false;
        }

        reportProgress(0.15, QString("Found %1 training samples. Predicting...").arg(samples.size()));

        // Create output raster
        Raster result(cols, rows, 1, DataType::Float64);
        result.setGeoTransform(depRaster->geoTransform());
        result.setProjection(depRaster->projection());
        result.setNoDataValue(noDataDep);
        auto& out = result.data(0);

        // For each pixel, compute distance in feature space to all training samples,
        // find K nearest, compute inverse-distance-weighted average
        for (int64_t px = 0; px < total; ++px) {
            // Build feature vector for this pixel
            std::vector<double> pixFeatures(numIndep);
            bool valid = true;
            for (int f = 0; f < numIndep; ++f) {
                const auto& fdata = indepRasters[f]->data(0);
                double noDataF = indepRasters[f]->noDataValue();
                bool hasNDF = indepRasters[f]->hasNoData();
                if (hasNDF && fdata[px] == noDataF) {
                    valid = false;
                    break;
                }
                pixFeatures[f] = fdata[px];
            }

            if (!valid) {
                out[px] = noDataDep;
                continue;
            }

            // Compute distances to all training samples
            std::vector<std::pair<double, int>> dists;
            dists.reserve(samples.size());
            for (int s = 0; s < static_cast<int>(samples.size()); ++s) {
                double dist2 = 0.0;
                for (int f = 0; f < numIndep; ++f) {
                    double diff = pixFeatures[f] - samples[s].features[f];
                    dist2 += diff * diff;
                }
                dists.emplace_back(dist2, s);
            }

            // Partial sort to find K nearest
            std::partial_sort(dists.begin(), dists.begin() + K, dists.end());

            // Inverse-distance weighted average
            double weightSum = 0.0;
            double valueSum = 0.0;

            for (int k = 0; k < K; ++k) {
                double dist = std::sqrt(dists[k].first);
                int idx = dists[k].second;

                if (dist < 1e-12) {
                    // Exact match - use this value directly
                    valueSum = samples[idx].dependent;
                    weightSum = 1.0;
                    break;
                }

                double w = 1.0 / dist;
                weightSum += w;
                valueSum += w * samples[idx].dependent;
            }

            out[px] = (weightSum > 0.0) ? valueSum / weightSum : noDataDep;

            if (px % 100000 == 0)
                reportProgress(0.15 + 0.8 * static_cast<double>(px) / total);
        }

        reportProgress(0.95, "Writing output raster...");

        if (!GdalIO::write(result, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0, QString("KNN regression complete (K=%1, %2 training samples).")
                                .arg(K).arg(samples.size()));
        return true;
    }
};

REGISTER_MODULE(KnnRegressModule)

} // namespace aplaceholder
