#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

class MceModule : public Module {
public:
    QString name() const override { return "MCE"; }
    QString description() const override {
        return "Multi-criteria evaluation using Boolean intersection, "
               "Weighted Linear Combination (WLC), or Ordered Weighted Average (OWA). "
               "Combines factor and constraint images into a suitability surface.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            {"factors", "Factor images (comma-separated paths)", ParameterDef::String,
             {}, {}, 0, 0, "Comma-separated paths to factor images (0-255 byte range)", true},
            {"constraints", "Constraint images (comma-separated paths)", ParameterDef::String,
             {}, {}, 0, 0, "Comma-separated paths to Boolean constraint images", false},
            {"weights", "Factor weights (comma-separated)", ParameterDef::String,
             {}, {}, 0, 0, "Weights for each factor, must sum to 1.0", false},
            ParameterDef::combo("method", "Evaluation method",
                {"Boolean", "WLC", "OWA"}, 1,
                "Boolean=AND intersection, WLC=Weighted Linear Combination, OWA=Ordered Weighted Average"),
            {"order_weights", "Order weights for OWA (comma-separated)", ParameterDef::String,
             {}, {}, 0, 0, "Order weights for OWA method, must sum to 1.0", false},
            ParameterDef::output("output", "Output suitability image"),
        };
    }

    bool execute() override {
        int method = parameter("method").toInt();
        QString outputPath = parameter("output").toString();

        // Parse factor file paths
        std::vector<QString> factorPaths = parsePaths(parameter("factors").toString());
        std::vector<QString> constraintPaths = parsePaths(parameter("constraints").toString());

        if (factorPaths.empty() && constraintPaths.empty()) {
            setError("At least one factor or constraint image is required");
            return false;
        }

        // Parse weights
        std::vector<double> weights = parseDoubles(parameter("weights").toString());
        std::vector<double> orderWeights = parseDoubles(parameter("order_weights").toString());

        // Load factor rasters
        std::vector<std::unique_ptr<Raster>> factors;
        for (const auto& path : factorPaths) {
            auto r = GdalIO::read(path.trimmed());
            if (!r) {
                setError("Failed to read factor image: " + path);
                return false;
            }
            factors.push_back(std::move(r));
        }

        // Load constraint rasters
        std::vector<std::unique_ptr<Raster>> constraints;
        for (const auto& path : constraintPaths) {
            auto r = GdalIO::read(path.trimmed());
            if (!r) {
                setError("Failed to read constraint image: " + path);
                return false;
            }
            constraints.push_back(std::move(r));
        }

        // Determine reference dimensions from the first available raster
        const Raster* ref = nullptr;
        if (!factors.empty()) ref = factors[0].get();
        else if (!constraints.empty()) ref = constraints[0].get();

        int cols = ref->cols(), rows = ref->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Validate that all rasters have matching dimensions
        for (size_t i = 0; i < factors.size(); ++i) {
            if (factors[i]->cols() != cols || factors[i]->rows() != rows) {
                setError("Factor image " + factorPaths[i] + " has mismatched dimensions");
                return false;
            }
        }
        for (size_t i = 0; i < constraints.size(); ++i) {
            if (constraints[i]->cols() != cols || constraints[i]->rows() != rows) {
                setError("Constraint image " + constraintPaths[i] + " has mismatched dimensions");
                return false;
            }
        }

        // Validate weights for WLC and OWA
        size_t nFactors = factors.size();
        if (method == 1 || method == 2) { // WLC or OWA
            if (weights.size() != nFactors) {
                setError("Number of weights (" + QString::number(weights.size()) +
                         ") must match number of factors (" + QString::number(nFactors) + ")");
                return false;
            }
            double wsum = std::accumulate(weights.begin(), weights.end(), 0.0);
            if (std::abs(wsum - 1.0) > 0.001) {
                setError("Factor weights must sum to 1.0 (got " + QString::number(wsum) + ")");
                return false;
            }
        }

        if (method == 2) { // OWA
            if (orderWeights.size() != nFactors) {
                setError("Number of order weights must match number of factors");
                return false;
            }
            double owsum = std::accumulate(orderWeights.begin(), orderWeights.end(), 0.0);
            if (std::abs(owsum - 1.0) > 0.001) {
                setError("Order weights must sum to 1.0 (got " + QString::number(owsum) + ")");
                return false;
            }
        }

        // Create output raster
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(ref->geoTransform());
        output.setProjection(ref->projection());
        output.setNoDataValue(-9999);
        auto& out = output.data(0);

        reportProgress(0.0, "Computing MCE...");

        for (int64_t i = 0; i < total; ++i) {
            // Compute constraint product (AND logic)
            double constraintProduct = 1.0;
            for (size_t c = 0; c < constraints.size(); ++c) {
                double cv = constraints[c]->data(0)[i];
                constraintProduct *= cv;
            }

            if (method == 0) {
                // Boolean: multiply all factor and constraint images together
                double val = constraintProduct;
                for (size_t f = 0; f < factors.size(); ++f) {
                    val *= factors[f]->data(0)[i];
                }
                out[i] = val;
            } else if (method == 1) {
                // WLC: S = sum(wi * xi) * product(cj)
                if (constraintProduct == 0.0) {
                    out[i] = 0.0;
                } else {
                    double wlcSum = 0.0;
                    for (size_t f = 0; f < nFactors; ++f) {
                        wlcSum += weights[f] * factors[f]->data(0)[i];
                    }
                    out[i] = wlcSum * constraintProduct;
                }
            } else if (method == 2) {
                // OWA: sort weighted factor scores, apply order weights
                if (constraintProduct == 0.0) {
                    out[i] = 0.0;
                } else {
                    // Step 1: compute weighted factor values (wi * xi)
                    std::vector<double> weightedValues(nFactors);
                    for (size_t f = 0; f < nFactors; ++f) {
                        weightedValues[f] = weights[f] * factors[f]->data(0)[i];
                    }

                    // Step 2: sort from lowest to highest
                    std::sort(weightedValues.begin(), weightedValues.end());

                    // Step 3: apply order weights
                    // Order weight 1 goes to lowest ranked, order weight n to highest
                    double owaSum = 0.0;
                    for (size_t k = 0; k < nFactors; ++k) {
                        owaSum += orderWeights[k] * weightedValues[k];
                    }

                    // Scale by number of factors to keep output in factor range
                    // OWA result = n * sum(ow_k * sorted_wv_k) * product(cj)
                    out[i] = static_cast<double>(nFactors) * owaSum * constraintProduct;
                }
            }

            if (i % 500000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, outputPath);
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

    std::vector<double> parseDoubles(const QString& str) const {
        std::vector<double> result;
        if (str.isEmpty()) return result;
        QStringList parts = str.split(",", Qt::SkipEmptyParts);
        for (const auto& p : parts) {
            bool ok;
            double val = p.trimmed().toDouble(&ok);
            if (ok) result.push_back(val);
        }
        return result;
    }
};

REGISTER_MODULE(MceModule)

} // namespace aplaceholder
