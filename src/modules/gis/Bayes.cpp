#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class BayesModule : public Module {
public:
    QString name() const override { return "BAYES"; }
    QString description() const override {
        return "Bayesian probability update. Computes posterior probability "
               "from prior probability, likelihood, and evidence (marginal) rasters.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("prior", "Prior probability image",
                "A priori probability raster (0-1)"),
            ParameterDef::file("likelihood", "Likelihood image",
                "Conditional probability P(evidence|hypothesis) raster (0-1)"),
            ParameterDef::file("evidence", "Evidence (marginal) image",
                "Marginal probability P(evidence) raster (0-1)"),
            ParameterDef::output("output", "Output posterior probability image"),
        };
    }

    bool execute() override {
        auto prior = GdalIO::read(parameter("prior").toString());
        auto likelihood = GdalIO::read(parameter("likelihood").toString());
        auto evidence = GdalIO::read(parameter("evidence").toString());
        if (!prior || !likelihood || !evidence) {
            setError("Failed to read one or more input rasters");
            return false;
        }

        int cols = prior->cols(), rows = prior->rows();
        if (likelihood->cols() != cols || likelihood->rows() != rows ||
            evidence->cols() != cols || evidence->rows() != rows) {
            setError("All input rasters must have the same dimensions");
            return false;
        }

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(prior->geoTransform());
        output.setProjection(prior->projection());
        output.setNoDataValue(prior->noDataValue());

        const auto& dPrior = prior->data(0);
        const auto& dLikelihood = likelihood->data(0);
        const auto& dEvidence = evidence->data(0);
        auto& out = output.data(0);
        double noData = prior->noDataValue();
        bool hasND = prior->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (dPrior[i] == noData || dLikelihood[i] == noData ||
                          dEvidence[i] == noData)) {
                out[i] = noData;
                continue;
            }

            if (dEvidence[i] == 0.0) {
                out[i] = noData;
            } else {
                // posterior = (prior * likelihood) / evidence
                out[i] = (dPrior[i] * dLikelihood[i]) / dEvidence[i];
            }

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(BayesModule)

} // namespace aplaceholder
