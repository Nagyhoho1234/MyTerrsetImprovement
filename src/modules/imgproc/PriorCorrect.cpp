#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace aplaceholder {

class PriorCorrectModule : public Module {
public:
    QString name() const override { return "PRIORCORRECT"; }
    QString description() const override {
        return "Prior probability correction for soft classifiers. Adjusts posterior "
               "probability maps using Bayesian updating with user-supplied prior "
               "probability surfaces. Corrects for unequal class prevalence.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input posterior probability image (multi-band)",
                "Multi-band image where each band is the posterior probability for one class"),
            ParameterDef::file("priors", "Prior probability image (multi-band)",
                "Multi-band image where each band is the prior probability for one class"),
            ParameterDef::output("output", "Output corrected probability image"),
            ParameterDef::boolean("normalize", "Normalize output probabilities", true,
                "Ensure corrected probabilities sum to 1.0 per pixel"),
            ParameterDef::real("min_prior", "Minimum prior value", 0.001, 0.0, 1.0,
                "Floor value for priors to avoid division by zero"),
        };
    }

    bool execute() override {
        auto posteriors = GdalIO::read(parameter("input").toString());
        auto priors = GdalIO::read(parameter("priors").toString());
        if (!posteriors || !priors) {
            setError("Failed to read input posterior or prior probability image");
            return false;
        }

        int cols = posteriors->cols(), rows = posteriors->rows();
        int nClasses = posteriors->bands();
        int64_t total = static_cast<int64_t>(cols) * rows;
        bool doNormalize = parameter("normalize").toBool();
        double minPrior = parameter("min_prior").toDouble();

        if (priors->cols() != cols || priors->rows() != rows) {
            setError("Prior and posterior images must have the same dimensions");
            return false;
        }
        if (priors->bands() != nClasses) {
            setError("Prior and posterior images must have the same number of bands (classes)");
            return false;
        }

        bool hasND = posteriors->hasNoData();
        double noData = posteriors->noDataValue();

        std::vector<const std::vector<double>*> postBands(nClasses);
        std::vector<const std::vector<double>*> priorBands(nClasses);
        for (int c = 0; c < nClasses; ++c) {
            postBands[c] = &posteriors->data(c);
            priorBands[c] = &priors->data(c);
        }

        double trainingPrior = 1.0 / nClasses;

        reportProgress(0.10, "Applying Bayesian prior correction...");

        Raster output(cols, rows, nClasses, DataType::Float64);
        output.setGeoTransform(posteriors->geoTransform());
        output.setProjection(posteriors->projection());
        output.setNoDataValue(-9999.0);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && (*postBands[0])[i] == noData) {
                for (int c = 0; c < nClasses; ++c)
                    output.data(c)[i] = -9999.0;
                continue;
            }

            double sumCorrected = 0.0;
            std::vector<double> corrected(nClasses);

            for (int c = 0; c < nClasses; ++c) {
                double post = (*postBands[c])[i];
                double prior = (*priorBands[c])[i];
                if (prior < minPrior) prior = minPrior;
                corrected[c] = post * prior / trainingPrior;
                sumCorrected += corrected[c];
            }

            for (int c = 0; c < nClasses; ++c) {
                if (doNormalize && sumCorrected > 1e-15)
                    output.data(c)[i] = corrected[c] / sumCorrected;
                else
                    output.data(c)[i] = corrected[c];
            }

            if (i % 500000 == 0)
                reportProgress(0.10 + 0.85 * static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(PriorCorrectModule)

} // namespace aplaceholder
