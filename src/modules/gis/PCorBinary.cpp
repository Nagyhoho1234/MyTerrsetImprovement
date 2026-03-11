#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class PCorBinaryModule : public Module {
public:
    QString name() const override { return "PCORBINARY"; }
    QString description() const override {
        return "Prior-corrected binary classification. Adjusts posterior probabilities "
               "using Bayes' theorem with a known prior proportion. Corrects for sampling "
               "bias when the training sample prevalence differs from the true population "
               "prevalence. Input: probability map and prior proportion. "
               "Output: corrected probability map.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input_probability", "Input posterior probability raster (0-1)"),
            ParameterDef::file("prior_proportion", "True prior proportion of the positive class (0-1)"),
            ParameterDef::file("sample_prevalence", "Prevalence of the positive class in the training sample (0-1)"),
            ParameterDef::output("output", "Output corrected probability raster"),
        };
    }

    bool execute() override {
        auto rProb = GdalIO::read(parameter("input_probability").toString());
        if (!rProb) {
            setError("Failed to read input probability raster");
            return false;
        }

        bool ok1 = false, ok2 = false;
        double priorProp = parameter("prior_proportion").toString().toDouble(&ok1);
        double samplePrev = parameter("sample_prevalence").toString().toDouble(&ok2);

        if (!ok1 || priorProp <= 0.0 || priorProp >= 1.0) {
            setError("Prior proportion must be a value between 0 and 1 (exclusive)");
            return false;
        }

        if (!ok2 || samplePrev <= 0.0 || samplePrev >= 1.0) {
            setError("Sample prevalence must be a value between 0 and 1 (exclusive)");
            return false;
        }

        int cols = rProb->cols(), rows = rProb->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        const auto& probData = rProb->data(0);
        double noData = rProb->noDataValue();
        bool hasND = rProb->hasNoData();

        reportProgress(0.0, "Applying prior correction via Bayes' theorem...");

        // Bayes' correction formula:
        // P_corrected = (P_model * prior / sample_prev) /
        //               (P_model * prior / sample_prev + (1 - P_model) * (1 - prior) / (1 - sample_prev))
        //
        // This adjusts the model's posterior probability from sample prevalence to
        // the true prior proportion.

        double outNoData = -9999.0;

        Raster corrected(cols, rows, 1, DataType::Float64);
        corrected.setGeoTransform(rProb->geoTransform());
        corrected.setProjection(rProb->projection());
        corrected.setNoDataValue(outNoData);
        auto& outData = corrected.data(0);

        double ratioPos = priorProp / samplePrev;
        double ratioNeg = (1.0 - priorProp) / (1.0 - samplePrev);

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && probData[i] == noData) {
                outData[i] = outNoData;
                continue;
            }

            double p = probData[i];

            // Clamp to valid probability range
            if (p < 0.0) p = 0.0;
            if (p > 1.0) p = 1.0;

            double numerator = p * ratioPos;
            double denominator = numerator + (1.0 - p) * ratioNeg;

            if (denominator <= 0.0) {
                outData[i] = outNoData;
            } else {
                double correctedP = numerator / denominator;
                // Clamp result to [0, 1]
                if (correctedP < 0.0) correctedP = 0.0;
                if (correctedP > 1.0) correctedP = 1.0;
                outData[i] = correctedP;
            }
        }

        reportProgress(0.8, "Writing corrected probability raster...");

        if (!GdalIO::write(corrected, parameter("output").toString())) {
            setError("Failed to write corrected probability raster");
            return false;
        }

        reportProgress(1.0,
            QString("Prior correction applied: prior=%1, sample_prevalence=%2")
                .arg(priorProp, 0, 'f', 4)
                .arg(samplePrev, 0, 'f', 4));

        return true;
    }
};

REGISTER_MODULE(PCorBinaryModule)

} // namespace aplaceholder
