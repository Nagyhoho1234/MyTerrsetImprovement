#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <QStringList>

namespace aplaceholder {

class BeliefModule : public Module {
public:
    QString name() const override { return "BELIEF"; }
    QString description() const override {
        return "Dempster-Shafer evidence combination. Combines multiple belief/plausibility "
               "raster pairs using Dempster's rule of combination and outputs belief, "
               "plausibility, and belief interval maps.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef({
                "belief_rasters", "Belief rasters (comma-separated)",
                ParameterDef::String, {}, {}, 0, 0,
                "Comma-separated list of belief raster file paths", true
            }),
            ParameterDef({
                "plausibility_rasters", "Plausibility rasters (comma-separated)",
                ParameterDef::String, {}, {}, 0, 0,
                "Comma-separated list of plausibility raster file paths", true
            }),
            ParameterDef::output("output_belief", "Output belief image"),
            ParameterDef::output("output_plausibility", "Output plausibility image"),
            ParameterDef::output("output_interval", "Output belief interval image"),
        };
    }

    bool execute() override {
        QStringList belFiles = parameter("belief_rasters").toString().split(",",
            Qt::SkipEmptyParts);
        QStringList plFiles = parameter("plausibility_rasters").toString().split(",",
            Qt::SkipEmptyParts);

        if (belFiles.size() < 2 || plFiles.size() < 2) {
            setError("At least two belief/plausibility raster pairs are required");
            return false;
        }
        if (belFiles.size() != plFiles.size()) {
            setError("Number of belief rasters must match number of plausibility rasters");
            return false;
        }

        // Read all belief and plausibility rasters
        std::vector<std::unique_ptr<Raster>> belRasters, plRasters;
        for (const auto& f : belFiles) {
            auto r = GdalIO::read(f.trimmed());
            if (!r) { setError("Failed to read belief raster: " + f.trimmed()); return false; }
            belRasters.push_back(std::move(r));
        }
        for (const auto& f : plFiles) {
            auto r = GdalIO::read(f.trimmed());
            if (!r) { setError("Failed to read plausibility raster: " + f.trimmed()); return false; }
            plRasters.push_back(std::move(r));
        }

        int cols = belRasters[0]->cols(), rows = belRasters[0]->rows();
        for (size_t k = 1; k < belRasters.size(); ++k) {
            if (belRasters[k]->cols() != cols || belRasters[k]->rows() != rows ||
                plRasters[k]->cols() != cols || plRasters[k]->rows() != rows) {
                setError("All input rasters must have the same dimensions");
                return false;
            }
        }

        Raster outBel(cols, rows, 1, DataType::Float64);
        Raster outPl(cols, rows, 1, DataType::Float64);
        Raster outInt(cols, rows, 1, DataType::Float64);
        outBel.setGeoTransform(belRasters[0]->geoTransform());
        outBel.setProjection(belRasters[0]->projection());
        outBel.setNoDataValue(belRasters[0]->noDataValue());
        outPl.setGeoTransform(belRasters[0]->geoTransform());
        outPl.setProjection(belRasters[0]->projection());
        outPl.setNoDataValue(belRasters[0]->noDataValue());
        outInt.setGeoTransform(belRasters[0]->geoTransform());
        outInt.setProjection(belRasters[0]->projection());
        outInt.setNoDataValue(belRasters[0]->noDataValue());

        auto& dOutBel = outBel.data(0);
        auto& dOutPl = outPl.data(0);
        auto& dOutInt = outInt.data(0);
        double noData = belRasters[0]->noDataValue();
        bool hasND = belRasters[0]->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;
        int nPairs = static_cast<int>(belRasters.size());

        for (int64_t i = 0; i < total; ++i) {
            bool isNoData = false;
            if (hasND) {
                for (int k = 0; k < nPairs; ++k) {
                    if (belRasters[k]->data(0)[i] == noData ||
                        plRasters[k]->data(0)[i] == noData) {
                        isNoData = true;
                        break;
                    }
                }
            }
            if (isNoData) {
                dOutBel[i] = noData;
                dOutPl[i] = noData;
                dOutInt[i] = noData;
                continue;
            }

            // Dempster's rule: iterative pairwise combination
            // For each evidence source, BPA for hypothesis = belief,
            // BPA for ignorance (theta) = 1 - plausibility is disbelief,
            // ignorance = plausibility - belief.
            // Simplified two-hypothesis model: hypothesis H and complement ~H.
            // m(H) = belief, m(~H) = 1 - plausibility, m(theta) = plausibility - belief
            double bel = belRasters[0]->data(0)[i];
            double pl = plRasters[0]->data(0)[i];

            for (int k = 1; k < nPairs; ++k) {
                double bel2 = belRasters[k]->data(0)[i];
                double pl2 = plRasters[k]->data(0)[i];

                // BPAs for source 1: m1(H), m1(~H), m1(theta)
                double m1H = bel;
                double m1notH = 1.0 - pl;
                double m1T = pl - bel;

                // BPAs for source 2: m2(H), m2(~H), m2(theta)
                double m2H = bel2;
                double m2notH = 1.0 - pl2;
                double m2T = pl2 - bel2;

                // Conflict: intersections that yield empty set
                double conflict = m1H * m2notH + m1notH * m2H;

                if (conflict >= 1.0) {
                    dOutBel[i] = noData;
                    dOutPl[i] = noData;
                    dOutInt[i] = noData;
                    goto next_pixel;
                }

                double norm = 1.0 / (1.0 - conflict);

                // Combined BPAs
                double combH = (m1H * m2H + m1H * m2T + m1T * m2H) * norm;
                double combNotH = (m1notH * m2notH + m1notH * m2T + m1T * m2notH) * norm;
                // double combT = (m1T * m2T) * norm;

                bel = combH;
                pl = 1.0 - combNotH;
            }

            dOutBel[i] = bel;
            dOutPl[i] = pl;
            dOutInt[i] = pl - bel;

            next_pixel:
            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(0.9, "Writing outputs...");
        bool ok = GdalIO::write(outBel, parameter("output_belief").toString());
        ok = ok && GdalIO::write(outPl, parameter("output_plausibility").toString());
        ok = ok && GdalIO::write(outInt, parameter("output_interval").toString());

        if (!ok) {
            setError("Failed to write one or more output rasters");
            return false;
        }

        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(BeliefModule)

} // namespace aplaceholder
