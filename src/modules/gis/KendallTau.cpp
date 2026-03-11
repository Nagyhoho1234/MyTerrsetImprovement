#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace aplaceholder {

class KendallTauModule : public Module {
public:
    QString name() const override { return "KENDALLTAU"; }
    QString description() const override {
        return "Kendall Tau-b rank correlation between two rasters. "
               "For each pixel pair, counts concordant and discordant pairs "
               "and outputs a tau value raster.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input1", "First input raster image"),
            ParameterDef::file("input2", "Second input raster image"),
            ParameterDef::output("output", "Output tau value raster"),
        };
    }

    bool execute() override {
        auto r1 = GdalIO::read(parameter("input1").toString());
        auto r2 = GdalIO::read(parameter("input2").toString());
        if (!r1 || !r2) {
            setError("Failed to read input rasters");
            return false;
        }

        if (r1->cols() != r2->cols() || r1->rows() != r2->rows()) {
            setError("Input rasters must have the same dimensions");
            return false;
        }

        int cols = r1->cols(), rows = r1->rows();
        int bands1 = r1->bands();
        int bands2 = r2->bands();
        int numBands = std::min(bands1, bands2);

        if (numBands < 2) {
            setError("Input rasters must have at least 2 bands for pairwise tau computation");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData1 = r1->noDataValue();
        double noData2 = r2->noDataValue();
        bool hasND1 = r1->hasNoData();
        bool hasND2 = r2->hasNoData();

        reportProgress(0.1, "Computing Kendall Tau-b per pixel...");

        Raster result(cols, rows, 1, DataType::Float64);
        result.setGeoTransform(r1->geoTransform());
        result.setProjection(r1->projection());
        result.setNoDataValue(noData1);
        auto& out = result.data(0);

        for (int64_t px = 0; px < total; ++px) {
            // Collect band values for this pixel from both rasters
            std::vector<double> x, y;
            bool isNoData = false;

            for (int b = 0; b < numBands; ++b) {
                double v1 = r1->data(b)[px];
                double v2 = r2->data(b)[px];

                if ((hasND1 && v1 == noData1) || (hasND2 && v2 == noData2)) {
                    isNoData = true;
                    break;
                }
                x.push_back(v1);
                y.push_back(v2);
            }

            if (isNoData || x.size() < 2) {
                out[px] = noData1;
                continue;
            }

            int64_t N = static_cast<int64_t>(x.size());
            int64_t concordant = 0, discordant = 0;
            int64_t tiedX = 0, tiedY = 0;

            for (int64_t i = 0; i < N - 1; ++i) {
                for (int64_t j = i + 1; j < N; ++j) {
                    double dx = x[j] - x[i];
                    double dy = y[j] - y[i];

                    if (dx == 0.0 && dy == 0.0) {
                        // tied on both — not counted
                    } else if (dx == 0.0) {
                        ++tiedX;
                    } else if (dy == 0.0) {
                        ++tiedY;
                    } else if ((dx > 0 && dy > 0) || (dx < 0 && dy < 0)) {
                        ++concordant;
                    } else {
                        ++discordant;
                    }
                }
            }

            double n1 = static_cast<double>(concordant + discordant + tiedX);
            double n2 = static_cast<double>(concordant + discordant + tiedY);
            double denom = std::sqrt(n1 * n2);

            if (denom > 0.0) {
                out[px] = static_cast<double>(concordant - discordant) / denom;
            } else {
                out[px] = 0.0;
            }

            if (px % 100000 == 0)
                reportProgress(0.1 + 0.8 * static_cast<double>(px) / total);
        }

        reportProgress(0.95, "Writing output raster...");

        if (!GdalIO::write(result, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0, "Kendall Tau-b computation complete.");
        return true;
    }
};

REGISTER_MODULE(KendallTauModule)

} // namespace aplaceholder
