#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>

namespace aplaceholder {

class SedimentModule : public Module {
public:
    QString name() const override { return "SEDIMENT_DELIVERY"; }
    QString description() const override {
        return "Sediment delivery ratio (SDR) and sediment export estimation. Computes "
               "SDR = exp(-beta * travel_time) where travel_time is derived from flow "
               "accumulation, then multiplies by RUSLE soil loss to estimate per-pixel "
               "sediment export (tons/pixel).";
    }
    QString category() const override { return "Ecosystem Services"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("rusle_raster", "RUSLE soil loss raster (tons/ha)",
                "Raster of potential soil loss from RUSLE (A = R*K*LS*C*P)"),
            ParameterDef::file("flow_accumulation", "Flow accumulation raster",
                "Flow accumulation raster used as proxy for travel time. "
                "Higher values indicate longer travel paths to stream."),
            ParameterDef::output("output", "Output sediment export raster",
                "Continuous raster of sediment export (RUSLE * SDR)"),
            ParameterDef::real("beta", "Beta coefficient", 1.0, 0.001, 100.0,
                "Exponential decay coefficient for SDR. Higher values mean "
                "sediment delivery decreases faster with travel time. Default 1.0."),
        };
    }

    bool execute() override {
        // 1. Load input rasters
        auto rusle = GdalIO::read(parameter("rusle_raster").toString());
        if (!rusle) {
            setError("Failed to read RUSLE raster.");
            return false;
        }

        auto flowAcc = GdalIO::read(parameter("flow_accumulation").toString());
        if (!flowAcc) {
            setError("Failed to read flow accumulation raster.");
            return false;
        }

        reportProgress(0.10, "Input rasters loaded.");

        // 2. Validate dimensions
        int cols = rusle->cols();
        int rows = rusle->rows();
        int64_t total = rusle->cellCount();

        if (flowAcc->cols() != cols || flowAcc->rows() != rows) {
            setError(QString("Dimension mismatch: RUSLE is %1x%2 but flow accumulation is %3x%4")
                         .arg(cols).arg(rows)
                         .arg(flowAcc->cols()).arg(flowAcc->rows()));
            return false;
        }

        // 3. Parameters
        double beta = parameter("beta").toDouble();

        // 4. NoData handling
        double noData = rusle->noDataValue();
        bool hasND_R = rusle->hasNoData();
        bool hasND_F = flowAcc->hasNoData();
        double ndR = rusle->noDataValue();
        double ndF = flowAcc->noDataValue();

        // 5. Normalize flow accumulation to get travel time proxy [0, 1]
        //    Higher flow accumulation = closer to stream = lower travel time
        const auto& fData = flowAcc->data(0);
        double maxFlowAcc = 0.0;
        for (int64_t i = 0; i < total; ++i) {
            if (hasND_F && fData[i] == ndF)
                continue;
            if (fData[i] > maxFlowAcc)
                maxFlowAcc = fData[i];
        }

        if (maxFlowAcc <= 0.0)
            maxFlowAcc = 1.0;

        reportProgress(0.20, QString("Max flow accumulation: %1").arg(maxFlowAcc));

        // 6. Compute SDR and sediment export
        //    travel_time_i = 1 - (flow_acc_i / max_flow_acc)
        //    SDR_i = exp(-beta * travel_time_i)
        //    sediment_export_i = RUSLE_i * SDR_i
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(rusle->geoTransform());
        output.setProjection(rusle->projection());
        output.setNoDataValue(noData);

        const auto& rData = rusle->data(0);
        auto& outData = output.data(0);

        int64_t validCount = 0;
        double totalSediment = 0.0;

        for (int64_t i = 0; i < total; ++i) {
            bool skip = false;
            if (hasND_R && rData[i] == ndR) skip = true;
            if (hasND_F && fData[i] == ndF) skip = true;

            if (skip) {
                outData[i] = noData;
                continue;
            }

            // Travel time: inverse of connectivity (high flow acc = close to stream = low travel time)
            double travelTime = 1.0 - (fData[i] / maxFlowAcc);

            // Sediment delivery ratio
            double sdr = std::exp(-beta * travelTime);

            // Sediment export
            outData[i] = rData[i] * sdr;
            totalSediment += outData[i];
            ++validCount;

            if (i % 500000 == 0)
                reportProgress(0.20 + 0.70 * static_cast<double>(i) / total);
        }

        reportProgress(0.90, "Writing output...");

        // 7. Write output
        QString outPath = parameter("output").toString();
        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster: " + outPath);
            return false;
        }

        reportProgress(1.0,
            QString("Done. Total sediment export: %1 tons, valid pixels: %2")
                .arg(totalSediment, 0, 'f', 2)
                .arg(validCount));

        return true;
    }
};

REGISTER_MODULE(SedimentModule)

} // namespace aplaceholder
