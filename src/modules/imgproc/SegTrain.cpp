#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <QFile>
#include <QTextStream>
#include <QStringList>

namespace aplaceholder {

class SegTrainModule : public Module {
public:
    QString name() const override { return "SEGTRAIN"; }
    QString description() const override {
        return "Generate training signatures from segments. Computes mean spectral "
               "signature per segment across input bands and outputs a signature file.";
    }
    QString category() const override { return "Image Processing"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("segment_raster", "Segment image (region IDs)",
                "Segmented raster from SEGMENT or SEGMENTATION"),
            ParameterDef::file("bands", "Band images (comma-separated)",
                "Comma-separated list of band raster file paths"),
            ParameterDef::output("output_sig_file", "Output signature file (.sig)"),
        };
    }

    bool execute() override {
        // Read segment raster
        auto segRaster = GdalIO::read(parameter("segment_raster").toString());
        if (!segRaster) {
            setError("Failed to read segment raster");
            return false;
        }

        int cols = segRaster->cols();
        int rows = segRaster->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double segND = segRaster->noDataValue();
        bool segHasND = segRaster->hasNoData();
        const auto& segData = segRaster->data(0);

        // Parse band file paths
        QStringList bandPaths = parameter("bands").toString().split(",", Qt::SkipEmptyParts);
        int numBands = bandPaths.size();
        if (numBands == 0) {
            setError("No band images specified");
            return false;
        }

        // Read all band rasters
        std::vector<std::unique_ptr<Raster>> bandRasters;
        for (int b = 0; b < numBands; ++b) {
            auto br = GdalIO::read(bandPaths[b].trimmed());
            if (!br) {
                setError("Failed to read band raster: " + bandPaths[b].trimmed());
                return false;
            }
            if (br->cols() != cols || br->rows() != rows) {
                setError("Band raster dimensions do not match segment raster: " +
                         bandPaths[b].trimmed());
                return false;
            }
            bandRasters.push_back(std::move(br));
        }

        reportProgress(0.2, "Computing spectral statistics per segment...");

        // Accumulate per-segment statistics
        // segStats[segID] = {sum per band, sum of squares per band, count}
        struct SegStats {
            std::vector<double> sum;
            std::vector<double> sumSq;
            int64_t count = 0;
        };
        std::unordered_map<int, SegStats> segStats;

        for (int64_t i = 0; i < total; ++i) {
            double sv = segData[i];
            if (segHasND && sv == segND) continue;

            int segID = static_cast<int>(sv);
            auto& stats = segStats[segID];
            if (stats.sum.empty()) {
                stats.sum.resize(numBands, 0.0);
                stats.sumSq.resize(numBands, 0.0);
            }

            bool skip = false;
            for (int b = 0; b < numBands; ++b) {
                double val = bandRasters[b]->data(0)[i];
                if (bandRasters[b]->hasNoData() && val == bandRasters[b]->noDataValue()) {
                    skip = true;
                    break;
                }
            }
            if (skip) continue;

            for (int b = 0; b < numBands; ++b) {
                double val = bandRasters[b]->data(0)[i];
                stats.sum[b] += val;
                stats.sumSq[b] += val * val;
            }
            stats.count++;

            if (i % (total / 20 + 1) == 0)
                reportProgress(0.2 + 0.5 * static_cast<double>(i) / total);
        }

        reportProgress(0.8, "Writing signature file...");

        // Write signature file
        QString outPath = parameter("output_sig_file").toString();
        QFile outFile(outPath);
        if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            setError("Failed to open output signature file: " + outPath);
            return false;
        }

        QTextStream out(&outFile);
        // Header
        out << "# Segment Training Signatures\n";
        out << "# Bands: " << numBands << "\n";
        out << "# Band files:";
        for (const auto& bp : bandPaths)
            out << " " << bp.trimmed();
        out << "\n";
        out << "# SegmentID\tCount";
        for (int b = 0; b < numBands; ++b)
            out << "\tMean_B" << (b + 1) << "\tStdDev_B" << (b + 1);
        out << "\n";

        // Write per-segment mean and stddev
        for (auto& [segID, stats] : segStats) {
            if (stats.count == 0) continue;

            out << segID << "\t" << stats.count;
            for (int b = 0; b < numBands; ++b) {
                double mean = stats.sum[b] / stats.count;
                double variance = (stats.sumSq[b] / stats.count) - (mean * mean);
                if (variance < 0.0) variance = 0.0;
                double stddev = std::sqrt(variance);
                out << "\t" << mean << "\t" << stddev;
            }
            out << "\n";
        }

        outFile.close();

        reportProgress(1.0, "Signature file written.");
        return true;
    }
};

REGISTER_MODULE(SegTrainModule)

} // namespace aplaceholder
