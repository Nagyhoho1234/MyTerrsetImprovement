#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>

namespace aplaceholder {

class StChoiceModule : public Module {
public:
    QString name() const override { return "STCHOICE"; }
    QString description() const override {
        return "Spatial-temporal choice model for land cover change allocation. "
               "Takes a transition probability matrix and a land cover map, "
               "allocates change spatially using stochastic selection weighted "
               "by transition potentials.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("landcover", "Input land cover raster"),
            ParameterDef::file("transition_matrix", "Transition probability matrix file"),
            ParameterDef::file("potential_raster", "Transition potential raster"),
            ParameterDef::output("output", "Output allocated land cover raster"),
            ParameterDef::integer("num_classes", "Number of land cover classes", 5, 2, 256,
                "Number of distinct land cover categories"),
            ParameterDef::integer("seed", "Random seed", 0, 0, 999999,
                "Seed for stochastic selection (0 = random)"),
        };
    }

    bool execute() override {
        auto landcover = GdalIO::read(parameter("landcover").toString());
        if (!landcover) { setError("Failed to read land cover raster"); return false; }

        auto potential = GdalIO::read(parameter("potential_raster").toString());
        if (!potential) { setError("Failed to read transition potential raster"); return false; }

        int cols = landcover->cols(), rows = landcover->rows();
        int numClasses = parameter("num_classes").toInt();
        int seed = parameter("seed").toInt();

        if (potential->cols() != cols || potential->rows() != rows) {
            setError("Potential raster dimensions must match land cover raster");
            return false;
        }

        const auto& lcData = landcover->data(0);
        const auto& potData = potential->data(0);
        double noData = landcover->noDataValue();
        bool hasND = landcover->hasNoData();
        int64_t total = static_cast<int64_t>(cols) * rows;

        // Parse transition probability matrix file
        // Format: each line is "from_class to_class probability"
        QString matrixPath = parameter("transition_matrix").toString();
        if (matrixPath.isEmpty()) {
            setError("Transition probability matrix file is required");
            return false;
        }

        // Initialize transition matrix with zeros
        std::vector<std::vector<double>> transMatrix(numClasses,
            std::vector<double>(numClasses, 0.0));

        std::ifstream fin(matrixPath.toStdString());
        if (!fin.is_open()) {
            setError("Failed to open transition probability matrix file: " + matrixPath);
            return false;
        }

        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            int fromClass, toClass;
            double prob;
            if (sscanf(line.c_str(), "%d %d %lf", &fromClass, &toClass, &prob) == 3) {
                if (fromClass >= 0 && fromClass < numClasses &&
                    toClass >= 0 && toClass < numClasses) {
                    transMatrix[fromClass][toClass] = prob;
                }
            }
        }
        fin.close();

        // Validate rows sum to ~1.0
        for (int i = 0; i < numClasses; ++i) {
            double rowSum = 0.0;
            for (int j = 0; j < numClasses; ++j) {
                rowSum += transMatrix[i][j];
            }
            if (rowSum > 0.0 && std::abs(rowSum - 1.0) > 0.05) {
                setError("Transition probabilities for class " +
                    QString::number(i) + " sum to " +
                    QString::number(rowSum) + " (expected ~1.0)");
                return false;
            }
        }

        // Set up random number generator
        std::mt19937 rng;
        if (seed > 0) {
            rng.seed(static_cast<unsigned>(seed));
        } else {
            std::random_device rd;
            rng.seed(rd());
        }
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(landcover->geoTransform());
        output.setProjection(landcover->projection());
        output.setNoDataValue(noData);
        auto& out = output.data(0);

        reportProgress(0.0, "Allocating land cover changes...");

        for (int64_t i = 0; i < total; ++i) {
            if (hasND && lcData[i] == noData) {
                out[i] = noData;
                continue;
            }

            int currentClass = static_cast<int>(lcData[i]);
            if (currentClass < 0 || currentClass >= numClasses) {
                out[i] = lcData[i]; // Keep original if out of range
                continue;
            }

            // Get transition potential weight for this pixel
            double potWeight = potData[i];
            if (potWeight <= 0.0) {
                out[i] = lcData[i]; // No change potential
                continue;
            }

            // Build cumulative probability distribution weighted by potential
            std::vector<double> cumProb(numClasses, 0.0);
            double cumSum = 0.0;
            for (int j = 0; j < numClasses; ++j) {
                double adjustedProb = transMatrix[currentClass][j];
                if (j != currentClass) {
                    // Weight non-self transitions by the spatial potential
                    adjustedProb *= potWeight;
                }
                cumSum += adjustedProb;
                cumProb[j] = cumSum;
            }

            // Normalize cumulative probabilities
            if (cumSum > 0.0) {
                for (int j = 0; j < numClasses; ++j) {
                    cumProb[j] /= cumSum;
                }
            }

            // Stochastic selection
            double rand = dist(rng);
            int newClass = currentClass; // Default: no change
            for (int j = 0; j < numClasses; ++j) {
                if (rand <= cumProb[j]) {
                    newClass = j;
                    break;
                }
            }

            out[i] = static_cast<double>(newClass);

            if (i % 1000000 == 0)
                reportProgress(static_cast<double>(i) / total);
        }

        reportProgress(1.0, "Writing output...");
        return GdalIO::write(output, parameter("output").toString());
    }
};

REGISTER_MODULE(StChoiceModule)

} // namespace aplaceholder
