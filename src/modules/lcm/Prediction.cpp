#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class PredictionModule : public Module {
public:
    QString name() const override { return "PREDICTION"; }
    QString description() const override {
        return "Generates hard or soft prediction maps from transition potentials "
               "and Markov chain transition probabilities. Hard predictions assign "
               "the most likely class; soft predictions output continuous probabilities.";
    }
    QString category() const override { return "Land Change Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("transition_potentials", "Transition potential surfaces (comma-separated)"),
            ParameterDef::file("markov_matrix", "Markov transition probability matrix (CSV)"),
            ParameterDef::output("output", "Output prediction map"),
            ParameterDef::combo("prediction_type", "Prediction type",
                {"Hard", "Soft"}, 0,
                "Hard assigns most likely class; soft outputs probabilities"),
        };
    }

    bool execute() override {
        // Read transition probability matrix
        std::vector<int> classes;
        std::vector<std::vector<double>> transMatrix;
        if (!readTransitionMatrix(parameter("markov_matrix").toString(),
                                  classes, transMatrix)) {
            return false;
        }
        int n = static_cast<int>(classes.size());

        std::map<int, int> classIndex;
        for (int i = 0; i < n; ++i)
            classIndex[classes[i]] = i;

        // Read transition potential rasters (comma-separated)
        QString potPaths = parameter("transition_potentials").toString();
        QStringList potFiles = potPaths.split(",", Qt::SkipEmptyParts);

        std::vector<std::shared_ptr<Raster>> potentials;
        for (const auto& path : potFiles) {
            auto pr = GdalIO::read(path.trimmed());
            if (!pr) {
                setError("Failed to read transition potential: " + path.trimmed());
                return false;
            }
            potentials.push_back(std::move(pr));
        }

        if (potentials.empty()) {
            setError("No transition potential surfaces provided");
            return false;
        }

        int cols = potentials[0]->cols(), rows = potentials[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = potentials[0]->noDataValue();
        bool hasND = potentials[0]->hasNoData();

        // Verify dimensions match
        for (size_t k = 1; k < potentials.size(); ++k) {
            if (potentials[k]->cols() != cols || potentials[k]->rows() != rows) {
                setError("Transition potential raster dimensions do not match");
                return false;
            }
        }

        // If fewer potential maps than classes, pad with zeros
        while (static_cast<int>(potentials.size()) < n) {
            auto zero = std::make_shared<Raster>(cols, rows, 1, DataType::Float64);
            auto& zd = zero->data(0);
            for (int64_t i = 0; i < total; ++i)
                zd[i] = 0.0;
            potentials.push_back(zero);
        }

        int predType = parameter("prediction_type").toInt();

        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(potentials[0]->geoTransform());
        output.setProjection(potentials[0]->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);

        reportProgress(0.0, predType == 0 ? "Computing hard prediction..." :
                                            "Computing soft prediction...");

        if (predType == 0) {
            // Hard prediction: MOLA-style competitive allocation
            // For each pixel, compute combined score for each target class:
            //   score[class_j] = sum_i(transition_probability[i][j] * potential[j])
            // Assign to class with highest score.
            for (int64_t i = 0; i < total; ++i) {
                if (hasND) {
                    bool anyNoData = false;
                    for (int k = 0; k < n; ++k) {
                        if (potentials[k]->data(0)[i] == noData) {
                            anyNoData = true;
                            break;
                        }
                    }
                    if (anyNoData) {
                        outData[i] = noData;
                        continue;
                    }
                }

                double bestScore = -1.0;
                int bestClass = classes[0];

                for (int j = 0; j < n; ++j) {
                    double potential = potentials[j]->data(0)[i];

                    // Combine with Markov probabilities: use average transition
                    // probability to this class weighted by transition potential
                    double avgProb = 0.0;
                    for (int k = 0; k < n; ++k)
                        avgProb += transMatrix[k][j];
                    avgProb /= n;

                    double score = potential * avgProb;

                    if (score > bestScore) {
                        bestScore = score;
                        bestClass = classes[j];
                    }
                }

                outData[i] = static_cast<double>(bestClass);

                if (i % 1000000 == 0)
                    reportProgress(0.9 * static_cast<double>(i) / total);
            }
        } else {
            // Soft prediction: maximum transition potential across all transitions
            for (int64_t i = 0; i < total; ++i) {
                if (hasND) {
                    bool anyNoData = false;
                    for (int k = 0; k < n; ++k) {
                        if (potentials[k]->data(0)[i] == noData) {
                            anyNoData = true;
                            break;
                        }
                    }
                    if (anyNoData) {
                        outData[i] = noData;
                        continue;
                    }
                }

                double maxPotential = 0.0;
                for (int j = 0; j < n; ++j) {
                    double val = potentials[j]->data(0)[i];
                    maxPotential = std::max(maxPotential, val);
                }

                outData[i] = maxPotential;

                if (i % 1000000 == 0)
                    reportProgress(0.9 * static_cast<double>(i) / total);
            }
        }

        reportProgress(0.95, "Writing output...");
        if (!GdalIO::write(output, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0, "Prediction complete.");
        return true;
    }

private:
    bool readTransitionMatrix(const QString& path,
                              std::vector<int>& classes,
                              std::vector<std::vector<double>>& matrix) {
        std::ifstream file(path.toStdString());
        if (!file.is_open()) {
            setError("Failed to open transition matrix file: " + path);
            return false;
        }

        std::string line;
        // Read header to get class labels
        if (!std::getline(file, line)) {
            setError("Empty transition matrix file");
            return false;
        }

        std::stringstream headerSS(line);
        std::string token;
        // Skip first column header (From\To)
        std::getline(headerSS, token, ',');
        while (std::getline(headerSS, token, ',')) {
            token.erase(0, token.find_first_not_of(" \t\r\n"));
            token.erase(token.find_last_not_of(" \t\r\n") + 1);
            if (!token.empty())
                classes.push_back(std::stoi(token));
        }

        int n = static_cast<int>(classes.size());
        matrix.resize(n, std::vector<double>(n, 0.0));

        int row = 0;
        while (std::getline(file, line) && row < n) {
            std::stringstream ss(line);
            // Skip row label
            std::getline(ss, token, ',');
            for (int col = 0; col < n && std::getline(ss, token, ','); ++col) {
                token.erase(0, token.find_first_not_of(" \t\r\n"));
                token.erase(token.find_last_not_of(" \t\r\n") + 1);
                if (!token.empty())
                    matrix[row][col] = std::stod(token);
            }
            ++row;
        }

        file.close();
        return true;
    }
};

REGISTER_MODULE(PredictionModule)

} // namespace aplaceholder
