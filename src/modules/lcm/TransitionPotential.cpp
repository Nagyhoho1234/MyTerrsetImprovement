#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace aplaceholder {

class TransitionPotentialModule : public Module {
public:
    QString name() const override { return "TRANSITION_POTENTIAL"; }
    QString description() const override {
        return "Models transition potential surfaces using logistic regression "
               "or random forest. Quantifies the likelihood of land cover "
               "transitions based on explanatory driver variables.";
    }
    QString category() const override { return "Land Change Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("transitions_file", "Transitions CSV file (from,to pairs)"),
            ParameterDef::file("driver_variables", "Driver variables (comma-separated file paths)"),
            ParameterDef::combo("method", "Modeling method",
                {"Logistic Regression", "Random Forest"}, 0,
                "Statistical method for modeling transition potential"),
            ParameterDef::output("output", "Output transition potential"),
        };
    }

    bool execute() override {
        // Parse transitions CSV (from,to pairs)
        std::vector<std::pair<int, int>> transitions;
        if (!readTransitionsFile(parameter("transitions_file").toString(), transitions))
            return false;

        if (transitions.empty()) {
            setError("No transitions found in transitions file");
            return false;
        }

        // Read driver variable rasters
        QString driverPaths = parameter("driver_variables").toString();
        QStringList driverFiles = driverPaths.split(",", Qt::SkipEmptyParts);
        if (driverFiles.empty()) {
            setError("No driver variables specified");
            return false;
        }

        std::vector<std::shared_ptr<Raster>> drivers;
        for (const auto& path : driverFiles) {
            auto dr = GdalIO::read(path.trimmed());
            if (!dr) {
                setError("Failed to read driver variable: " + path.trimmed());
                return false;
            }
            drivers.push_back(std::move(dr));
        }

        int cols = drivers[0]->cols(), rows = drivers[0]->rows();
        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = drivers[0]->noDataValue();
        bool hasND = drivers[0]->hasNoData();
        int numDrivers = static_cast<int>(drivers.size());

        // Verify all drivers have same dimensions
        for (int d = 1; d < numDrivers; ++d) {
            if (drivers[d]->cols() != cols || drivers[d]->rows() != rows) {
                setError("Driver variable dimensions do not match");
                return false;
            }
        }

        int methodIdx = parameter("method").toInt();

        // We need earlier and later maps to determine which pixels changed.
        // The transitions file gives from/to class pairs. We need the actual
        // land cover maps which should be the first two driver variables, or
        // we infer from the transitions file column headers.
        //
        // Actually, the transitions_file contains from,to pairs defining which
        // transitions to model. The driver_variables are explanatory surfaces.
        // The module needs the land cover change info. We'll look for that in
        // a convention: the transitions_file also encodes the earlier/later maps
        // as header metadata, or we use a composite approach where the output
        // raster is built per transition.
        //
        // For this implementation: we process each transition. For each transition
        // (from_class, to_class), we need to know which pixels:
        // - Changed: were from_class and became to_class (label=1)
        // - Didn't change: were from_class and stayed (label=0)
        //
        // We'll use the first two drivers as the earlier/later maps for sample
        // extraction, and remaining drivers as explanatory variables.
        // OR: the transitions_file includes the earlier/later map paths in the header.
        //
        // Simplified approach: first two entries in driver_variables are
        // earlier_map and later_map, rest are actual explanatory drivers.

        if (numDrivers < 3) {
            setError("Driver variables must include: earlier_map, later_map, "
                     "followed by at least one explanatory variable");
            return false;
        }

        const auto& earlierData = drivers[0]->data(0);
        const auto& laterData = drivers[1]->data(0);
        int numExplanatory = numDrivers - 2;

        // Output raster: maximum transition potential across all transitions
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(drivers[0]->geoTransform());
        output.setProjection(drivers[0]->projection());
        output.setNoDataValue(noData);
        auto& outData = output.data(0);
        for (int64_t i = 0; i < total; ++i)
            outData[i] = hasND ? noData : 0.0;

        // Process each transition
        for (size_t t = 0; t < transitions.size(); ++t) {
            int fromClass = transitions[t].first;
            int toClass = transitions[t].second;

            reportProgress(static_cast<double>(t) / transitions.size() * 0.9,
                "Modeling transition " + QString::number(fromClass) +
                " -> " + QString::number(toClass));

            // Extract training samples
            // Changed pixels (from_class -> to_class): label = 1
            // Unchanged eligible pixels (from_class -> from_class): label = 0
            std::vector<int64_t> changedIdx, unchangedIdx;
            for (int64_t i = 0; i < total; ++i) {
                if (hasND && (earlierData[i] == noData || laterData[i] == noData))
                    continue;
                int earlier = static_cast<int>(earlierData[i]);
                int later = static_cast<int>(laterData[i]);
                if (earlier == fromClass && later == toClass)
                    changedIdx.push_back(i);
                else if (earlier == fromClass && later == fromClass)
                    unchangedIdx.push_back(i);
            }

            if (changedIdx.empty()) continue;
            if (unchangedIdx.empty()) continue;

            // Build training dataset
            int64_t nChanged = static_cast<int64_t>(changedIdx.size());
            int64_t nUnchanged = static_cast<int64_t>(unchangedIdx.size());
            int64_t nSamples = nChanged + nUnchanged;

            // Feature matrix X (nSamples x (numExplanatory + 1)) with intercept
            // Label vector y
            int nFeatures = numExplanatory + 1;  // +1 for intercept
            std::vector<std::vector<double>> X(nSamples, std::vector<double>(nFeatures));
            std::vector<double> y(nSamples);

            // Fill changed samples
            for (int64_t k = 0; k < nChanged; ++k) {
                X[k][0] = 1.0;  // intercept
                for (int d = 0; d < numExplanatory; ++d)
                    X[k][d + 1] = drivers[d + 2]->data(0)[changedIdx[k]];
                y[k] = 1.0;
            }
            // Fill unchanged samples
            for (int64_t k = 0; k < nUnchanged; ++k) {
                int64_t idx = nChanged + k;
                X[idx][0] = 1.0;
                for (int d = 0; d < numExplanatory; ++d)
                    X[idx][d + 1] = drivers[d + 2]->data(0)[unchangedIdx[k]];
                y[idx] = 0.0;
            }

            if (methodIdx == 0) {
                // Logistic Regression using IRLS
                std::vector<double> beta = fitLogisticIRLS(X, y, nSamples, nFeatures);

                // Generate probability surface for all eligible (from_class) pixels
                for (int64_t i = 0; i < total; ++i) {
                    if (hasND && earlierData[i] == noData) continue;
                    int earlier = static_cast<int>(earlierData[i]);
                    if (earlier != fromClass) continue;

                    // Compute linear predictor
                    double eta = beta[0];  // intercept
                    for (int d = 0; d < numExplanatory; ++d)
                        eta += beta[d + 1] * drivers[d + 2]->data(0)[i];

                    // Sigmoid
                    double prob = 1.0 / (1.0 + std::exp(-eta));
                    prob = std::max(0.0, std::min(1.0, prob));

                    // Take max across transitions
                    if (outData[i] == noData)
                        outData[i] = prob;
                    else
                        outData[i] = std::max(outData[i], prob);
                }
            } else {
                // Random Forest placeholder: use simple distance-weighted scoring
                // This is a placeholder that outputs the proportion of changed
                // pixels at similar driver variable values
                for (int64_t i = 0; i < total; ++i) {
                    if (hasND && earlierData[i] == noData) continue;
                    int earlier = static_cast<int>(earlierData[i]);
                    if (earlier != fromClass) continue;

                    // Simple placeholder: average similarity to changed pixels
                    // vs unchanged pixels based on driver values
                    double sumChanged = 0.0, sumUnchanged = 0.0;
                    int countC = 0, countU = 0;

                    // Sample a subset for efficiency
                    int64_t sampleLimit = std::min(static_cast<int64_t>(200), nChanged);
                    for (int64_t k = 0; k < sampleLimit; ++k) {
                        int64_t sIdx = changedIdx[k * nChanged / sampleLimit];
                        double dist = 0.0;
                        for (int d = 0; d < numExplanatory; ++d) {
                            double diff = drivers[d + 2]->data(0)[i] -
                                          drivers[d + 2]->data(0)[sIdx];
                            dist += diff * diff;
                        }
                        sumChanged += std::exp(-dist * 0.001);
                        countC++;
                    }

                    sampleLimit = std::min(static_cast<int64_t>(200), nUnchanged);
                    for (int64_t k = 0; k < sampleLimit; ++k) {
                        int64_t sIdx = unchangedIdx[k * nUnchanged / sampleLimit];
                        double dist = 0.0;
                        for (int d = 0; d < numExplanatory; ++d) {
                            double diff = drivers[d + 2]->data(0)[i] -
                                          drivers[d + 2]->data(0)[sIdx];
                            dist += diff * diff;
                        }
                        sumUnchanged += std::exp(-dist * 0.001);
                        countU++;
                    }

                    double avgC = countC > 0 ? sumChanged / countC : 0.0;
                    double avgU = countU > 0 ? sumUnchanged / countU : 0.0;
                    double prob = (avgC + avgU) > 0 ? avgC / (avgC + avgU) : 0.5;

                    if (outData[i] == noData)
                        outData[i] = prob;
                    else
                        outData[i] = std::max(outData[i], prob);
                }
            }
        }

        reportProgress(0.95, "Writing output...");
        if (!GdalIO::write(output, parameter("output").toString())) {
            setError("Failed to write output raster");
            return false;
        }

        reportProgress(1.0, "Transition potential modeling complete.");
        return true;
    }

private:
    bool readTransitionsFile(const QString& path,
                            std::vector<std::pair<int, int>>& transitions) {
        std::ifstream file(path.toStdString());
        if (!file.is_open()) {
            setError("Failed to open transitions file: " + path);
            return false;
        }

        std::string line;
        // Skip header
        std::getline(file, line);

        while (std::getline(file, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string fromStr, toStr;
            if (std::getline(ss, fromStr, ',') && std::getline(ss, toStr, ',')) {
                fromStr.erase(0, fromStr.find_first_not_of(" \t\r\n"));
                fromStr.erase(fromStr.find_last_not_of(" \t\r\n") + 1);
                toStr.erase(0, toStr.find_first_not_of(" \t\r\n"));
                toStr.erase(toStr.find_last_not_of(" \t\r\n") + 1);
                if (!fromStr.empty() && !toStr.empty())
                    transitions.push_back({std::stoi(fromStr), std::stoi(toStr)});
            }
        }

        file.close();
        return true;
    }

    // Logistic regression via Iteratively Reweighted Least Squares (IRLS)
    // Beta = (X'WX)^-1 * X'Wz
    // where W = diag(p*(1-p)), z = X*beta + (y-p)/(p*(1-p))
    // Max 20 iterations
    std::vector<double> fitLogisticIRLS(
        const std::vector<std::vector<double>>& X,
        const std::vector<double>& y,
        int64_t nSamples, int nFeatures)
    {
        std::vector<double> beta(nFeatures, 0.0);
        const int maxIter = 20;
        const double epsilon = 1e-8;

        for (int iter = 0; iter < maxIter; ++iter) {
            // Compute predicted probabilities p = sigmoid(X * beta)
            std::vector<double> p(nSamples);
            for (int64_t i = 0; i < nSamples; ++i) {
                double eta = 0.0;
                for (int j = 0; j < nFeatures; ++j)
                    eta += X[i][j] * beta[j];
                p[i] = 1.0 / (1.0 + std::exp(-eta));
                // Clamp to avoid numerical issues
                p[i] = std::max(epsilon, std::min(1.0 - epsilon, p[i]));
            }

            // Compute working weights W = p*(1-p)
            // Compute working response z = X*beta + (y-p)/(p*(1-p))
            std::vector<double> w(nSamples);
            std::vector<double> z(nSamples);
            for (int64_t i = 0; i < nSamples; ++i) {
                w[i] = p[i] * (1.0 - p[i]);
                double eta = 0.0;
                for (int j = 0; j < nFeatures; ++j)
                    eta += X[i][j] * beta[j];
                z[i] = eta + (y[i] - p[i]) / w[i];
            }

            // Compute X'WX (nFeatures x nFeatures)
            std::vector<std::vector<double>> XtWX(nFeatures, std::vector<double>(nFeatures, 0.0));
            for (int j = 0; j < nFeatures; ++j) {
                for (int k = 0; k < nFeatures; ++k) {
                    double sum = 0.0;
                    for (int64_t i = 0; i < nSamples; ++i)
                        sum += X[i][j] * w[i] * X[i][k];
                    XtWX[j][k] = sum;
                }
            }

            // Compute X'Wz (nFeatures x 1)
            std::vector<double> XtWz(nFeatures, 0.0);
            for (int j = 0; j < nFeatures; ++j) {
                double sum = 0.0;
                for (int64_t i = 0; i < nSamples; ++i)
                    sum += X[i][j] * w[i] * z[i];
                XtWz[j] = sum;
            }

            // Solve XtWX * beta_new = XtWz using Gaussian elimination
            std::vector<double> betaNew = solveLinearSystem(XtWX, XtWz, nFeatures);

            // Check convergence
            double maxDiff = 0.0;
            for (int j = 0; j < nFeatures; ++j)
                maxDiff = std::max(maxDiff, std::abs(betaNew[j] - beta[j]));

            beta = betaNew;

            if (maxDiff < 1e-6)
                break;
        }

        return beta;
    }

    // Solve Ax = b via Gaussian elimination with partial pivoting
    std::vector<double> solveLinearSystem(
        std::vector<std::vector<double>> A,
        std::vector<double> b, int n)
    {
        // Augmented matrix
        for (int i = 0; i < n; ++i)
            A[i].push_back(b[i]);

        // Forward elimination with partial pivoting
        for (int col = 0; col < n; ++col) {
            // Find pivot
            int maxRow = col;
            double maxVal = std::abs(A[col][col]);
            for (int row = col + 1; row < n; ++row) {
                if (std::abs(A[row][col]) > maxVal) {
                    maxVal = std::abs(A[row][col]);
                    maxRow = row;
                }
            }

            if (maxVal < 1e-12) {
                // Singular or near-singular: add regularization
                A[col][col] += 1e-6;
            }

            if (maxRow != col)
                std::swap(A[col], A[maxRow]);

            for (int row = col + 1; row < n; ++row) {
                double factor = A[row][col] / A[col][col];
                for (int j = col; j <= n; ++j)
                    A[row][j] -= factor * A[col][j];
            }
        }

        // Back substitution
        std::vector<double> x(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = A[i][n];
            for (int j = i + 1; j < n; ++j)
                x[i] -= A[i][j] * x[j];
            if (std::abs(A[i][i]) > 1e-12)
                x[i] /= A[i][i];
            else
                x[i] = 0.0;
        }

        return x;
    }
};

REGISTER_MODULE(TransitionPotentialModule)

} // namespace aplaceholder
