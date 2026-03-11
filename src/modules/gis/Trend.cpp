#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>

namespace aplaceholder {

struct TrendPoint {
    double x, y, value;
};

class TrendModule : public Module {
public:
    QString name() const override { return "TREND"; }
    QString description() const override {
        return "Polynomial Trend Surface Analysis. Fits a polynomial trend surface "
               "(orders 1-6) to point data and evaluates it over a raster grid. "
               "Reports R-squared and polynomial coefficients.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("points_file", "Input point data (CSV: x,y,value)",
                "CSV file with columns x, y, value"),
            ParameterDef::file("reference_raster", "Reference raster (extent/dimensions)",
                "Raster defining the output extent and resolution"),
            ParameterDef::output("output", "Output trend surface"),
            ParameterDef::integer("order", "Polynomial order", 2, 1, 6,
                "Polynomial order: 1=plane, 2=quadratic, 3=cubic, etc."),
            ParameterDef::output("output_report", "Output report file (optional)",
                "Text file with R-squared and coefficients"),
        };
    }

    bool execute() override {
        QString pointsPath = parameter("points_file").toString();
        QString refPath = parameter("reference_raster").toString();
        QString outPath = parameter("output").toString();
        int order = parameter("order").toInt();
        QString reportPath = parameter("output_report").toString();

        if (order < 1 || order > 6) {
            setError("Polynomial order must be between 1 and 6");
            return false;
        }

        // Read points from CSV
        std::vector<TrendPoint> points;
        {
            std::ifstream ifs(pointsPath.toStdString());
            if (!ifs.is_open()) {
                setError("Failed to open points file: " + pointsPath);
                return false;
            }
            std::string line;
            std::getline(ifs, line); // Skip header
            while (std::getline(ifs, line)) {
                if (line.empty()) continue;
                std::istringstream ss(line);
                TrendPoint pt;
                char sep;
                if (ss >> pt.x >> sep >> pt.y >> sep >> pt.value) {
                    points.push_back(pt);
                }
            }
        }

        if (points.empty()) {
            setError("No valid points read from CSV file");
            return false;
        }

        // Build list of term labels and count terms for given order
        // Terms for order N: all x^i * y^j where i+j <= N, ordered by total degree
        std::vector<std::pair<int,int>> terms; // (power of x, power of y)
        std::vector<QString> termLabels;
        for (int deg = 0; deg <= order; ++deg) {
            for (int i = deg; i >= 0; --i) {
                int j = deg - i;
                terms.push_back({i, j});
                if (i == 0 && j == 0)
                    termLabels.push_back("1");
                else if (j == 0)
                    termLabels.push_back(QString("x^%1").arg(i));
                else if (i == 0)
                    termLabels.push_back(QString("y^%1").arg(j));
                else
                    termLabels.push_back(QString("x^%1*y^%2").arg(i).arg(j));
            }
        }

        int nTerms = static_cast<int>(terms.size());
        int nPts = static_cast<int>(points.size());

        if (nPts < nTerms) {
            setError(QString("Not enough points (%1) for polynomial order %2 "
                             "(requires at least %3 points)")
                     .arg(nPts).arg(order).arg(nTerms));
            return false;
        }

        reportProgress(0.1, "Building design matrix...");

        // Build design matrix X (nPts x nTerms)
        std::vector<double> X(nPts * nTerms);
        std::vector<double> Y(nPts);

        for (int p = 0; p < nPts; ++p) {
            Y[p] = points[p].value;
            for (int t = 0; t < nTerms; ++t) {
                double val = 1.0;
                for (int k = 0; k < terms[t].first; ++k)
                    val *= points[p].x;
                for (int k = 0; k < terms[t].second; ++k)
                    val *= points[p].y;
                X[p * nTerms + t] = val;
            }
        }

        reportProgress(0.2, "Computing normal equations...");

        // Compute X'X (nTerms x nTerms)
        std::vector<double> XtX(nTerms * nTerms, 0.0);
        for (int i = 0; i < nTerms; ++i) {
            for (int j = i; j < nTerms; ++j) {
                double sum = 0.0;
                for (int p = 0; p < nPts; ++p)
                    sum += X[p * nTerms + i] * X[p * nTerms + j];
                XtX[i * nTerms + j] = sum;
                XtX[j * nTerms + i] = sum;
            }
        }

        // Compute X'Y (nTerms x 1)
        std::vector<double> XtY(nTerms, 0.0);
        for (int i = 0; i < nTerms; ++i) {
            double sum = 0.0;
            for (int p = 0; p < nPts; ++p)
                sum += X[p * nTerms + i] * Y[p];
            XtY[i] = sum;
        }

        reportProgress(0.3, "Solving via Cholesky decomposition...");

        // Cholesky decomposition: XtX = L * L^T
        std::vector<double> L(nTerms * nTerms, 0.0);
        for (int i = 0; i < nTerms; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = XtX[i * nTerms + j];
                for (int k = 0; k < j; ++k)
                    sum -= L[i * nTerms + k] * L[j * nTerms + k];
                if (i == j) {
                    if (sum <= 0.0) {
                        setError("Cholesky decomposition failed - matrix not positive definite. "
                                 "Try a lower polynomial order or check for collinear points.");
                        return false;
                    }
                    L[i * nTerms + j] = std::sqrt(sum);
                } else {
                    L[i * nTerms + j] = sum / L[j * nTerms + j];
                }
            }
        }

        // Forward substitution: L * w = XtY
        std::vector<double> w(nTerms);
        for (int i = 0; i < nTerms; ++i) {
            double sum = XtY[i];
            for (int j = 0; j < i; ++j)
                sum -= L[i * nTerms + j] * w[j];
            w[i] = sum / L[i * nTerms + i];
        }

        // Back substitution: L^T * b = w
        std::vector<double> coeffs(nTerms);
        for (int i = nTerms - 1; i >= 0; --i) {
            double sum = w[i];
            for (int j = i + 1; j < nTerms; ++j)
                sum -= L[j * nTerms + i] * coeffs[j];
            coeffs[i] = sum / L[i * nTerms + i];
        }

        reportProgress(0.4, "Computing R-squared...");

        // Compute R-squared
        double yMean = 0.0;
        for (int p = 0; p < nPts; ++p)
            yMean += Y[p];
        yMean /= nPts;

        double ssTot = 0.0, ssRes = 0.0;
        for (int p = 0; p < nPts; ++p) {
            double predicted = 0.0;
            for (int t = 0; t < nTerms; ++t)
                predicted += coeffs[t] * X[p * nTerms + t];
            double residual = Y[p] - predicted;
            ssRes += residual * residual;
            double diff = Y[p] - yMean;
            ssTot += diff * diff;
        }
        double rSquared = (ssTot > 0.0) ? (1.0 - ssRes / ssTot) : 0.0;

        // Read reference raster for extent and dimensions
        auto ref = GdalIO::readMetadata(refPath);
        if (!ref) {
            setError("Failed to read reference raster: " + refPath);
            return false;
        }

        reportProgress(0.5, "Evaluating polynomial on grid...");

        int cols = ref->cols();
        int rows = ref->rows();
        Raster output(cols, rows, 1, DataType::Float64);
        output.setGeoTransform(ref->geoTransform());
        output.setProjection(ref->projection());
        output.setNoDataValue(-9999.0);

        auto& outData = output.data(0);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                double px, py;
                output.colRowToXY(c, r, px, py);

                double val = 0.0;
                for (int t = 0; t < nTerms; ++t) {
                    double termVal = coeffs[t];
                    for (int k = 0; k < terms[t].first; ++k)
                        termVal *= px;
                    for (int k = 0; k < terms[t].second; ++k)
                        termVal *= py;
                    val += termVal;
                }
                outData[static_cast<size_t>(r) * cols + c] = val;
            }
            if (r % 50 == 0)
                reportProgress(0.5 + 0.4 * static_cast<double>(r) / rows);
        }

        reportProgress(0.9, "Writing output...");

        if (!GdalIO::write(output, outPath)) {
            setError("Failed to write output raster");
            return false;
        }

        // Write report if requested
        if (!reportPath.isEmpty()) {
            std::ofstream rpt(reportPath.toStdString());
            if (rpt.is_open()) {
                rpt << "Polynomial Trend Surface Analysis Report\n";
                rpt << "========================================\n\n";
                rpt << "Polynomial order: " << order << "\n";
                rpt << "Number of data points: " << nPts << "\n";
                rpt << "Number of terms: " << nTerms << "\n";
                rpt << "R-squared: " << rSquared << "\n\n";
                rpt << "Coefficients:\n";
                for (int t = 0; t < nTerms; ++t) {
                    rpt << "  " << termLabels[t].toStdString()
                        << " = " << coeffs[t] << "\n";
                }
                rpt << "\nResidual sum of squares: " << ssRes << "\n";
                rpt << "Total sum of squares: " << ssTot << "\n";
            }
        }

        reportProgress(1.0, "Complete");
        return true;
    }
};

REGISTER_MODULE(TrendModule)

} // namespace aplaceholder
