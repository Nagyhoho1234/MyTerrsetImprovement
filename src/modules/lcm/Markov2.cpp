#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

namespace aplaceholder {

class Markov2Module : public Module {
public:
    QString name() const override { return "MARKOV2"; }
    QString description() const override {
        return "Enhanced Markov chain with multiple time steps. Computes transition probability "
               "matrix from two time periods, projects forward N steps, and outputs conditional "
               "probability surfaces.";
    }
    QString category() const override { return "Land Change Modeler"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::raster("time1_raster", "Time 1 land cover raster",
                "Earlier land cover classification raster"),
            ParameterDef::raster("time2_raster", "Time 2 land cover raster",
                "Later land cover classification raster"),
            ParameterDef::output("output_prefix", "Output prefix"),
            ParameterDef::integer("projection_steps", "Projection steps", 1, 1, 100,
                "Number of Markov chain steps to project forward"),
        };
    }

    bool execute() override {
        QString t1Path = parameter("time1_raster").toString();
        QString t2Path = parameter("time2_raster").toString();
        QString outPrefix = parameter("output_prefix").toString();
        int projSteps = parameter("projection_steps").toInt();

        reportProgress(0.0, "Reading input rasters...");

        auto t1Raster = GdalIO::read(t1Path);
        if (!t1Raster) { setError("Failed to read time 1 raster"); return false; }

        auto t2Raster = GdalIO::read(t2Path);
        if (!t2Raster) { setError("Failed to read time 2 raster"); return false; }

        int cols = t1Raster->cols(), rows = t1Raster->rows();
        if (t2Raster->cols() != cols || t2Raster->rows() != rows) {
            setError("Time 1 and Time 2 rasters must have the same dimensions");
            return false;
        }

        int64_t total = static_cast<int64_t>(cols) * rows;
        double noData = t1Raster->noDataValue();
        bool hasND = t1Raster->hasNoData();

        reportProgress(0.1, "Identifying land cover classes...");

        // Identify all unique classes
        std::set<int> classSet;
        for (int64_t i = 0; i < total; ++i) {
            double v1 = t1Raster->data(0)[i];
            double v2 = t2Raster->data(0)[i];
            if (hasND && (v1 == noData || v2 == noData)) continue;
            classSet.insert(static_cast<int>(v1));
            classSet.insert(static_cast<int>(v2));
        }

        std::vector<int> classes(classSet.begin(), classSet.end());
        std::sort(classes.begin(), classes.end());
        int nClasses = (int)classes.size();

        if (nClasses < 2) {
            setError("At least 2 land cover classes required");
            return false;
        }

        // Map class value to index
        std::map<int, int> classIndex;
        for (int c = 0; c < nClasses; ++c)
            classIndex[classes[c]] = c;

        reportProgress(0.2, "Computing transition matrix...");

        // Count transitions
        std::vector<std::vector<int64_t>> transCount(nClasses, std::vector<int64_t>(nClasses, 0));
        std::vector<int64_t> fromCount(nClasses, 0);

        for (int64_t i = 0; i < total; ++i) {
            double v1 = t1Raster->data(0)[i];
            double v2 = t2Raster->data(0)[i];
            if (hasND && (v1 == noData || v2 == noData)) continue;
            int from = classIndex[static_cast<int>(v1)];
            int to = classIndex[static_cast<int>(v2)];
            transCount[from][to]++;
            fromCount[from]++;
        }

        // Compute transition probability matrix
        std::vector<std::vector<double>> transMat(nClasses, std::vector<double>(nClasses, 0.0));
        for (int from = 0; from < nClasses; ++from) {
            if (fromCount[from] > 0) {
                for (int to = 0; to < nClasses; ++to)
                    transMat[from][to] = static_cast<double>(transCount[from][to]) / fromCount[from];
            } else {
                // No pixels in this class; identity row
                transMat[from][from] = 1.0;
            }
        }

        reportProgress(0.3, "Projecting transition matrix forward...");

        // Raise transition matrix to the power of projSteps via repeated matrix multiplication
        auto projMat = transMat;
        for (int step = 1; step < projSteps; ++step) {
            std::vector<std::vector<double>> result(nClasses, std::vector<double>(nClasses, 0.0));
            for (int i = 0; i < nClasses; ++i)
                for (int j = 0; j < nClasses; ++j)
                    for (int k = 0; k < nClasses; ++k)
                        result[i][j] += projMat[i][k] * transMat[k][j];
            projMat = result;

            reportProgress(0.3 + 0.2 * step / projSteps);
        }

        reportProgress(0.5, "Computing conditional probability surfaces...");

        // For each target class, output a conditional probability surface
        // P(class_j at time T+N | class at time 2)
        for (int targetClass = 0; targetClass < nClasses; ++targetClass) {
            Raster probOut(cols, rows, 1, DataType::Float64);
            probOut.setGeoTransform(t1Raster->geoTransform());
            probOut.setProjection(t1Raster->projection());
            probOut.setNoDataValue(noData);
            auto& probData = probOut.data(0);

            for (int64_t i = 0; i < total; ++i) {
                double v2 = t2Raster->data(0)[i];
                if (hasND && v2 == noData) {
                    probData[i] = noData;
                    continue;
                }
                int fromIdx = classIndex[static_cast<int>(v2)];
                probData[i] = projMat[fromIdx][targetClass];
            }

            QString path = QString("%1_prob_class_%2.tif").arg(outPrefix).arg(classes[targetClass]);
            if (!GdalIO::write(probOut, path)) {
                setError(QString("Failed to write: %1").arg(path));
                return false;
            }

            reportProgress(0.5 + 0.45 * (targetClass + 1.0) / nClasses);
        }

        reportProgress(1.0);
        return true;
    }
};

REGISTER_MODULE(Markov2Module)

} // namespace aplaceholder
