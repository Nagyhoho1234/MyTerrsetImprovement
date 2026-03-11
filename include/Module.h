#pragma once

#include <QString>
#include <QVariant>
#include <QMap>
#include <vector>
#include <memory>
#include <functional>
#include "Parameter.h"

namespace aplaceholder {

class Raster;

// Simple color for chart series (avoids QColor dependency in core)
struct ChartColor {
    int r = 31, g = 119, b = 180;  // default: blue
    ChartColor() = default;
    ChartColor(int r_, int g_, int b_) : r(r_), g(g_), b(b_) {}
};

// Chart result data — modules populate this for visual output
struct ChartSeries {
    QString label;
    std::vector<double> x;
    std::vector<double> y;
    ChartColor color;
};

struct ChartResult {
    enum Type { None, Histogram, Scatter, Line, Bar, Pie };
    Type type = None;
    QString title;
    QString xLabel;
    QString yLabel;
    std::vector<ChartSeries> series;
    std::vector<QString> categoryLabels;  // For bar/pie charts

    bool hasData() const { return type != None && !series.empty(); }
};

// Base class for all analytical modules
class Module {
public:
    virtual ~Module() = default;

    // Module identity
    virtual QString name() const = 0;
    virtual QString description() const = 0;
    virtual QString category() const = 0;

    // Parameters — defines what the module dialog shows
    virtual std::vector<ParameterDef> parameterDefs() const = 0;

    // Set parameters before execution
    void setParameter(const QString& key, const QVariant& value);
    QVariant parameter(const QString& key) const;

    // Execution
    virtual bool execute() = 0;

    // Chart result — modules set this during execute() for visual output
    const ChartResult& chartResult() const { return m_chartResult; }

    // Progress callback (0.0 to 1.0)
    using ProgressCallback = std::function<void(double progress, const QString& message)>;
    void setProgressCallback(ProgressCallback cb) { m_progressCb = std::move(cb); }

    // Error handling
    QString lastError() const { return m_lastError; }

protected:
    void reportProgress(double progress, const QString& msg = {});
    void setError(const QString& error) { m_lastError = error; }
    void setChartResult(ChartResult result) { m_chartResult = std::move(result); }

    QMap<QString, QVariant> m_params;
    ProgressCallback m_progressCb;
    QString m_lastError;
    ChartResult m_chartResult;
};

} // namespace aplaceholder
