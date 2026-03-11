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

    // Progress callback (0.0 to 1.0)
    using ProgressCallback = std::function<void(double progress, const QString& message)>;
    void setProgressCallback(ProgressCallback cb) { m_progressCb = std::move(cb); }

    // Error handling
    QString lastError() const { return m_lastError; }

protected:
    void reportProgress(double progress, const QString& msg = {});
    void setError(const QString& error) { m_lastError = error; }

    QMap<QString, QVariant> m_params;
    ProgressCallback m_progressCb;
    QString m_lastError;
};

} // namespace aplaceholder
