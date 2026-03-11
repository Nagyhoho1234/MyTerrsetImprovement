#include "Pipeline.h"
#include "Module.h"
#include "ModuleRegistry.h"
#include <QJsonDocument>
#include <QJsonArray>
#include <QJsonObject>
#include <QFile>

namespace aplaceholder {

void Pipeline::addStep(const QString& moduleName, const QMap<QString, QVariant>& params)
{
    m_steps.push_back({moduleName, params});
}

void Pipeline::clear()
{
    m_steps.clear();
}

bool Pipeline::execute()
{
    for (size_t i = 0; i < m_steps.size(); ++i) {
        const auto& step = m_steps[i];
        auto module = ModuleRegistry::instance().create(step.moduleName);
        if (!module) {
            m_lastError = QString("Unknown module: %1").arg(step.moduleName);
            return false;
        }

        for (auto it = step.params.begin(); it != step.params.end(); ++it) {
            module->setParameter(it.key(), it.value());
        }

        if (!module->execute()) {
            m_lastError = QString("Step %1 (%2) failed: %3")
                .arg(i + 1).arg(step.moduleName).arg(module->lastError());
            return false;
        }
    }
    return true;
}

bool Pipeline::saveToFile(const QString& path) const
{
    QJsonArray steps;
    for (const auto& step : m_steps) {
        QJsonObject obj;
        obj["module"] = step.moduleName;
        QJsonObject params;
        for (auto it = step.params.begin(); it != step.params.end(); ++it) {
            params[it.key()] = QJsonValue::fromVariant(it.value());
        }
        obj["params"] = params;
        steps.append(obj);
    }

    QJsonDocument doc(steps);
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
        return false;
    }
    file.write(doc.toJson());
    return true;
}

Pipeline Pipeline::loadFromFile(const QString& path)
{
    Pipeline pipeline;
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
        return pipeline;

    auto doc = QJsonDocument::fromJson(file.readAll());
    auto steps = doc.array();
    for (const auto& val : steps) {
        auto obj = val.toObject();
        QMap<QString, QVariant> params;
        auto paramsObj = obj["params"].toObject();
        for (auto it = paramsObj.begin(); it != paramsObj.end(); ++it) {
            params[it.key()] = it.value().toVariant();
        }
        pipeline.addStep(obj["module"].toString(), params);
    }
    return pipeline;
}

} // namespace aplaceholder
