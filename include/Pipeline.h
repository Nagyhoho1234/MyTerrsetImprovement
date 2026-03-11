#pragma once

#include <QString>
#include <QMap>
#include <QVariant>
#include <vector>
#include <memory>

namespace aplaceholder {

class Module;

// Pipeline — chain modules together (like TerrSet's Macro Modeler)
class Pipeline {
public:
    struct Step {
        QString moduleName;
        QMap<QString, QVariant> params;
    };

    void addStep(const QString& moduleName, const QMap<QString, QVariant>& params = {});
    void clear();
    int stepCount() const { return static_cast<int>(m_steps.size()); }

    // Execute all steps sequentially
    bool execute();

    // Save/load pipeline to JSON
    bool saveToFile(const QString& path) const;
    static Pipeline loadFromFile(const QString& path);

    QString lastError() const { return m_lastError; }

private:
    std::vector<Step> m_steps;
    QString m_lastError;
};

} // namespace aplaceholder
