#include "ModelerTab.h"
#include "ModelerSession.h"
#include "Module.h"
#include "ModuleRegistry.h"

namespace aplaceholder {

ModelerTab::ModelerTab(ModelerSession* session, QWidget* parent)
    : QWidget(parent)
    , m_session(session)
{}

bool ModelerTab::runModule(const QString& moduleName,
                           const QMap<QString, QVariant>& params)
{
    auto module = ModuleRegistry::instance().create(moduleName);
    if (!module) {
        emit statusMessage(QString("Module '%1' is not available.").arg(moduleName));
        return false;
    }

    for (auto it = params.cbegin(); it != params.cend(); ++it)
        module->setParameter(it.key(), it.value());

    module->setProgressCallback([this](double frac, const QString& msg) {
        emit progressUpdated(frac, msg);
    });

    emit statusMessage(QString("Running %1...").arg(moduleName));

    if (!module->execute()) {
        emit statusMessage(QString("%1 failed: %2").arg(moduleName, module->lastError()));
        return false;
    }

    emit statusMessage(QString("%1 completed.").arg(moduleName));
    return true;
}

} // namespace aplaceholder
