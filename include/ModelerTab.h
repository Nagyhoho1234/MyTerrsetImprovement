#pragma once

#include <QWidget>
#include <QString>
#include <QVariant>
#include <QMap>
#include <memory>

namespace aplaceholder {

class ModelerSession;
class Module;

/// Abstract base class for a single tab inside a modeler wizard.
class ModelerTab : public QWidget {
    Q_OBJECT
public:
    explicit ModelerTab(ModelerSession* session, QWidget* parent = nullptr);

    virtual QString tabName() const = 0;
    virtual QString tabDescription() const = 0;

    /// Called when this tab becomes the active tab — refresh UI from session.
    virtual void onActivated() {}

    /// Can the user proceed past this tab?
    virtual bool isComplete() const = 0;

signals:
    void progressUpdated(double fraction, const QString& message);
    void statusMessage(const QString& msg);
    void completionChanged();

protected:
    ModelerSession* m_session;

    /// Run a registered module by name with the given parameters.
    /// Progress is forwarded to progressUpdated(). Returns true on success.
    bool runModule(const QString& moduleName,
                   const QMap<QString, QVariant>& params);
};

} // namespace aplaceholder
