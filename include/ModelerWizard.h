#pragma once

#include "ModelerSession.h"
#include <QDialog>

class QTabWidget;
class QProgressBar;
class QLabel;

namespace aplaceholder {

class ModelerTab;

/// Abstract base dialog for modeler wizards.
/// Provides a tabbed interface with progress bar, tab gating, and shared session.
class ModelerWizard : public QDialog {
    Q_OBJECT
public:
    explicit ModelerWizard(const QString& title, QWidget* parent = nullptr);

    void addTab(ModelerTab* tab);

protected:
    ModelerSession m_session;
    QTabWidget*    m_tabWidget   = nullptr;
    QProgressBar*  m_progressBar = nullptr;
    QLabel*        m_statusLabel = nullptr;

private slots:
    void onTabChangeRequested(int index);
    void onProgressUpdated(double fraction, const QString& msg);
};

} // namespace aplaceholder
