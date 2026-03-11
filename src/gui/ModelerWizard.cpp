#include "ModelerWizard.h"
#include "ModelerTab.h"
#include <QTabWidget>
#include <QProgressBar>
#include <QLabel>
#include <QVBoxLayout>
#include <QDialogButtonBox>
#include <QMessageBox>

namespace aplaceholder {

ModelerWizard::ModelerWizard(const QString& title, QWidget* parent)
    : QDialog(parent)
{
    setWindowTitle(title);
    resize(900, 650);

    auto* layout = new QVBoxLayout(this);

    m_tabWidget = new QTabWidget;
    layout->addWidget(m_tabWidget, 1);

    // Progress area
    m_statusLabel = new QLabel("Ready");
    layout->addWidget(m_statusLabel);

    m_progressBar = new QProgressBar;
    m_progressBar->setRange(0, 1000);
    m_progressBar->setValue(0);
    m_progressBar->setTextVisible(true);
    layout->addWidget(m_progressBar);

    // Close button
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    layout->addWidget(buttons);

    connect(m_tabWidget, &QTabWidget::currentChanged,
            this, &ModelerWizard::onTabChangeRequested);
}

void ModelerWizard::addTab(ModelerTab* tab)
{
    m_tabWidget->addTab(tab, tab->tabName());

    connect(tab, &ModelerTab::progressUpdated,
            this, &ModelerWizard::onProgressUpdated);
    connect(tab, &ModelerTab::statusMessage,
            m_statusLabel, &QLabel::setText);
}

void ModelerWizard::onTabChangeRequested(int index)
{
    // Enforce sequential gating: cannot skip ahead past an incomplete tab
    for (int i = 0; i < index; ++i) {
        auto* tab = qobject_cast<ModelerTab*>(m_tabWidget->widget(i));
        if (tab && !tab->isComplete()) {
            QMessageBox::information(this, "Step Required",
                QString("Please complete the \"%1\" step first.")
                    .arg(tab->tabName()));
            m_tabWidget->setCurrentIndex(i);
            return;
        }
    }

    // Notify the new tab it is active
    auto* tab = qobject_cast<ModelerTab*>(m_tabWidget->widget(index));
    if (tab)
        tab->onActivated();
}

void ModelerWizard::onProgressUpdated(double fraction, const QString& msg)
{
    m_progressBar->setValue(static_cast<int>(fraction * 1000));
    if (!msg.isEmpty())
        m_statusLabel->setText(msg);
}

} // namespace aplaceholder
