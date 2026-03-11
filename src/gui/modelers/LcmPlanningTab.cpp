#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QTextEdit>
#include <QFileDialog>
#include <QMessageBox>

namespace aplaceholder {

class LcmPlanningTab : public ModelerTab {
    Q_OBJECT
public:
    explicit LcmPlanningTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Constraints and Incentives ---
        auto* constraintBox = new QGroupBox("Constraints and Incentives");
        auto* cLayout = new QGridLayout(constraintBox);

        cLayout->addWidget(new QLabel(
            "Specify constraint/incentive maps per transition.\n"
            "Values: 0 = absolute constraint, 1 = unconstrained, "
            "<1 = disincentive, >1 = incentive."), 0, 0, 1, 3);

        cLayout->addWidget(new QLabel("Constraint/incentive map:"), 1, 0);
        m_constraintEdit = new QLineEdit;
        cLayout->addWidget(m_constraintEdit, 1, 1);
        auto* browseConstraint = new QPushButton("Browse...");
        cLayout->addWidget(browseConstraint, 1, 2);

        layout->addWidget(constraintBox);

        // --- Planned Infrastructure Changes ---
        auto* infraBox = new QGroupBox("Planned Infrastructure Changes");
        auto* iLayout = new QGridLayout(infraBox);

        iLayout->addWidget(new QLabel(
            "Specify future infrastructure modifications and their effective dates.\n"
            "Infrastructure is integrated at each prediction stage when the stage date "
            "reaches the infrastructure date."), 0, 0, 1, 3);

        iLayout->addWidget(new QLabel("Infrastructure map:"), 1, 0);
        m_infraEdit = new QLineEdit;
        iLayout->addWidget(m_infraEdit, 1, 1);
        auto* browseInfra = new QPushButton("Browse...");
        iLayout->addWidget(browseInfra, 1, 2);

        iLayout->addWidget(new QLabel("Effective date (year):"), 2, 0);
        m_infraDateEdit = new QLineEdit;
        iLayout->addWidget(m_infraDateEdit, 2, 1);

        layout->addWidget(infraBox);

        // --- Apply ---
        m_applyBtn = new QPushButton("Apply Planning Parameters");
        layout->addWidget(m_applyBtn);

        // --- Notes ---
        auto* notesBox = new QGroupBox("Notes");
        auto* notesLayout = new QVBoxLayout(notesBox);
        m_notesEdit = new QTextEdit;
        m_notesEdit->setPlaceholderText(
            "Planning notes: document assumptions, scenarios, and interventions...");
        m_notesEdit->setMaximumHeight(100);
        notesLayout->addWidget(m_notesEdit);
        layout->addWidget(notesBox);

        layout->addStretch();

        // Connections
        connect(browseConstraint, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select Constraint Map",
                QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
            if (!path.isEmpty()) m_constraintEdit->setText(path);
        });
        connect(browseInfra, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select Infrastructure Map",
                QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
            if (!path.isEmpty()) m_infraEdit->setText(path);
        });
        connect(m_applyBtn, &QPushButton::clicked, this, [this]() {
            if (!m_constraintEdit->text().trimmed().isEmpty())
                m_session->set("lcm/constraint_map", m_constraintEdit->text().trimmed());
            if (!m_infraEdit->text().trimmed().isEmpty())
                m_session->set("lcm/infrastructure_map", m_infraEdit->text().trimmed());
            if (!m_infraDateEdit->text().trimmed().isEmpty())
                m_session->set("lcm/infrastructure_date", m_infraDateEdit->text().trimmed());
            emit statusMessage("Planning parameters applied.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "Planning"; }
    QString tabDescription() const override {
        return "Define constraints, incentives, and planned infrastructure changes.";
    }

    // Planning is optional — always complete
    bool isComplete() const override { return true; }

private:
    QLineEdit*   m_constraintEdit = nullptr;
    QLineEdit*   m_infraEdit      = nullptr;
    QLineEdit*   m_infraDateEdit  = nullptr;
    QPushButton* m_applyBtn       = nullptr;
    QTextEdit*   m_notesEdit      = nullptr;
};

ModelerTab* createLcmPlanningTab(ModelerSession* s, QWidget* p) {
    return new LcmPlanningTab(s, p);
}

} // namespace aplaceholder

#include "LcmPlanningTab.moc"
