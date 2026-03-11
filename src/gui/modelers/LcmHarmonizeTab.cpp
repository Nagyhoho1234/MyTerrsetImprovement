#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QTextEdit>
#include <QMessageBox>
#include <QFile>

namespace aplaceholder {

class LcmHarmonizeTab : public ModelerTab {
    Q_OBJECT
public:
    explicit LcmHarmonizeTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Requirements ---
        auto* reqBox = new QGroupBox("Input Map Requirements");
        auto* reqLayout = new QVBoxLayout(reqBox);

        reqLayout->addWidget(new QLabel(
            "LCM requires that input land cover maps meet the following conditions:\n\n"
            "1. Legends in both maps must be the same\n"
            "2. Categories must be the same and sequential\n"
            "3. Backgrounds must be the same with a value of zero\n"
            "4. Spatial dimensions (resolution, coordinates) must be the same\n\n"
            "Use the tools below to fix any mismatches."));

        layout->addWidget(reqBox);

        // --- Validation Status ---
        auto* statusBox = new QGroupBox("Validation Status");
        auto* statusLayout = new QVBoxLayout(statusBox);

        m_statusText = new QTextEdit;
        m_statusText->setReadOnly(true);
        m_statusText->setMaximumHeight(150);
        m_statusText->setPlaceholderText("Click 'Validate Inputs' to check map compatibility.");
        statusLayout->addWidget(m_statusText);

        m_validateBtn = new QPushButton("Validate Inputs");
        statusLayout->addWidget(m_validateBtn);

        layout->addWidget(statusBox);

        // --- Harmonization Tools ---
        auto* toolsBox = new QGroupBox("Harmonization Tools");
        auto* toolsLayout = new QGridLayout(toolsBox);

        m_harmonizeLegendBtn = new QPushButton("Harmonize Legends");
        m_harmonizeLegendBtn->setToolTip("Ensure both maps use the same legend categories.");
        m_harmonizeLegendBtn->setEnabled(false);
        toolsLayout->addWidget(m_harmonizeLegendBtn, 0, 0);

        m_harmonizeBgBtn = new QPushButton("Harmonize Backgrounds");
        m_harmonizeBgBtn->setToolTip("Set background value to zero in both maps.");
        m_harmonizeBgBtn->setEnabled(false);
        toolsLayout->addWidget(m_harmonizeBgBtn, 0, 1);

        m_resampleBtn = new QPushButton("Resample to Match");
        m_resampleBtn->setToolTip("Resample one map to match the spatial extent and "
                                   "resolution of the other.");
        m_resampleBtn->setEnabled(false);
        toolsLayout->addWidget(m_resampleBtn, 1, 0);

        layout->addWidget(toolsBox);

        layout->addStretch();

        // Connections
        connect(m_validateBtn, &QPushButton::clicked, this, &LcmHarmonizeTab::validateInputs);
        connect(m_harmonizeLegendBtn, &QPushButton::clicked, this, [this]() {
            emit statusMessage("Legend harmonization: use RECLASS to remap categories.");
        });
        connect(m_harmonizeBgBtn, &QPushButton::clicked, this, [this]() {
            emit statusMessage("Background harmonization: use ASSIGN to set background to 0.");
        });
        connect(m_resampleBtn, &QPushButton::clicked, this, [this]() {
            emit statusMessage("Resampling: use RESAMPLE to match spatial dimensions.");
        });
    }

    QString tabName() const override { return "Harmonize"; }
    QString tabDescription() const override {
        return "Validate and harmonize input land cover maps for LCM compatibility.";
    }

    // Harmonize is optional/conditional
    bool isComplete() const override { return true; }

    void onActivated() override {
        validateInputs();
    }

private slots:
    void validateInputs() {
        if (!m_session->contains("lcm/earlier_map") || !m_session->contains("lcm/later_map")) {
            m_statusText->setText("No input maps specified yet. "
                                   "Set session parameters in the Change Analysis tab first.");
            return;
        }

        // Attempt to read both maps and compare metadata
        QString earlier = m_session->get("lcm/earlier_map").toString();
        QString later   = m_session->get("lcm/later_map").toString();

        QString report;
        report += "Earlier map: " + earlier + "\n";
        report += "Later map: " + later + "\n\n";

        // Basic file existence check
        bool ok = true;
        if (!QFile::exists(earlier)) {
            report += "ERROR: Earlier map file not found.\n";
            ok = false;
        }
        if (!QFile::exists(later)) {
            report += "ERROR: Later map file not found.\n";
            ok = false;
        }

        if (ok) {
            report += "Both input files exist.\n";
            report += "Detailed spatial/legend validation requires loading the rasters.\n";
            report += "Run Change Analysis to verify full compatibility.";
        }

        m_statusText->setText(report);

        bool needsFix = !ok;
        m_harmonizeLegendBtn->setEnabled(needsFix);
        m_harmonizeBgBtn->setEnabled(needsFix);
        m_resampleBtn->setEnabled(needsFix);
    }

private:
    QTextEdit*   m_statusText         = nullptr;
    QPushButton* m_validateBtn        = nullptr;
    QPushButton* m_harmonizeLegendBtn = nullptr;
    QPushButton* m_harmonizeBgBtn     = nullptr;
    QPushButton* m_resampleBtn        = nullptr;
};

ModelerTab* createLcmHarmonizeTab(ModelerSession* s, QWidget* p) {
    return new LcmHarmonizeTab(s, p);
}

} // namespace aplaceholder

#include "LcmHarmonizeTab.moc"
