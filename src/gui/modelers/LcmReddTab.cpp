#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QTableWidget>
#include <QHeaderView>
#include <QFileDialog>
#include <QMessageBox>

namespace aplaceholder {

class LcmReddTab : public ModelerTab {
    Q_OBJECT
public:
    explicit LcmReddTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- CO2 Emissions ---
        auto* co2Box = new QGroupBox("Calculate CO2 Emissions");
        auto* co2Layout = new QGridLayout(co2Box);

        co2Layout->addWidget(new QLabel(
            "Specify carbon pool densities (tC/ha) as constant values or raster images."),
            0, 0, 1, 3);

        const QStringList pools = {
            "Above-ground biomass",
            "Below-ground biomass",
            "Dead wood",
            "Harvest wood products",
            "Litter",
            "Soil organic carbon"
        };

        m_poolTable = new QTableWidget(pools.size(), 3);
        m_poolTable->setHorizontalHeaderLabels({"Carbon Pool", "Type", "Value / Path"});
        m_poolTable->horizontalHeader()->setStretchLastSection(true);
        m_poolTable->setMaximumHeight(200);

        for (int i = 0; i < pools.size(); ++i) {
            m_poolTable->setItem(i, 0, new QTableWidgetItem(pools[i]));
            auto* typeCombo = new QComboBox;
            typeCombo->addItems({"Constant (tC/ha)", "Raster image"});
            m_poolTable->setCellWidget(i, 1, typeCombo);
            m_poolTable->setItem(i, 2, new QTableWidgetItem("0.0"));
        }

        co2Layout->addWidget(m_poolTable, 1, 0, 1, 3);
        layout->addWidget(co2Box);

        // --- Non-CO2 Emissions ---
        auto* nonCo2Box = new QGroupBox("Calculate Non-CO2 Emissions");
        auto* nonLayout = new QGridLayout(nonCo2Box);

        nonLayout->addWidget(new QLabel(
            "Estimate CH4 and N2O emissions from fire-based deforestation."), 0, 0, 1, 2);

        nonLayout->addWidget(new QLabel("Proportion of forest burned:"), 1, 0);
        m_burnFraction = new QDoubleSpinBox;
        m_burnFraction->setRange(0.0, 1.0);
        m_burnFraction->setSingleStep(0.05);
        m_burnFraction->setValue(0.5);
        nonLayout->addWidget(m_burnFraction, 1, 1);

        nonLayout->addWidget(new QLabel("Combustion efficiency:"), 2, 0);
        m_combustionEff = new QDoubleSpinBox;
        m_combustionEff->setRange(0.0, 1.0);
        m_combustionEff->setSingleStep(0.05);
        m_combustionEff->setValue(0.5);
        nonLayout->addWidget(m_combustionEff, 2, 1);

        layout->addWidget(nonCo2Box);

        // --- Net GHG Emissions ---
        auto* netBox = new QGroupBox("Calculate Net GHG Emissions");
        auto* netLayout = new QGridLayout(netBox);

        netLayout->addWidget(new QLabel(
            "Estimate net GHG emission reductions from REDD project.\n"
            "C-REDD = C-Baseline - C-Actual - C-Leakage"), 0, 0, 1, 2);

        netLayout->addWidget(new QLabel("Leakage rate:"), 1, 0);
        m_leakageRate = new QDoubleSpinBox;
        m_leakageRate->setRange(0.0, 1.0);
        m_leakageRate->setSingleStep(0.05);
        m_leakageRate->setValue(0.1);
        netLayout->addWidget(m_leakageRate, 1, 1);

        netLayout->addWidget(new QLabel("Success rate:"), 2, 0);
        m_successRate = new QDoubleSpinBox;
        m_successRate->setRange(0.0, 1.0);
        m_successRate->setSingleStep(0.05);
        m_successRate->setValue(0.8);
        netLayout->addWidget(m_successRate, 2, 1);

        m_effectivenessLabel = new QLabel;
        netLayout->addWidget(m_effectivenessLabel, 3, 0, 1, 2);

        layout->addWidget(netBox);

        // --- Calculate ---
        m_calcBtn = new QPushButton("Calculate REDD Emissions");
        layout->addWidget(m_calcBtn);

        layout->addStretch();

        // Connections
        auto updateEffectiveness = [this]() {
            double eff = m_successRate->value() - m_leakageRate->value();
            m_effectivenessLabel->setText(
                QString("Effectiveness = Success - Leakage = %1").arg(eff, 0, 'f', 2));
        };
        connect(m_leakageRate, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, updateEffectiveness);
        connect(m_successRate, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, updateEffectiveness);
        updateEffectiveness();

        connect(m_calcBtn, &QPushButton::clicked, this, [this]() {
            // Store parameters in session for potential future use
            m_session->set("lcm/redd_leakage", m_leakageRate->value());
            m_session->set("lcm/redd_success", m_successRate->value());
            m_session->set("lcm/redd_burn_fraction", m_burnFraction->value());
            m_session->set("lcm/redd_combustion_eff", m_combustionEff->value());
            emit statusMessage("REDD parameters saved. Full REDD calculation requires "
                               "prediction results and carbon pool data.");
        });
    }

    QString tabName() const override { return "REDD Project"; }
    QString tabDescription() const override {
        return "Estimate GHG emissions from REDD (Reducing Emissions from "
               "Deforestation and Forest Degradation) projects.";
    }

    // REDD is optional
    bool isComplete() const override { return true; }

private:
    QTableWidget*   m_poolTable       = nullptr;
    QDoubleSpinBox* m_burnFraction    = nullptr;
    QDoubleSpinBox* m_combustionEff   = nullptr;
    QDoubleSpinBox* m_leakageRate     = nullptr;
    QDoubleSpinBox* m_successRate     = nullptr;
    QLabel*         m_effectivenessLabel = nullptr;
    QPushButton*    m_calcBtn         = nullptr;
};

ModelerTab* createLcmReddTab(ModelerSession* s, QWidget* p) {
    return new LcmReddTab(s, p);
}

} // namespace aplaceholder

#include "LcmReddTab.moc"
