#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QPushButton>

namespace aplaceholder {

class GeoModelParamsTab : public ModelerTab {
    Q_OBJECT
public:
    explicit GeoModelParamsTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        auto* box = new QGroupBox("Model Parameters");
        auto* gLayout = new QGridLayout(box);

        gLayout->addWidget(new QLabel(
            "Economic model parameters controlling the relationship\n"
            "between deforestation change and agricultural prices."), 0, 0, 1, 2);

        gLayout->addWidget(new QLabel("Price elasticity:"), 1, 0);
        m_elasticity = new QDoubleSpinBox;
        m_elasticity->setRange(-10, 10);
        m_elasticity->setValue(0.5);
        m_elasticity->setDecimals(4);
        m_elasticity->setSingleStep(0.01);
        gLayout->addWidget(m_elasticity, 1, 1);

        gLayout->addWidget(new QLabel("Exogenous increase in ag price:"), 2, 0);
        m_exoPrice = new QDoubleSpinBox;
        m_exoPrice->setRange(0, 100);
        m_exoPrice->setValue(0.0);
        m_exoPrice->setDecimals(4);
        m_exoPrice->setSingleStep(0.01);
        gLayout->addWidget(m_exoPrice, 2, 1);

        gLayout->addWidget(new QLabel("Soil carbon fraction:"), 3, 0);
        m_soilFrac = new QDoubleSpinBox;
        m_soilFrac->setRange(0, 1);
        m_soilFrac->setValue(0.25);
        m_soilFrac->setDecimals(3);
        gLayout->addWidget(m_soilFrac, 3, 1);

        gLayout->addWidget(new QLabel("Peat emission factor (tCO2/ha):"), 4, 0);
        m_peatFactor = new QDoubleSpinBox;
        m_peatFactor->setRange(0, 10000);
        m_peatFactor->setValue(55.0);
        m_peatFactor->setDecimals(1);
        gLayout->addWidget(m_peatFactor, 4, 1);

        auto* applyBtn = new QPushButton("Apply Model Parameters");
        gLayout->addWidget(applyBtn, 5, 0, 1, 2);
        layout->addWidget(box);

        layout->addStretch();

        connect(applyBtn, &QPushButton::clicked, this, [this]() {
            m_session->set("geo/price_elasticity", m_elasticity->value());
            m_session->set("geo/exo_price_increase", m_exoPrice->value());
            m_session->set("geo/soil_carbon_fraction", m_soilFrac->value());
            m_session->set("geo/peat_emission_factor", m_peatFactor->value());
            m_applied = true;
            emit statusMessage("Model parameters applied.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "Model Parameters"; }
    QString tabDescription() const override { return "Economic model parameters for REDD+ analysis."; }
    bool isComplete() const override { return m_applied; }

private:
    QDoubleSpinBox* m_elasticity = nullptr;
    QDoubleSpinBox* m_exoPrice = nullptr;
    QDoubleSpinBox* m_soilFrac = nullptr;
    QDoubleSpinBox* m_peatFactor = nullptr;
    bool m_applied = false;
};

ModelerTab* createGeoModelParamsTab(ModelerSession* s, QWidget* p) {
    return new GeoModelParamsTab(s, p);
}

} // namespace aplaceholder

#include "GeoModelParamsTab.moc"
