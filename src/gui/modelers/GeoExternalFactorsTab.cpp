#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QPushButton>

namespace aplaceholder {

class GeoExternalFactorsTab : public ModelerTab {
    Q_OBJECT
public:
    explicit GeoExternalFactorsTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        auto* box = new QGroupBox("External Factors");
        auto* gLayout = new QGridLayout(box);

        gLayout->addWidget(new QLabel("Carbon price ($/tCO2e):"), 0, 0);
        m_carbonPrice = new QDoubleSpinBox;
        m_carbonPrice->setRange(0, 1000);
        m_carbonPrice->setValue(5.0);
        m_carbonPrice->setDecimals(2);
        m_carbonPrice->setPrefix("$ ");
        gLayout->addWidget(m_carbonPrice, 0, 1);

        gLayout->addWidget(new QLabel(
            "Country national reference level\n"
            "(proportion of BAU emissions):"), 1, 0);
        m_refLevel = new QDoubleSpinBox;
        m_refLevel->setRange(0, 1);
        m_refLevel->setValue(1.0);
        m_refLevel->setDecimals(3);
        m_refLevel->setSingleStep(0.01);
        gLayout->addWidget(m_refLevel, 1, 1);

        auto* applyBtn = new QPushButton("Apply External Factors");
        gLayout->addWidget(applyBtn, 2, 0, 1, 2);
        layout->addWidget(box);

        layout->addStretch();

        connect(applyBtn, &QPushButton::clicked, this, [this]() {
            m_session->set("geo/carbon_price", m_carbonPrice->value());
            m_session->set("geo/ref_level", m_refLevel->value());
            m_applied = true;
            emit statusMessage("External factors applied.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "External Factors"; }
    QString tabDescription() const override { return "Carbon price and national reference level."; }
    bool isComplete() const override { return m_applied; }

private:
    QDoubleSpinBox* m_carbonPrice = nullptr;
    QDoubleSpinBox* m_refLevel = nullptr;
    bool m_applied = false;
};

ModelerTab* createGeoExternalFactorsTab(ModelerSession* s, QWidget* p) {
    return new GeoExternalFactorsTab(s, p);
}

} // namespace aplaceholder

#include "GeoExternalFactorsTab.moc"
