#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QPushButton>

namespace aplaceholder {

class GeoReddRulesTab : public ModelerTab {
    Q_OBJECT
public:
    explicit GeoReddRulesTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        auto* box = new QGroupBox("Decision on National REDD+ Rules and Incentives");
        auto* gLayout = new QGridLayout(box);

        gLayout->addWidget(new QLabel(
            "Configure benefit and cost sharing between national\n"
            "and sub-national levels."), 0, 0, 1, 2);

        gLayout->addWidget(new QLabel("Benefit sharing (proportion to sub-national):"), 1, 0);
        m_benefitShare = new QDoubleSpinBox;
        m_benefitShare->setRange(0, 1);
        m_benefitShare->setValue(0.5);
        m_benefitShare->setDecimals(3);
        m_benefitShare->setSingleStep(0.01);
        gLayout->addWidget(m_benefitShare, 1, 1);

        gLayout->addWidget(new QLabel("Cost sharing (proportion to sub-national):"), 2, 0);
        m_costShare = new QDoubleSpinBox;
        m_costShare->setRange(0, 1);
        m_costShare->setValue(0.5);
        m_costShare->setDecimals(3);
        m_costShare->setSingleStep(0.01);
        gLayout->addWidget(m_costShare, 2, 1);

        auto* applyBtn = new QPushButton("Apply REDD+ Rules");
        gLayout->addWidget(applyBtn, 3, 0, 1, 2);
        layout->addWidget(box);

        layout->addStretch();

        connect(applyBtn, &QPushButton::clicked, this, [this]() {
            m_session->set("geo/benefit_share", m_benefitShare->value());
            m_session->set("geo/cost_share", m_costShare->value());
            m_applied = true;
            emit statusMessage("REDD+ rules applied.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "REDD+ Rules"; }
    QString tabDescription() const override { return "National REDD+ benefit and cost sharing rules."; }
    bool isComplete() const override { return m_applied; }

private:
    QDoubleSpinBox* m_benefitShare = nullptr;
    QDoubleSpinBox* m_costShare = nullptr;
    bool m_applied = false;
};

ModelerTab* createGeoReddRulesTab(ModelerSession* s, QWidget* p) {
    return new GeoReddRulesTab(s, p);
}

} // namespace aplaceholder

#include "GeoReddRulesTab.moc"
