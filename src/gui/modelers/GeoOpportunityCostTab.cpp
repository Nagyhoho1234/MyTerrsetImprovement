#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QComboBox>
#include <QTextEdit>

namespace aplaceholder {

class GeoOpportunityCostTab : public ModelerTab {
    Q_OBJECT
public:
    explicit GeoOpportunityCostTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        auto* box = new QGroupBox("Effective Opportunity Cost");
        auto* gLayout = new QGridLayout(box);

        gLayout->addWidget(new QLabel(
            "Run regression between deforestation and driver variables\n"
            "to compute effective opportunity cost scores."), 0, 0, 1, 2);

        gLayout->addWidget(new QLabel("Regression type:"), 1, 0);
        m_regType = new QComboBox;
        m_regType->addItems({"Logistic Regression (Boolean data)",
                             "Poisson Regression (Fractional data)",
                             "External Model (user coefficients)"});
        gLayout->addWidget(m_regType, 1, 1);

        gLayout->addWidget(new QLabel("Stratification:"), 2, 0);
        m_stratification = new QComboBox;
        m_stratification->addItems({"No stratification (single class)",
                                     "Quantile classification",
                                     "User-defined classification",
                                     "Geographic stratification (admin regions)"});
        gLayout->addWidget(m_stratification, 2, 1);

        auto* runBtn = new QPushButton("Run Regression");
        gLayout->addWidget(runBtn, 3, 0, 1, 2);

        gLayout->addWidget(new QLabel("Results:"), 4, 0);
        m_results = new QTextEdit;
        m_results->setReadOnly(true);
        m_results->setMaximumHeight(150);
        m_results->setPlaceholderText("Regression coefficients and opportunity cost scores will appear here...");
        gLayout->addWidget(m_results, 4, 1);

        layout->addWidget(box);
        layout->addStretch();

        connect(runBtn, &QPushButton::clicked, this, [this]() {
            m_session->set("geo/regression_type", m_regType->currentText());
            m_session->set("geo/stratification", m_stratification->currentText());
            m_results->setText("Regression type: " + m_regType->currentText() +
                "\nStratification: " + m_stratification->currentText() +
                "\n\nOC = (B0 + B1*X1 + ... + BN*XN) / B1\n\n"
                "Regression execution ready. Requires input images to be configured.");
            m_ran = true;
            emit statusMessage("Opportunity cost regression configured.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "Opportunity Cost"; }
    QString tabDescription() const override { return "Regression analysis for effective opportunity cost."; }
    bool isComplete() const override { return m_ran; }

private:
    QComboBox* m_regType = nullptr;
    QComboBox* m_stratification = nullptr;
    QTextEdit* m_results = nullptr;
    bool m_ran = false;
};

ModelerTab* createGeoOpportunityCostTab(ModelerSession* s, QWidget* p) {
    return new GeoOpportunityCostTab(s, p);
}

} // namespace aplaceholder

#include "GeoOpportunityCostTab.moc"
