#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>
#include <QSpinBox>
#include <QCheckBox>
#include <QFileDialog>

namespace aplaceholder {

class EtmAnalysisTab : public ModelerTab {
    Q_OBJECT
public:
    explicit EtmAnalysisTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Series Trend Analysis ---
        auto* trendBox = new QGroupBox("Series Trend Analysis");
        auto* tLayout = new QGridLayout(trendBox);

        tLayout->addWidget(new QLabel("Input series (deseasoned recommended):"), 0, 0);
        m_trendInput = new QLineEdit;
        tLayout->addWidget(m_trendInput, 0, 1);
        auto* browseTrend = new QPushButton("Browse...");
        tLayout->addWidget(browseTrend, 0, 2);

        tLayout->addWidget(new QLabel("Trend type:"), 1, 0);
        m_trendType = new QComboBox;
        m_trendType->addItems({"Linearity (r-squared)", "Linear Correlation (Pearson)",
                               "Linear Trend (OLS slope)", "Median Trend (Theil-Sen)",
                               "Monotonic Trend (Mann-Kendall)", "Mann-Kendall Significance"});
        tLayout->addWidget(m_trendType, 1, 1);

        auto* runTrend = new QPushButton("Run Trend Analysis");
        tLayout->addWidget(runTrend, 2, 0, 1, 3);
        layout->addWidget(trendBox);

        // --- STA ---
        auto* staBox = new QGroupBox("Seasonal Trend Analysis (STA)");
        auto* sLayout = new QGridLayout(staBox);

        sLayout->addWidget(new QLabel(
            "Detect trends in seasonal progression character.\n"
            "Must NOT use deseasoned data."), 0, 0, 1, 3);

        sLayout->addWidget(new QLabel("Input series:"), 1, 0);
        m_staInput = new QLineEdit;
        sLayout->addWidget(m_staInput, 1, 1);
        auto* browseSta = new QPushButton("Browse...");
        sLayout->addWidget(browseSta, 1, 2);

        m_greenUpDown = new QCheckBox("Compute green up/down dates (phenology)");
        sLayout->addWidget(m_greenUpDown, 2, 0, 1, 3);

        auto* runSta = new QPushButton("Run STA");
        sLayout->addWidget(runSta, 3, 0, 1, 3);
        layout->addWidget(staBox);

        // --- Decomposition Methods ---
        auto* decompBox = new QGroupBox("Decomposition Methods");
        auto* dLayout = new QGridLayout(decompBox);

        dLayout->addWidget(new QLabel("Method:"), 0, 0);
        m_decompMethod = new QComboBox;
        m_decompMethod->addItems({"PCA/EOF (S-mode)", "PCA/EOF (T-mode)",
                                   "EPCA/EEOF (Extended PCA)", "MSSA",
                                   "EOT", "Cross-EOT (CEOT)", "Extended EOT (EEOT)",
                                   "MEOT (Multichannel EOT)", "TEOT (T-mode EOT)",
                                   "Fourier PCA", "CCA"});
        dLayout->addWidget(m_decompMethod, 0, 1);

        dLayout->addWidget(new QLabel("Embedding dimension (MSSA/MEOT):"), 1, 0);
        m_embedding = new QSpinBox;
        m_embedding->setRange(1, 999);
        m_embedding->setValue(13);
        dLayout->addWidget(m_embedding, 1, 1);

        dLayout->addWidget(new QLabel("Sampling rate (EOT family):"), 2, 0);
        m_sampling = new QSpinBox;
        m_sampling->setRange(1, 99);
        m_sampling->setValue(1);
        dLayout->addWidget(m_sampling, 2, 1);

        m_standardized = new QCheckBox("Standardized (correlation matrix)");
        m_standardized->setChecked(true);
        dLayout->addWidget(m_standardized, 3, 0, 1, 2);

        auto* runDecomp = new QPushButton("Run Decomposition");
        dLayout->addWidget(runDecomp, 4, 0, 1, 2);
        layout->addWidget(decompBox);

        // --- Linear Modeling ---
        auto* lmBox = new QGroupBox("Linear Modeling");
        auto* lLayout = new QGridLayout(lmBox);

        lLayout->addWidget(new QLabel("Dependent series (image):"), 0, 0);
        m_depInput = new QLineEdit;
        lLayout->addWidget(m_depInput, 0, 1);

        lLayout->addWidget(new QLabel("Independent series:"), 1, 0);
        m_indepInput = new QLineEdit;
        lLayout->addWidget(m_indepInput, 1, 1);

        m_partialCorr = new QCheckBox("Partial correlation");
        lLayout->addWidget(m_partialCorr, 2, 0);

        dLayout->addWidget(new QLabel("Lag (negative=lead):"), 3, 0);
        m_lag = new QSpinBox;
        m_lag->setRange(-100, 100);
        m_lag->setValue(0);
        lLayout->addWidget(m_lag, 3, 1);

        auto* runLm = new QPushButton("Run Linear Model");
        lLayout->addWidget(runLm, 4, 0, 1, 2);
        layout->addWidget(lmBox);

        layout->addStretch();

        // Connections
        connect(browseTrend, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select Series",
                QString(), "Group Files (*.rgf);;All Files (*)");
            if (!path.isEmpty()) m_trendInput->setText(path);
        });
        connect(browseSta, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select Series",
                QString(), "Group Files (*.rgf);;All Files (*)");
            if (!path.isEmpty()) m_staInput->setText(path);
        });

        connect(runTrend, &QPushButton::clicked, this, [this]() {
            m_session->set("etm/trend_input", m_trendInput->text().trimmed());
            m_session->set("etm/trend_type", m_trendType->currentText());
            emit statusMessage("Trend analysis configured: " + m_trendType->currentText());
        });

        connect(runDecomp, &QPushButton::clicked, this, [this]() {
            m_session->set("etm/decomp_method", m_decompMethod->currentText());
            m_session->set("etm/embedding_dim", m_embedding->value());
            m_session->set("etm/sampling_rate", m_sampling->value());
            m_session->set("etm/standardized", m_standardized->isChecked());
            emit statusMessage("Decomposition configured: " + m_decompMethod->currentText());
        });
    }

    QString tabName() const override { return "Analysis"; }
    QString tabDescription() const override {
        return "Trend analysis, STA, decomposition methods, and linear modeling.";
    }
    bool isComplete() const override { return true; } // optional

private:
    QLineEdit* m_trendInput = nullptr;
    QComboBox* m_trendType = nullptr;
    QLineEdit* m_staInput = nullptr;
    QCheckBox* m_greenUpDown = nullptr;
    QComboBox* m_decompMethod = nullptr;
    QSpinBox* m_embedding = nullptr;
    QSpinBox* m_sampling = nullptr;
    QCheckBox* m_standardized = nullptr;
    QLineEdit* m_depInput = nullptr;
    QLineEdit* m_indepInput = nullptr;
    QCheckBox* m_partialCorr = nullptr;
    QSpinBox* m_lag = nullptr;
};

ModelerTab* createEtmAnalysisTab(ModelerSession* s, QWidget* p) {
    return new EtmAnalysisTab(s, p);
}

} // namespace aplaceholder

#include "EtmAnalysisTab.moc"
