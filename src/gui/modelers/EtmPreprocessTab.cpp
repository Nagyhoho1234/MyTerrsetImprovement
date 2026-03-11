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
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QFileDialog>

namespace aplaceholder {

class EtmPreprocessTab : public ModelerTab {
    Q_OBJECT
public:
    explicit EtmPreprocessTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Missing Data Interpolation ---
        auto* interpBox = new QGroupBox("Missing Data Interpolation");
        auto* iLayout = new QGridLayout(interpBox);

        iLayout->addWidget(new QLabel("Input series:"), 0, 0);
        m_interpInput = new QLineEdit;
        iLayout->addWidget(m_interpInput, 0, 1);
        auto* browseInterp = new QPushButton("Browse...");
        iLayout->addWidget(browseInterp, 0, 2);

        iLayout->addWidget(new QLabel("Method:"), 1, 0);
        m_interpMethod = new QComboBox;
        m_interpMethod->addItems({"Harmonic (HANTS-based)", "Linear Temporal",
                                   "Spatial (3x3 median)", "Climatology (temporal median)"});
        iLayout->addWidget(m_interpMethod, 1, 1);

        iLayout->addWidget(new QLabel("Number of harmonics (HANTS):"), 2, 0);
        m_harmonics = new QSpinBox;
        m_harmonics->setRange(1, 20);
        m_harmonics->setValue(2);
        iLayout->addWidget(m_harmonics, 2, 1);

        iLayout->addWidget(new QLabel("Max data gap:"), 3, 0);
        m_maxGap = new QSpinBox;
        m_maxGap->setRange(1, 365);
        m_maxGap->setValue(30);
        iLayout->addWidget(m_maxGap, 3, 1);

        m_replaceAll = new QCheckBox("Replace ALL values (smoothing mode)");
        iLayout->addWidget(m_replaceAll, 4, 0, 1, 3);

        auto* runInterp = new QPushButton("Run Interpolation");
        iLayout->addWidget(runInterp, 5, 0, 1, 3);
        layout->addWidget(interpBox);

        // --- Denoise ---
        auto* denoiseBox = new QGroupBox("Denoise");
        auto* dLayout = new QGridLayout(denoiseBox);

        dLayout->addWidget(new QLabel("Method:"), 0, 0);
        m_denoiseMethod = new QComboBox;
        m_denoiseMethod->addItems({"Temporal Filter (Mean)", "Temporal Filter (Gaussian)",
                                    "Temporal Filter (Maximum)", "Temporal Filter (Cumulative Sum)",
                                    "Temporal Filter (Cumulative Mean)",
                                    "Maximum Value Composite",
                                    "Inverse PCA", "Inverse Fourier"});
        dLayout->addWidget(m_denoiseMethod, 0, 1);

        dLayout->addWidget(new QLabel("Filter window size:"), 1, 0);
        m_filterWindow = new QSpinBox;
        m_filterWindow->setRange(3, 365);
        m_filterWindow->setValue(5);
        dLayout->addWidget(m_filterWindow, 1, 1);

        dLayout->addWidget(new QLabel("Components to keep (Inverse PCA):"), 2, 0);
        m_pcaComponents = new QSpinBox;
        m_pcaComponents->setRange(1, 100);
        m_pcaComponents->setValue(20);
        dLayout->addWidget(m_pcaComponents, 2, 1);

        auto* runDenoise = new QPushButton("Run Denoise");
        dLayout->addWidget(runDenoise, 3, 0, 1, 2);
        layout->addWidget(denoiseBox);

        // --- Deseason ---
        auto* deseasonBox = new QGroupBox("Deseason");
        auto* dsLayout = new QGridLayout(deseasonBox);

        dsLayout->addWidget(new QLabel("Method:"), 0, 0);
        m_deseasonMethod = new QComboBox;
        m_deseasonMethod->addItems({"Anomalies", "Standardized Anomalies", "Temporal Filter"});
        dsLayout->addWidget(m_deseasonMethod, 0, 1);

        auto* runDeseason = new QPushButton("Run Deseason");
        dsLayout->addWidget(runDeseason, 1, 0, 1, 2);
        layout->addWidget(deseasonBox);

        // --- Detrend / Prewhiten ---
        auto* detrendBox = new QGroupBox("Detrend / Prewhiten");
        auto* dtLayout = new QGridLayout(detrendBox);

        dtLayout->addWidget(new QLabel("Method:"), 0, 0);
        m_detrendMethod = new QComboBox;
        m_detrendMethod->addItems({"Detrend (Linear)", "Detrend (Difference Series)",
                                    "Trend-Preserving Prewhitening (Wang-Swail)",
                                    "Durbin-Watson Statistic", "Cochrane-Orcutt"});
        dtLayout->addWidget(m_detrendMethod, 0, 1);

        auto* runDetrend = new QPushButton("Run Detrend/Prewhiten");
        dtLayout->addWidget(runDetrend, 1, 0, 1, 2);
        layout->addWidget(detrendBox);

        // --- Generate / Edit Series ---
        auto* genBox = new QGroupBox("Generate / Edit Series");
        auto* gLayout = new QGridLayout(genBox);

        gLayout->addWidget(new QLabel("Utility:"), 0, 0);
        m_genUtility = new QComboBox;
        m_genUtility->addItems({"Linear Index Series", "Sin Index Series", "Cos Index Series",
                                "Lagged Series", "Truncated Series", "Supplemented Series",
                                "Skip Factor Series", "Rename Series Images",
                                "Aggregate Series", "Subset Series (WINDOW)"});
        gLayout->addWidget(m_genUtility, 0, 1);

        auto* runGen = new QPushButton("Run Utility");
        gLayout->addWidget(runGen, 1, 0, 1, 2);
        layout->addWidget(genBox);

        layout->addStretch();

        // Connections
        connect(browseInterp, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select Series",
                QString(), "Group Files (*.rgf);;All Files (*)");
            if (!path.isEmpty()) m_interpInput->setText(path);
        });

        connect(runInterp, &QPushButton::clicked, this, [this]() {
            m_session->set("etm/interp_input", m_interpInput->text().trimmed());
            m_session->set("etm/interp_method", m_interpMethod->currentText());
            m_session->set("etm/harmonics", m_harmonics->value());
            m_session->set("etm/max_gap", m_maxGap->value());
            emit statusMessage("Interpolation configured: " + m_interpMethod->currentText());
        });

        connect(runDenoise, &QPushButton::clicked, this, [this]() {
            m_session->set("etm/denoise_method", m_denoiseMethod->currentText());
            m_session->set("etm/filter_window", m_filterWindow->value());
            m_session->set("etm/pca_components", m_pcaComponents->value());
            emit statusMessage("Denoise configured: " + m_denoiseMethod->currentText());
        });
    }

    QString tabName() const override { return "Preprocess"; }
    QString tabDescription() const override {
        return "Missing data interpolation, denoising, deseasoning, and detrending.";
    }
    bool isComplete() const override { return true; } // optional

private:
    QLineEdit* m_interpInput = nullptr;
    QComboBox* m_interpMethod = nullptr;
    QSpinBox* m_harmonics = nullptr;
    QSpinBox* m_maxGap = nullptr;
    QCheckBox* m_replaceAll = nullptr;
    QComboBox* m_denoiseMethod = nullptr;
    QSpinBox* m_filterWindow = nullptr;
    QSpinBox* m_pcaComponents = nullptr;
    QComboBox* m_deseasonMethod = nullptr;
    QComboBox* m_detrendMethod = nullptr;
    QComboBox* m_genUtility = nullptr;
};

ModelerTab* createEtmPreprocessTab(ModelerSession* s, QWidget* p) {
    return new EtmPreprocessTab(s, p);
}

} // namespace aplaceholder

#include "EtmPreprocessTab.moc"
