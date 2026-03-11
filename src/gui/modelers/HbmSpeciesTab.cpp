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
#include <QFileDialog>
#include <QMessageBox>

namespace aplaceholder {

class HbmSpeciesTab : public ModelerTab {
    Q_OBJECT
public:
    explicit HbmSpeciesTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Habitat Assessment ---
        auto* habitatBox = new QGroupBox("Habitat Assessment");
        auto* hLayout = new QGridLayout(habitatBox);

        hLayout->addWidget(new QLabel("Land cover map:"), 0, 0);
        m_lcMapEdit = new QLineEdit;
        hLayout->addWidget(m_lcMapEdit, 0, 1);
        auto* browseLc = new QPushButton("Browse...");
        hLayout->addWidget(browseLc, 0, 2);

        hLayout->addWidget(new QLabel("Habitat suitability map (optional):"), 1, 0);
        m_suitEdit = new QLineEdit;
        hLayout->addWidget(m_suitEdit, 1, 1);
        auto* browseSuit = new QPushButton("Browse...");
        hLayout->addWidget(browseSuit, 1, 2);

        hLayout->addWidget(new QLabel("Gap distance within range:"), 2, 0);
        m_gapWithin = new QDoubleSpinBox;
        m_gapWithin->setRange(0, 100000);
        m_gapWithin->setValue(500);
        hLayout->addWidget(m_gapWithin, 2, 1);

        hLayout->addWidget(new QLabel("Gap distance outside range:"), 3, 0);
        m_gapOutside = new QDoubleSpinBox;
        m_gapOutside->setRange(0, 100000);
        m_gapOutside->setValue(2000);
        hLayout->addWidget(m_gapOutside, 3, 1);

        hLayout->addWidget(new QLabel("Minimum core area:"), 4, 0);
        m_minCore = new QDoubleSpinBox;
        m_minCore->setRange(0, 1e9);
        m_minCore->setValue(1000);
        hLayout->addWidget(m_minCore, 4, 1);

        hLayout->addWidget(new QLabel("Minimum edge buffer:"), 5, 0);
        m_minBuffer = new QDoubleSpinBox;
        m_minBuffer->setRange(0, 100000);
        m_minBuffer->setValue(100);
        hLayout->addWidget(m_minBuffer, 5, 1);

        auto* runHabitat = new QPushButton("Run Habitat Assessment");
        hLayout->addWidget(runHabitat, 6, 0, 1, 3);
        layout->addWidget(habitatBox);

        // --- Species Distribution Modeling ---
        auto* sdmBox = new QGroupBox("Species Distribution / Habitat Suitability");
        auto* sLayout = new QGridLayout(sdmBox);

        sLayout->addWidget(new QLabel("Training data type:"), 0, 0);
        m_trainingType = new QComboBox;
        m_trainingType->addItems({"None (MCE)", "Presence Only", "Presence/Absence", "Abundance"});
        sLayout->addWidget(m_trainingType, 0, 1);

        sLayout->addWidget(new QLabel("Method:"), 1, 0);
        m_method = new QComboBox;
        m_method->addItems({"MCE (Weighted Linear Combination)", "MaxEnt", "Mahalanobis Typicality",
                           "Weighted Mahalanobis", "MLP Neural Network", "Logistic Regression",
                           "Multiple Regression"});
        sLayout->addWidget(m_method, 1, 1);

        sLayout->addWidget(new QLabel("Training sites:"), 2, 0);
        m_trainEdit = new QLineEdit;
        sLayout->addWidget(m_trainEdit, 2, 1);
        auto* browseTrain = new QPushButton("Browse...");
        sLayout->addWidget(browseTrain, 2, 2);

        auto* runSdm = new QPushButton("Run Species Modeling");
        sLayout->addWidget(runSdm, 3, 0, 1, 3);
        layout->addWidget(sdmBox);

        // --- Bioclimatic Variables ---
        auto* bioBox = new QGroupBox("Bioclimatic Variables (19 BioClim)");
        auto* bLayout = new QGridLayout(bioBox);

        bLayout->addWidget(new QLabel(
            "Generate 19 biologically meaningful variables from monthly\n"
            "temperature and precipitation data."), 0, 0, 1, 3);

        bLayout->addWidget(new QLabel("Monthly temperature series:"), 1, 0);
        m_tempEdit = new QLineEdit;
        bLayout->addWidget(m_tempEdit, 1, 1);
        auto* browseTemp = new QPushButton("Browse...");
        bLayout->addWidget(browseTemp, 1, 2);

        bLayout->addWidget(new QLabel("Monthly precipitation series:"), 2, 0);
        m_precipEdit = new QLineEdit;
        bLayout->addWidget(m_precipEdit, 2, 1);
        auto* browsePrec = new QPushButton("Browse...");
        bLayout->addWidget(browsePrec, 2, 2);

        auto* runBio = new QPushButton("Generate Bioclimatic Variables");
        bLayout->addWidget(runBio, 3, 0, 1, 3);
        layout->addWidget(bioBox);

        layout->addStretch();

        // Connections
        auto browse = [this](QLineEdit* edit, const QString& title) {
            QString path = QFileDialog::getOpenFileName(this, title,
                QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
            if (!path.isEmpty()) edit->setText(path);
        };
        connect(browseLc, &QPushButton::clicked, this, [=]() { browse(m_lcMapEdit, "Select Land Cover Map"); });
        connect(browseSuit, &QPushButton::clicked, this, [=]() { browse(m_suitEdit, "Select Suitability Map"); });
        connect(browseTrain, &QPushButton::clicked, this, [=]() { browse(m_trainEdit, "Select Training Sites"); });
        connect(browseTemp, &QPushButton::clicked, this, [=]() { browse(m_tempEdit, "Select Temperature Series"); });
        connect(browsePrec, &QPushButton::clicked, this, [=]() { browse(m_precipEdit, "Select Precipitation Series"); });

        connect(runHabitat, &QPushButton::clicked, this, [this]() {
            if (m_lcMapEdit->text().trimmed().isEmpty()) {
                QMessageBox::warning(this, "Missing Input", "Please specify a land cover map.");
                return;
            }
            m_session->set("hbm/landcover_map", m_lcMapEdit->text().trimmed());
            m_session->set("hbm/suitability_map", m_suitEdit->text().trimmed());
            m_session->set("hbm/gap_within", m_gapWithin->value());
            m_session->set("hbm/gap_outside", m_gapOutside->value());
            m_session->set("hbm/min_core_area", m_minCore->value());
            m_session->set("hbm/min_edge_buffer", m_minBuffer->value());
            m_assessed = true;
            emit statusMessage("Habitat assessment parameters saved.");
            emit completionChanged();
        });

        connect(runSdm, &QPushButton::clicked, this, [this]() {
            m_session->set("hbm/training_type", m_trainingType->currentText());
            m_session->set("hbm/sdm_method", m_method->currentText());
            emit statusMessage("Species distribution modeling parameters saved.");
        });

        connect(runBio, &QPushButton::clicked, this, [this]() {
            if (m_tempEdit->text().trimmed().isEmpty() || m_precipEdit->text().trimmed().isEmpty()) {
                QMessageBox::warning(this, "Missing Input", "Please specify both temperature and precipitation series.");
                return;
            }
            m_session->set("hbm/temp_series", m_tempEdit->text().trimmed());
            m_session->set("hbm/precip_series", m_precipEdit->text().trimmed());
            emit statusMessage("Bioclimatic variable generation parameters saved.");
        });
    }

    QString tabName() const override { return "Species"; }
    QString tabDescription() const override {
        return "Habitat assessment, species modeling, and bioclimatic variables.";
    }
    bool isComplete() const override { return m_assessed; }

private:
    QLineEdit* m_lcMapEdit = nullptr;
    QLineEdit* m_suitEdit = nullptr;
    QDoubleSpinBox* m_gapWithin = nullptr;
    QDoubleSpinBox* m_gapOutside = nullptr;
    QDoubleSpinBox* m_minCore = nullptr;
    QDoubleSpinBox* m_minBuffer = nullptr;
    QComboBox* m_trainingType = nullptr;
    QComboBox* m_method = nullptr;
    QLineEdit* m_trainEdit = nullptr;
    QLineEdit* m_tempEdit = nullptr;
    QLineEdit* m_precipEdit = nullptr;
    bool m_assessed = false;
};

ModelerTab* createHbmSpeciesTab(ModelerSession* s, QWidget* p) {
    return new HbmSpeciesTab(s, p);
}

} // namespace aplaceholder

#include "HbmSpeciesTab.moc"
