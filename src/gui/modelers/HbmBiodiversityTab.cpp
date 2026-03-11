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
#include <QMessageBox>

namespace aplaceholder {

class HbmBiodiversityTab : public ModelerTab {
    Q_OBJECT
public:
    explicit HbmBiodiversityTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- IUCN Species Ranges ---
        auto* iucnBox = new QGroupBox("Subset IUCN Species Ranges");
        auto* iLayout = new QGridLayout(iucnBox);

        iLayout->addWidget(new QLabel("IUCN species data file:"), 0, 0);
        m_iucnEdit = new QLineEdit;
        iLayout->addWidget(m_iucnEdit, 0, 1);
        auto* browseIucn = new QPushButton("Browse...");
        iLayout->addWidget(browseIucn, 0, 2);

        iLayout->addWidget(new QLabel("Filter by Red List status:"), 1, 0);
        m_statusFilter = new QComboBox;
        m_statusFilter->addItems({"All", "Least Concerned", "Near Threatened", "Vulnerable",
                                  "Endangered", "Critically Endangered", "Extinct in Wild",
                                  "Extinct", "Data Deficient"});
        iLayout->addWidget(m_statusFilter, 1, 1);

        auto* loadIucn = new QPushButton("Load Species Ranges");
        iLayout->addWidget(loadIucn, 2, 0, 1, 3);
        layout->addWidget(iucnBox);

        // --- Biodiversity Analysis ---
        auto* bioBox = new QGroupBox("Biodiversity Analysis");
        auto* bLayout = new QGridLayout(bioBox);

        bLayout->addWidget(new QLabel("Species range input:"), 0, 0);
        m_rangeEdit = new QLineEdit;
        bLayout->addWidget(m_rangeEdit, 0, 1);
        auto* browseRange = new QPushButton("Browse...");
        bLayout->addWidget(browseRange, 0, 2);

        bLayout->addWidget(new QLabel("Input format:"), 1, 0);
        m_inputFormat = new QComboBox;
        m_inputFormat->addItems({"Vector composite polygon", "Vector group file", "Raster group file"});
        bLayout->addWidget(m_inputFormat, 1, 1);

        bLayout->addWidget(new QLabel("Metrics:"), 2, 0);
        m_alphaDiv = new QCheckBox("Alpha Diversity (species richness)");
        m_alphaDiv->setChecked(true);
        bLayout->addWidget(m_alphaDiv, 2, 1);
        m_gammaDiv = new QCheckBox("Gamma Diversity (regional richness)");
        bLayout->addWidget(m_gammaDiv, 3, 1);
        m_betaDiv = new QCheckBox("Beta Diversity (Whittaker's)");
        bLayout->addWidget(m_betaDiv, 4, 1);
        m_sorensen = new QCheckBox("Sorensen Dissimilarity");
        bLayout->addWidget(m_sorensen, 5, 1);
        m_rri = new QCheckBox("Range Restriction Index (endemism)");
        bLayout->addWidget(m_rri, 6, 1);

        bLayout->addWidget(new QLabel("Regional definition:"), 7, 0);
        m_regionType = new QComboBox;
        m_regionType->addItems({"Vector region polygon", "Raster region polygon",
                                "Focal zone (circular)"});
        bLayout->addWidget(m_regionType, 7, 1);

        bLayout->addWidget(new QLabel("Focal zone diameter:"), 8, 0);
        m_focalDiam = new QSpinBox;
        m_focalDiam->setRange(3, 999);
        m_focalDiam->setValue(11);
        m_focalDiam->setSuffix(" pixels");
        bLayout->addWidget(m_focalDiam, 8, 1);

        auto* runBio = new QPushButton("Run Biodiversity Analysis");
        bLayout->addWidget(runBio, 9, 0, 1, 3);
        layout->addWidget(bioBox);

        layout->addStretch();

        // Connections
        connect(browseIucn, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select IUCN Data",
                QString(), "All Files (*)");
            if (!path.isEmpty()) m_iucnEdit->setText(path);
        });
        connect(browseRange, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select Species Range Data",
                QString(), "All Files (*)");
            if (!path.isEmpty()) m_rangeEdit->setText(path);
        });

        connect(runBio, &QPushButton::clicked, this, [this]() {
            m_session->set("hbm/range_input", m_rangeEdit->text().trimmed());
            m_session->set("hbm/input_format", m_inputFormat->currentText());
            m_session->set("hbm/alpha_div", m_alphaDiv->isChecked());
            m_session->set("hbm/gamma_div", m_gammaDiv->isChecked());
            m_session->set("hbm/beta_div", m_betaDiv->isChecked());
            m_session->set("hbm/sorensen", m_sorensen->isChecked());
            m_session->set("hbm/rri", m_rri->isChecked());
            m_session->set("hbm/region_type", m_regionType->currentText());
            m_session->set("hbm/focal_diameter", m_focalDiam->value());
            emit statusMessage("Biodiversity analysis parameters saved.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "Biodiversity"; }
    QString tabDescription() const override {
        return "IUCN species range subsetting and biodiversity metrics.";
    }
    bool isComplete() const override { return true; } // optional

private:
    QLineEdit* m_iucnEdit = nullptr;
    QComboBox* m_statusFilter = nullptr;
    QLineEdit* m_rangeEdit = nullptr;
    QComboBox* m_inputFormat = nullptr;
    QCheckBox* m_alphaDiv = nullptr;
    QCheckBox* m_gammaDiv = nullptr;
    QCheckBox* m_betaDiv = nullptr;
    QCheckBox* m_sorensen = nullptr;
    QCheckBox* m_rri = nullptr;
    QComboBox* m_regionType = nullptr;
    QSpinBox* m_focalDiam = nullptr;
};

ModelerTab* createHbmBiodiversityTab(ModelerSession* s, QWidget* p) {
    return new HbmBiodiversityTab(s, p);
}

} // namespace aplaceholder

#include "HbmBiodiversityTab.moc"
