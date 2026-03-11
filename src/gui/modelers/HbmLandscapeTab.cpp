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
#include <QMessageBox>

namespace aplaceholder {

class HbmLandscapeTab : public ModelerTab {
    Q_OBJECT
public:
    explicit HbmLandscapeTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Landscape Pattern Analysis ---
        auto* patternBox = new QGroupBox("Landscape Pattern Analysis");
        auto* pLayout = new QGridLayout(patternBox);

        pLayout->addWidget(new QLabel("Land cover map:"), 0, 0);
        m_lcEdit = new QLineEdit;
        pLayout->addWidget(m_lcEdit, 0, 1);
        auto* browseLc = new QPushButton("Browse...");
        pLayout->addWidget(browseLc, 0, 2);

        pLayout->addWidget(new QLabel("Neighborhood size:"), 1, 0);
        m_neighborhood = new QComboBox;
        m_neighborhood->addItems({"3x3", "5x5", "7x7"});
        pLayout->addWidget(m_neighborhood, 1, 1);

        pLayout->addWidget(new QLabel("Metrics:"), 2, 0);
        m_entropy = new QCheckBox("Normalized Entropy (Shannon's)");
        m_entropy->setChecked(true);
        pLayout->addWidget(m_entropy, 2, 1);
        m_richness = new QCheckBox("Relative Richness");
        pLayout->addWidget(m_richness, 3, 1);
        m_edgeDensity = new QCheckBox("Edge Density");
        pLayout->addWidget(m_edgeDensity, 4, 1);
        m_patchArea = new QCheckBox("Patch Area");
        pLayout->addWidget(m_patchArea, 5, 1);
        m_compactness = new QCheckBox("Patch Compactness");
        pLayout->addWidget(m_compactness, 6, 1);

        auto* runPattern = new QPushButton("Run Pattern Analysis");
        pLayout->addWidget(runPattern, 7, 0, 1, 3);
        layout->addWidget(patternBox);

        // --- Landscape Change Process ---
        auto* changeBox = new QGroupBox("Landscape Change Process Analysis");
        auto* cLayout = new QGridLayout(changeBox);

        cLayout->addWidget(new QLabel(
            "Compare two land cover maps to classify change processes\n"
            "(deformation, shift, perforation, shrinkage, enlargement,\n"
            "attrition, aggregation, creation, dissection, fragmentation)."), 0, 0, 1, 3);

        cLayout->addWidget(new QLabel("Earlier land cover:"), 1, 0);
        m_lcEarlier = new QLineEdit;
        cLayout->addWidget(m_lcEarlier, 1, 1);
        auto* browseEarlier = new QPushButton("Browse...");
        cLayout->addWidget(browseEarlier, 1, 2);

        cLayout->addWidget(new QLabel("Later land cover:"), 2, 0);
        m_lcLater = new QLineEdit;
        cLayout->addWidget(m_lcLater, 2, 1);
        auto* browseLater = new QPushButton("Browse...");
        cLayout->addWidget(browseLater, 2, 2);

        auto* runChange = new QPushButton("Run Change Process Analysis");
        cLayout->addWidget(runChange, 3, 0, 1, 3);
        layout->addWidget(changeBox);

        // --- Sea Level Rise ---
        auto* slrBox = new QGroupBox("Sea Level Rise Impact");
        auto* sLayout = new QGridLayout(slrBox);

        sLayout->addWidget(new QLabel("DEM raster:"), 0, 0);
        m_demEdit = new QLineEdit;
        sLayout->addWidget(m_demEdit, 0, 1);
        auto* browseDem = new QPushButton("Browse...");
        sLayout->addWidget(browseDem, 0, 2);

        sLayout->addWidget(new QLabel("Projected sea level rise (m):"), 1, 0);
        m_slrThreshold = new QDoubleSpinBox;
        m_slrThreshold->setRange(0.01, 100);
        m_slrThreshold->setValue(1.0);
        m_slrThreshold->setDecimals(2);
        sLayout->addWidget(m_slrThreshold, 1, 1);

        sLayout->addWidget(new QLabel("RMSE of projection (m):"), 2, 0);
        m_slrRmse = new QDoubleSpinBox;
        m_slrRmse->setRange(0.01, 50);
        m_slrRmse->setValue(0.5);
        m_slrRmse->setDecimals(2);
        sLayout->addWidget(m_slrRmse, 2, 1);

        auto* runSlr = new QPushButton("Calculate Submersion Probability");
        sLayout->addWidget(runSlr, 3, 0, 1, 3);
        layout->addWidget(slrBox);

        layout->addStretch();

        // Browse connections
        auto browse = [this](QLineEdit* edit, const QString& title) {
            QString path = QFileDialog::getOpenFileName(this, title,
                QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
            if (!path.isEmpty()) edit->setText(path);
        };
        connect(browseLc, &QPushButton::clicked, this, [=]() { browse(m_lcEdit, "Select Land Cover Map"); });
        connect(browseEarlier, &QPushButton::clicked, this, [=]() { browse(m_lcEarlier, "Select Earlier LC Map"); });
        connect(browseLater, &QPushButton::clicked, this, [=]() { browse(m_lcLater, "Select Later LC Map"); });
        connect(browseDem, &QPushButton::clicked, this, [=]() { browse(m_demEdit, "Select DEM"); });

        connect(runPattern, &QPushButton::clicked, this, [this]() {
            if (m_lcEdit->text().trimmed().isEmpty()) {
                QMessageBox::warning(this, "Missing Input", "Please specify a land cover map.");
                return;
            }
            m_session->set("hbm/pattern_lc_map", m_lcEdit->text().trimmed());
            m_session->set("hbm/pattern_neighborhood", m_neighborhood->currentText());
            emit statusMessage("Pattern analysis parameters saved.");
        });

        connect(runSlr, &QPushButton::clicked, this, [this]() {
            if (m_demEdit->text().trimmed().isEmpty()) {
                QMessageBox::warning(this, "Missing Input", "Please specify a DEM raster.");
                return;
            }
            m_session->set("hbm/slr_dem", m_demEdit->text().trimmed());
            m_session->set("hbm/slr_threshold", m_slrThreshold->value());
            m_session->set("hbm/slr_rmse", m_slrRmse->value());
            emit statusMessage("Sea level rise parameters saved.");
        });
    }

    QString tabName() const override { return "Landscape Analysis"; }
    QString tabDescription() const override {
        return "Landscape pattern, change process analysis, and sea level rise impact.";
    }
    bool isComplete() const override { return true; } // optional

private:
    QLineEdit* m_lcEdit = nullptr;
    QComboBox* m_neighborhood = nullptr;
    QCheckBox* m_entropy = nullptr;
    QCheckBox* m_richness = nullptr;
    QCheckBox* m_edgeDensity = nullptr;
    QCheckBox* m_patchArea = nullptr;
    QCheckBox* m_compactness = nullptr;
    QLineEdit* m_lcEarlier = nullptr;
    QLineEdit* m_lcLater = nullptr;
    QLineEdit* m_demEdit = nullptr;
    QDoubleSpinBox* m_slrThreshold = nullptr;
    QDoubleSpinBox* m_slrRmse = nullptr;
};

ModelerTab* createHbmLandscapeTab(ModelerSession* s, QWidget* p) {
    return new HbmLandscapeTab(s, p);
}

} // namespace aplaceholder

#include "HbmLandscapeTab.moc"
