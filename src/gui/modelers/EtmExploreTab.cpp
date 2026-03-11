#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>
#include <QListWidget>
#include <QFileDialog>

namespace aplaceholder {

class EtmExploreTab : public ModelerTab {
    Q_OBJECT
public:
    explicit EtmExploreTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Session Parameters ---
        auto* sessBox = new QGroupBox("ETM Session Parameters");
        auto* sLayout = new QGridLayout(sessBox);

        sLayout->addWidget(new QLabel("Session file (.etm):"), 0, 0);
        m_sessionEdit = new QLineEdit;
        sLayout->addWidget(m_sessionEdit, 0, 1);
        auto* browseSession = new QPushButton("Browse...");
        sLayout->addWidget(browseSession, 0, 2);

        sLayout->addWidget(new QLabel("Add series:"), 1, 0);
        auto* addSeries = new QPushButton("Add Image/Index Series...");
        sLayout->addWidget(addSeries, 1, 1, 1, 2);

        sLayout->addWidget(new QLabel("Loaded series:"), 2, 0);
        m_seriesList = new QListWidget;
        m_seriesList->setMaximumHeight(100);
        sLayout->addWidget(m_seriesList, 2, 1, 1, 2);

        sLayout->addWidget(new QLabel("Mask (optional):"), 3, 0);
        m_maskEdit = new QLineEdit;
        sLayout->addWidget(m_maskEdit, 3, 1);
        auto* browseMask = new QPushButton("Browse...");
        sLayout->addWidget(browseMask, 3, 2);

        auto* loadBtn = new QPushButton("Load/Create Session");
        sLayout->addWidget(loadBtn, 4, 0, 1, 3);
        layout->addWidget(sessBox);

        // --- Explore Space/Time ---
        auto* exploreBox = new QGroupBox("Explore Space/Time Dynamics");
        auto* eLayout = new QGridLayout(exploreBox);

        eLayout->addWidget(new QLabel(
            "Image series: 4D visualization (X, Y, data value, time)\n"
            "Index series: Displayed as graphs\n"
            "Supports animation, Hovmoller plots, cross-correlation"), 0, 0, 1, 2);

        eLayout->addWidget(new QLabel("Display mode:"), 1, 0);
        m_displayMode = new QComboBox;
        m_displayMode->addItems({"Cube", "Plane", "Sphere"});
        eLayout->addWidget(m_displayMode, 1, 1);

        eLayout->addWidget(new QLabel("Animate over:"), 2, 0);
        m_animateAxis = new QComboBox;
        m_animateAxis->addItems({"Time", "X", "Y"});
        eLayout->addWidget(m_animateAxis, 2, 1);

        auto* visualize = new QPushButton("Visualize Series");
        eLayout->addWidget(visualize, 3, 0, 1, 2);
        layout->addWidget(exploreBox);

        // --- Temporal Profiles ---
        auto* profBox = new QGroupBox("Explore Temporal Profiles");
        auto* pLayout = new QGridLayout(profBox);

        pLayout->addWidget(new QLabel("Region of interest:"), 0, 0);
        m_roiType = new QComboBox;
        m_roiType->addItems({"Circular (specify radius)", "Vector polygon"});
        pLayout->addWidget(m_roiType, 0, 1);

        pLayout->addWidget(new QLabel("Summary statistic:"), 1, 0);
        m_statType = new QComboBox;
        m_statType->addItems({"Mean", "Median", "Min", "Max", "Range", "Sum", "Std Dev"});
        pLayout->addWidget(m_statType, 1, 1);

        auto* extractProfile = new QPushButton("Extract Temporal Profile");
        pLayout->addWidget(extractProfile, 2, 0, 1, 2);
        layout->addWidget(profBox);

        layout->addStretch();

        // Connections
        connect(browseSession, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select ETM Session",
                QString(), "ETM Sessions (*.etm);;All Files (*)");
            if (!path.isEmpty()) m_sessionEdit->setText(path);
        });
        connect(browseMask, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select Mask",
                QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
            if (!path.isEmpty()) m_maskEdit->setText(path);
        });
        connect(addSeries, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Add Series",
                QString(), "Group Files (*.rgf *.avl);;All Files (*)");
            if (!path.isEmpty()) m_seriesList->addItem(path);
        });
        connect(loadBtn, &QPushButton::clicked, this, [this]() {
            m_session->set("etm/session_file", m_sessionEdit->text().trimmed());
            m_loaded = true;
            emit statusMessage("ETM session loaded.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "Explore"; }
    QString tabDescription() const override {
        return "Session parameters and space/time exploration tools.";
    }
    bool isComplete() const override { return m_loaded; }

private:
    QLineEdit* m_sessionEdit = nullptr;
    QListWidget* m_seriesList = nullptr;
    QLineEdit* m_maskEdit = nullptr;
    QComboBox* m_displayMode = nullptr;
    QComboBox* m_animateAxis = nullptr;
    QComboBox* m_roiType = nullptr;
    QComboBox* m_statType = nullptr;
    bool m_loaded = false;
};

ModelerTab* createEtmExploreTab(ModelerSession* s, QWidget* p) {
    return new EtmExploreTab(s, p);
}

} // namespace aplaceholder

#include "EtmExploreTab.moc"
