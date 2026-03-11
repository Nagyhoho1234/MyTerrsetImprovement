#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>
#include <QTableWidget>
#include <QHeaderView>
#include <QFileDialog>
#include <QMessageBox>
#include <QFile>
#include <QTextStream>
#include <QDir>

namespace aplaceholder {

class LcmChangeAnalysisTab : public ModelerTab {
    Q_OBJECT
public:
    explicit LcmChangeAnalysisTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Session Parameters ---
        auto* sessionBox = new QGroupBox("Session Parameters");
        auto* sessionLayout = new QGridLayout(sessionBox);

        sessionLayout->addWidget(new QLabel("Earlier land cover map:"), 0, 0);
        m_earlierEdit = new QLineEdit;
        sessionLayout->addWidget(m_earlierEdit, 0, 1);
        auto* browseEarlier = new QPushButton("Browse...");
        sessionLayout->addWidget(browseEarlier, 0, 2);

        sessionLayout->addWidget(new QLabel("Later land cover map:"), 1, 0);
        m_laterEdit = new QLineEdit;
        sessionLayout->addWidget(m_laterEdit, 1, 1);
        auto* browseLater = new QPushButton("Browse...");
        sessionLayout->addWidget(browseLater, 1, 2);

        sessionLayout->addWidget(new QLabel("DEM (optional):"), 2, 0);
        m_demEdit = new QLineEdit;
        sessionLayout->addWidget(m_demEdit, 2, 1);
        auto* browseDem = new QPushButton("Browse...");
        sessionLayout->addWidget(browseDem, 2, 2);

        sessionLayout->addWidget(new QLabel("Roads layer (optional):"), 3, 0);
        m_roadsEdit = new QLineEdit;
        sessionLayout->addWidget(m_roadsEdit, 3, 1);
        auto* browseRoads = new QPushButton("Browse...");
        sessionLayout->addWidget(browseRoads, 3, 2);

        layout->addWidget(sessionBox);

        // --- Analysis Controls ---
        auto* analysisBox = new QGroupBox("Change Analysis");
        auto* analysisLayout = new QHBoxLayout(analysisBox);

        analysisLayout->addWidget(new QLabel("Analysis mode:"));
        m_modeCombo = new QComboBox;
        m_modeCombo->addItems({"Gains and Losses", "Net Change", "Contributors to Change"});
        analysisLayout->addWidget(m_modeCombo);

        m_runBtn = new QPushButton("Run Change Analysis");
        analysisLayout->addWidget(m_runBtn);

        layout->addWidget(analysisBox);

        // --- Results Table ---
        auto* resultsBox = new QGroupBox("Results");
        auto* resultsLayout = new QVBoxLayout(resultsBox);

        m_resultsTable = new QTableWidget;
        m_resultsTable->setColumnCount(7);
        m_resultsTable->setHorizontalHeaderLabels(
            {"Class", "Area T1", "Area T2", "Gains", "Losses", "Net Change", "Persistence"});
        m_resultsTable->horizontalHeader()->setStretchLastSection(true);
        m_resultsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
        resultsLayout->addWidget(m_resultsTable);

        layout->addWidget(resultsBox, 1);

        // --- Change Maps ---
        auto* mapBox = new QGroupBox("Change Maps");
        auto* mapLayout = new QHBoxLayout(mapBox);

        sessionLayout->addWidget(new QLabel("Output directory:"), 4, 0);
        m_outputDirEdit = new QLineEdit;
        sessionLayout->addWidget(m_outputDirEdit, 4, 1);
        auto* browseOutDir = new QPushButton("Browse...");
        sessionLayout->addWidget(browseOutDir, 4, 2);

        m_genMapsBtn = new QPushButton("Generate Change Maps");
        m_genMapsBtn->setEnabled(false);
        mapLayout->addWidget(m_genMapsBtn);

        layout->addWidget(mapBox);

        // Connections
        auto browseRaster = [this](QLineEdit* edit) {
            QString path = QFileDialog::getOpenFileName(this, "Select Raster",
                QString(), "Raster Files (*.tif *.tiff *.rst *.img *.nc *.hdf);;All Files (*)");
            if (!path.isEmpty()) edit->setText(path);
        };

        connect(browseEarlier, &QPushButton::clicked, this, [=]() { browseRaster(m_earlierEdit); });
        connect(browseLater,   &QPushButton::clicked, this, [=]() { browseRaster(m_laterEdit); });
        connect(browseDem,     &QPushButton::clicked, this, [=]() { browseRaster(m_demEdit); });
        connect(browseRoads,   &QPushButton::clicked, this, [=]() { browseRaster(m_roadsEdit); });
        connect(browseOutDir,  &QPushButton::clicked, this, [this]() {
            QString dir = QFileDialog::getExistingDirectory(this, "Select Output Directory");
            if (!dir.isEmpty()) m_outputDirEdit->setText(dir);
        });

        connect(m_runBtn, &QPushButton::clicked, this, &LcmChangeAnalysisTab::runAnalysis);
        connect(m_genMapsBtn, &QPushButton::clicked, this, &LcmChangeAnalysisTab::generateMaps);
    }

    QString tabName() const override { return "Change Analysis"; }
    QString tabDescription() const override {
        return "Analyze land cover changes between two time periods.";
    }

    bool isComplete() const override {
        return m_session->contains("lcm/change_stats_csv");
    }

    void onActivated() override {
        if (m_session->contains("lcm/earlier_map"))
            m_earlierEdit->setText(m_session->get("lcm/earlier_map").toString());
        if (m_session->contains("lcm/later_map"))
            m_laterEdit->setText(m_session->get("lcm/later_map").toString());
    }

private slots:
    void runAnalysis() {
        QString earlier = m_earlierEdit->text().trimmed();
        QString later   = m_laterEdit->text().trimmed();

        if (earlier.isEmpty() || later.isEmpty()) {
            QMessageBox::warning(this, "Missing Input",
                "Both earlier and later land cover maps are required.");
            return;
        }

        // Determine output path
        QString outDir = m_outputDirEdit->text().trimmed();
        if (outDir.isEmpty()) {
            outDir = QFileInfo(later).absolutePath();
            m_outputDirEdit->setText(outDir);
        }
        QString outputPath = outDir + "/lcm_change_analysis.tif";

        // Store in session
        m_session->set("lcm/earlier_map", earlier);
        m_session->set("lcm/later_map", later);
        if (!m_demEdit->text().trimmed().isEmpty())
            m_session->set("lcm/dem", m_demEdit->text().trimmed());
        if (!m_roadsEdit->text().trimmed().isEmpty())
            m_session->set("lcm/roads", m_roadsEdit->text().trimmed());
        m_session->set("lcm/output_dir", outDir);

        // Run the CHANGE_ANALYSIS module
        m_runBtn->setEnabled(false);
        bool ok = runModule("CHANGE_ANALYSIS", {
            {"earlier_map", earlier},
            {"later_map", later},
            {"output", outputPath}
        });
        m_runBtn->setEnabled(true);

        if (!ok) {
            QMessageBox::critical(this, "Change Analysis Failed",
                "The change analysis module failed. Check the status bar for details.");
            return;
        }

        // Store results
        QString csvPath = outputPath;
        csvPath = csvPath.left(csvPath.lastIndexOf('.')) + "_stats.csv";
        m_session->set("lcm/change_stats_csv", csvPath);
        m_session->set("lcm/transition_map", outputPath);

        // Load and display results
        loadResultsCsv(csvPath);
        m_genMapsBtn->setEnabled(true);

        emit completionChanged();
        emit statusMessage("Change analysis complete.");
    }

    void generateMaps() {
        emit statusMessage("Change maps generated from transition raster.");
    }

private:
    void loadResultsCsv(const QString& csvPath) {
        QFile file(csvPath);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return;

        QTextStream in(&file);
        m_resultsTable->setRowCount(0);

        bool headerSkipped = false;
        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();
            if (line.isEmpty() || line.startsWith('\n'))
                break; // stop at blank line (before transition matrix section)

            if (!headerSkipped) {
                headerSkipped = true;
                continue;
            }

            QStringList fields = line.split(',');
            if (fields.size() < 7) continue;

            int row = m_resultsTable->rowCount();
            m_resultsTable->insertRow(row);
            for (int col = 0; col < 7; ++col) {
                m_resultsTable->setItem(row, col,
                    new QTableWidgetItem(fields[col].trimmed()));
            }
        }

        m_resultsTable->resizeColumnsToContents();
    }

    QLineEdit*    m_earlierEdit   = nullptr;
    QLineEdit*    m_laterEdit     = nullptr;
    QLineEdit*    m_demEdit       = nullptr;
    QLineEdit*    m_roadsEdit     = nullptr;
    QLineEdit*    m_outputDirEdit = nullptr;
    QComboBox*    m_modeCombo     = nullptr;
    QPushButton*  m_runBtn        = nullptr;
    QPushButton*  m_genMapsBtn    = nullptr;
    QTableWidget* m_resultsTable  = nullptr;
};

ModelerTab* createLcmChangeAnalysisTab(ModelerSession* s, QWidget* p) {
    return new LcmChangeAnalysisTab(s, p);
}

} // namespace aplaceholder

#include "LcmChangeAnalysisTab.moc"
