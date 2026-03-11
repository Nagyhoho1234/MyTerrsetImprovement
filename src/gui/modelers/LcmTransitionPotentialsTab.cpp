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
#include <QSpinBox>

namespace aplaceholder {

class LcmTransitionPotentialsTab : public ModelerTab {
    Q_OBJECT
public:
    explicit LcmTransitionPotentialsTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Transition Sub-Models ---
        auto* subModelBox = new QGroupBox("Transition Sub-Models");
        auto* subLayout = new QVBoxLayout(subModelBox);

        subLayout->addWidget(new QLabel(
            "Group transitions into sub-models. Transitions sharing the same "
            "driving forces should be grouped together."));

        m_transitionsTable = new QTableWidget;
        m_transitionsTable->setColumnCount(3);
        m_transitionsTable->setHorizontalHeaderLabels({"From", "To", "Sub-Model"});
        m_transitionsTable->horizontalHeader()->setStretchLastSection(true);
        subLayout->addWidget(m_transitionsTable);

        layout->addWidget(subModelBox);

        // --- Driver Variables ---
        auto* driversBox = new QGroupBox("Explanatory (Driver) Variables");
        auto* driversLayout = new QGridLayout(driversBox);

        driversLayout->addWidget(new QLabel("Driver rasters (comma-separated):"), 0, 0);
        m_driversEdit = new QLineEdit;
        driversLayout->addWidget(m_driversEdit, 0, 1);
        auto* browseDrivers = new QPushButton("Browse...");
        driversLayout->addWidget(browseDrivers, 0, 2);

        layout->addWidget(driversBox);

        // --- Modeling Method ---
        auto* methodBox = new QGroupBox("Modeling Method");
        auto* methodLayout = new QHBoxLayout(methodBox);

        methodLayout->addWidget(new QLabel("Method:"));
        m_methodCombo = new QComboBox;
        m_methodCombo->addItems({
            "Logistic Regression",
            "Random Forest (DecisionForest)",
            "MLP Neural Network",
            "Weighted Normalized Likelihood (WNL)",
            "Support Vector Machine (SVM)",
            "SimWeight (KNN)"
        });
        methodLayout->addWidget(m_methodCombo, 1);

        layout->addWidget(methodBox);

        // --- Output ---
        auto* outputBox = new QGroupBox("Output");
        auto* outputLayout = new QGridLayout(outputBox);

        outputLayout->addWidget(new QLabel("Output transition potential:"), 0, 0);
        m_outputEdit = new QLineEdit;
        outputLayout->addWidget(m_outputEdit, 0, 1);
        auto* browseOutput = new QPushButton("Browse...");
        outputLayout->addWidget(browseOutput, 0, 2);

        layout->addWidget(outputBox);

        // --- Run ---
        m_runBtn = new QPushButton("Run Transition Sub-Model");
        layout->addWidget(m_runBtn);

        // Connections
        connect(browseDrivers, &QPushButton::clicked, this, [this]() {
            QStringList files = QFileDialog::getOpenFileNames(this, "Select Driver Rasters",
                QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
            if (!files.isEmpty())
                m_driversEdit->setText(files.join(","));
        });

        connect(browseOutput, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getSaveFileName(this, "Output Transition Potential",
                QString(), "GeoTIFF (*.tif);;All Files (*)");
            if (!path.isEmpty()) m_outputEdit->setText(path);
        });

        connect(m_runBtn, &QPushButton::clicked, this, &LcmTransitionPotentialsTab::runSubModel);
    }

    QString tabName() const override { return "Transition Potentials"; }
    QString tabDescription() const override {
        return "Model the potential for land cover transitions using driver variables.";
    }

    bool isComplete() const override {
        return m_session->contains("lcm/transition_potentials");
    }

    void onActivated() override {
        // Populate transitions table from change analysis results
        if (m_session->contains("lcm/change_stats_csv") && m_transitionsTable->rowCount() == 0) {
            populateTransitionsFromStats();
        }

        // Default output path
        if (m_outputEdit->text().isEmpty() && m_session->contains("lcm/output_dir")) {
            m_outputEdit->setText(
                m_session->get("lcm/output_dir").toString() + "/lcm_transition_potential.tif");
        }
    }

private slots:
    void runSubModel() {
        QString drivers = m_driversEdit->text().trimmed();
        QString output  = m_outputEdit->text().trimmed();

        if (drivers.isEmpty() || output.isEmpty()) {
            QMessageBox::warning(this, "Missing Input",
                "Driver variables and output path are required.");
            return;
        }

        // Map combo index to module name
        int methodIdx = m_methodCombo->currentIndex();
        QString moduleName;
        QMap<QString, QVariant> params;

        switch (methodIdx) {
            case 0: // Logistic Regression
                moduleName = "LOGISTIC_REG";
                params = {
                    {"change_raster", m_session->get("lcm/transition_map").toString()},
                    {"drivers", drivers},
                    {"output_potential", output},
                    {"max_iterations", 20}
                };
                break;
            case 1: // Random Forest — use TRANSITION_POTENTIAL with Random Forest method
                moduleName = "TRANSITION_POTENTIAL";
                params = {
                    {"transitions_file", m_session->get("lcm/change_stats_csv").toString()},
                    {"driver_variables", drivers},
                    {"method", 1}, // Random Forest
                    {"output", output}
                };
                break;
            default:
                // For other methods, use TRANSITION_POTENTIAL with Logistic Regression
                moduleName = "TRANSITION_POTENTIAL";
                params = {
                    {"transitions_file", m_session->get("lcm/change_stats_csv").toString()},
                    {"driver_variables", drivers},
                    {"method", 0},
                    {"output", output}
                };
                break;
        }

        m_runBtn->setEnabled(false);
        bool ok = runModule(moduleName, params);
        m_runBtn->setEnabled(true);

        if (ok) {
            m_session->set("lcm/transition_potentials", output);
            emit completionChanged();
        } else {
            QMessageBox::critical(this, "Sub-Model Failed",
                "Transition potential modeling failed. Check the status bar.");
        }
    }

private:
    void populateTransitionsFromStats() {
        // Parse the change stats CSV to extract class pairs with non-zero transitions
        QString csvPath = m_session->get("lcm/change_stats_csv").toString();
        QFile file(csvPath);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return;

        QTextStream in(&file);
        // Skip to transition matrix section (after blank line)
        bool inMatrix = false;
        QStringList classes;
        int row = 0;

        while (!in.atEnd()) {
            QString line = in.readLine().trimmed();

            if (line.startsWith("Transition Matrix")) {
                // Parse header to get class IDs
                QStringList parts = line.split(',');
                for (int i = 1; i < parts.size(); ++i)
                    classes.append(parts[i].trimmed());
                inMatrix = true;
                continue;
            }

            if (inMatrix && !line.isEmpty()) {
                QStringList parts = line.split(',');
                if (parts.isEmpty()) continue;
                QString fromClass = parts[0].trimmed();

                for (int j = 1; j < parts.size() && j - 1 < classes.size(); ++j) {
                    int count = parts[j].trimmed().toInt();
                    if (count > 0 && fromClass != classes[j - 1]) {
                        int r = m_transitionsTable->rowCount();
                        m_transitionsTable->insertRow(r);
                        m_transitionsTable->setItem(r, 0, new QTableWidgetItem(fromClass));
                        m_transitionsTable->setItem(r, 1, new QTableWidgetItem(classes[j - 1]));
                        m_transitionsTable->setItem(r, 2, new QTableWidgetItem("Sub-Model 1"));
                    }
                }
            }
        }

        m_transitionsTable->resizeColumnsToContents();
    }

    QTableWidget* m_transitionsTable = nullptr;
    QLineEdit*    m_driversEdit      = nullptr;
    QComboBox*    m_methodCombo      = nullptr;
    QLineEdit*    m_outputEdit       = nullptr;
    QPushButton*  m_runBtn           = nullptr;
};

ModelerTab* createLcmTransitionPotentialsTab(ModelerSession* s, QWidget* p) {
    return new LcmTransitionPotentialsTab(s, p);
}

} // namespace aplaceholder

#include "LcmTransitionPotentialsTab.moc"
