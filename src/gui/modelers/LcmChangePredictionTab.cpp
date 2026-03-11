#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>
#include <QSpinBox>
#include <QFileDialog>
#include <QMessageBox>
#include <QTableWidget>
#include <QHeaderView>

namespace aplaceholder {

class LcmChangePredictionTab : public ModelerTab {
    Q_OBJECT
public:
    explicit LcmChangePredictionTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Change Demand Modeling (Markov Chain) ---
        auto* demandBox = new QGroupBox("Change Demand Modeling (Markov Chain)");
        auto* demandLayout = new QGridLayout(demandBox);

        demandLayout->addWidget(new QLabel("Prediction date (years from later map):"), 0, 0);
        m_predictionYears = new QSpinBox;
        m_predictionYears->setRange(1, 500);
        m_predictionYears->setValue(10);
        demandLayout->addWidget(m_predictionYears, 0, 1);

        demandLayout->addWidget(new QLabel("Output matrix file:"), 1, 0);
        m_matrixEdit = new QLineEdit;
        demandLayout->addWidget(m_matrixEdit, 1, 1);
        auto* browseMatrix = new QPushButton("Browse...");
        demandLayout->addWidget(browseMatrix, 1, 2);

        m_runMarkovBtn = new QPushButton("Run Markov Chain");
        demandLayout->addWidget(m_runMarkovBtn, 2, 0, 1, 3);

        layout->addWidget(demandBox);

        // --- Change Allocation ---
        auto* allocBox = new QGroupBox("Change Allocation");
        auto* allocLayout = new QGridLayout(allocBox);

        allocLayout->addWidget(new QLabel("Prediction type:"), 0, 0);
        m_predTypeCombo = new QComboBox;
        m_predTypeCombo->addItems({"Hard (specific land cover)", "Soft (vulnerability map)"});
        allocLayout->addWidget(m_predTypeCombo, 0, 1);

        allocLayout->addWidget(new QLabel("Method:"), 1, 0);
        m_allocMethodCombo = new QComboBox;
        m_allocMethodCombo->addItems({"MOLA (Change Allocation)", "CA-Markov (Cellular Automata)"});
        allocLayout->addWidget(m_allocMethodCombo, 1, 1);

        allocLayout->addWidget(new QLabel("Recalculation stages:"), 2, 0);
        m_stagesSpin = new QSpinBox;
        m_stagesSpin->setRange(1, 100);
        m_stagesSpin->setValue(1);
        m_stagesSpin->setToolTip("Stage=1: allocate all change at once. "
                                  "Stage=N: divide change linearly across N stages.");
        allocLayout->addWidget(m_stagesSpin, 2, 1);

        allocLayout->addWidget(new QLabel("Output prediction map:"), 3, 0);
        m_predOutputEdit = new QLineEdit;
        allocLayout->addWidget(m_predOutputEdit, 3, 1);
        auto* browsePred = new QPushButton("Browse...");
        allocLayout->addWidget(browsePred, 3, 2);

        m_runPredBtn = new QPushButton("Run Change Prediction");
        m_runPredBtn->setEnabled(false);
        allocLayout->addWidget(m_runPredBtn, 4, 0, 1, 3);

        layout->addWidget(allocBox);

        // --- Validation ---
        auto* validBox = new QGroupBox("Validation");
        auto* validLayout = new QGridLayout(validBox);

        validLayout->addWidget(new QLabel("Reality map (actual later date):"), 0, 0);
        m_realityEdit = new QLineEdit;
        validLayout->addWidget(m_realityEdit, 0, 1);
        auto* browseReality = new QPushButton("Browse...");
        validLayout->addWidget(browseReality, 0, 2);

        m_validateBtn = new QPushButton("Validate Prediction");
        m_validateBtn->setEnabled(false);
        validLayout->addWidget(m_validateBtn, 1, 0, 1, 3);

        m_validationLabel = new QLabel;
        validLayout->addWidget(m_validationLabel, 2, 0, 1, 3);

        layout->addWidget(validBox);

        layout->addStretch();

        // Connections
        connect(browseMatrix, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getSaveFileName(this, "Output Matrix File",
                QString(), "CSV (*.csv);;All Files (*)");
            if (!path.isEmpty()) m_matrixEdit->setText(path);
        });
        connect(browsePred, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getSaveFileName(this, "Output Prediction Map",
                QString(), "GeoTIFF (*.tif);;All Files (*)");
            if (!path.isEmpty()) m_predOutputEdit->setText(path);
        });
        connect(browseReality, &QPushButton::clicked, this, [this]() {
            QString path = QFileDialog::getOpenFileName(this, "Select Reality Map",
                QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
            if (!path.isEmpty()) m_realityEdit->setText(path);
        });

        connect(m_runMarkovBtn, &QPushButton::clicked, this, &LcmChangePredictionTab::runMarkov);
        connect(m_runPredBtn,   &QPushButton::clicked, this, &LcmChangePredictionTab::runPrediction);
        connect(m_validateBtn,  &QPushButton::clicked, this, &LcmChangePredictionTab::runValidation);
    }

    QString tabName() const override { return "Change Prediction"; }
    QString tabDescription() const override {
        return "Predict future land cover change using Markov chain and change allocation.";
    }

    bool isComplete() const override {
        return m_session->contains("lcm/prediction_map");
    }

    void onActivated() override {
        if (m_matrixEdit->text().isEmpty() && m_session->contains("lcm/output_dir")) {
            QString dir = m_session->get("lcm/output_dir").toString();
            m_matrixEdit->setText(dir + "/lcm_markov_matrix.csv");
            m_predOutputEdit->setText(dir + "/lcm_prediction.tif");
        }
    }

private slots:
    void runMarkov() {
        QString matrixPath = m_matrixEdit->text().trimmed();
        if (matrixPath.isEmpty()) {
            QMessageBox::warning(this, "Missing Input", "Output matrix path is required.");
            return;
        }

        m_runMarkovBtn->setEnabled(false);
        bool ok = runModule("MARKOV_CHAIN", {
            {"earlier_image", m_session->get("lcm/earlier_map").toString()},
            {"later_image",   m_session->get("lcm/later_map").toString()},
            {"output_matrix", matrixPath},
            {"prediction_date", m_predictionYears->value()}
        });
        m_runMarkovBtn->setEnabled(true);

        if (ok) {
            m_session->set("lcm/markov_matrix", matrixPath);
            m_runPredBtn->setEnabled(true);
            emit statusMessage("Markov chain completed.");
        } else {
            QMessageBox::critical(this, "Markov Chain Failed",
                "Markov chain analysis failed.");
        }
    }

    void runPrediction() {
        QString output = m_predOutputEdit->text().trimmed();
        if (output.isEmpty()) {
            QMessageBox::warning(this, "Missing Input", "Output prediction path is required.");
            return;
        }

        bool useCA = (m_allocMethodCombo->currentIndex() == 1);
        QString moduleName = useCA ? "CELLULAR_AUTOMATA" : "CHGALLOC";

        m_runPredBtn->setEnabled(false);
        bool ok = false;

        if (useCA) {
            ok = runModule("CELLULAR_AUTOMATA", {
                {"base_image", m_session->get("lcm/later_map").toString()},
                {"transition_matrix", m_session->get("lcm/markov_matrix").toString()},
                {"suitability_maps", m_session->get("lcm/transition_potentials").toString()},
                {"iterations", m_stagesSpin->value()},
                {"output", output}
            });
        } else {
            ok = runModule("CHGALLOC", {
                {"base_landcover", m_session->get("lcm/later_map").toString()},
                {"transition_potentials", m_session->get("lcm/transition_potentials").toString()},
                {"change_demand_file", m_session->get("lcm/markov_matrix").toString()},
                {"output", output}
            });
        }

        m_runPredBtn->setEnabled(true);

        if (ok) {
            m_session->set("lcm/prediction_map", output);
            m_session->set("lcm/prediction_type", m_predTypeCombo->currentIndex());
            m_validateBtn->setEnabled(true);
            emit completionChanged();
            emit statusMessage("Change prediction completed.");
        } else {
            QMessageBox::critical(this, "Prediction Failed",
                "Change prediction failed.");
        }
    }

    void runValidation() {
        QString reality = m_realityEdit->text().trimmed();
        if (reality.isEmpty()) {
            QMessageBox::warning(this, "Missing Input", "Reality map is required for validation.");
            return;
        }

        // Use CROSSTAB for 3-way validation
        QString validOutput = m_session->get("lcm/output_dir").toString() + "/lcm_validation.tif";
        bool ok = runModule("CROSSTAB", {
            {"input1", m_session->get("lcm/prediction_map").toString()},
            {"input2", reality},
            {"output", validOutput}
        });

        if (ok) {
            m_validationLabel->setText("Validation complete. Results saved to: " + validOutput);
        } else {
            m_validationLabel->setText("Validation failed.");
        }
    }

private:
    QSpinBox*    m_predictionYears  = nullptr;
    QLineEdit*   m_matrixEdit       = nullptr;
    QPushButton* m_runMarkovBtn     = nullptr;
    QComboBox*   m_predTypeCombo    = nullptr;
    QComboBox*   m_allocMethodCombo = nullptr;
    QSpinBox*    m_stagesSpin       = nullptr;
    QLineEdit*   m_predOutputEdit   = nullptr;
    QPushButton* m_runPredBtn       = nullptr;
    QLineEdit*   m_realityEdit      = nullptr;
    QPushButton* m_validateBtn      = nullptr;
    QLabel*      m_validationLabel  = nullptr;
};

ModelerTab* createLcmChangePredictionTab(ModelerSession* s, QWidget* p) {
    return new LcmChangePredictionTab(s, p);
}

} // namespace aplaceholder

#include "LcmChangePredictionTab.moc"
