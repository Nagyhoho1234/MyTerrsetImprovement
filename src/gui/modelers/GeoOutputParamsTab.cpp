#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QTextEdit>
#include <QFileDialog>

namespace aplaceholder {

class GeoOutputParamsTab : public ModelerTab {
    Q_OBJECT
public:
    explicit GeoOutputParamsTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        auto* box = new QGroupBox("Output Parameters");
        auto* gLayout = new QGridLayout(box);

        gLayout->addWidget(new QLabel(
            "Calculate proportional change in agricultural price via\n"
            "iterative convergence loop."), 0, 0, 1, 2);

        gLayout->addWidget(new QLabel("Max iterations:"), 1, 0);
        m_maxIter = new QSpinBox;
        m_maxIter->setRange(10, 10000);
        m_maxIter->setValue(1000);
        gLayout->addWidget(m_maxIter, 1, 1);

        gLayout->addWidget(new QLabel("Model precision:"), 2, 0);
        m_precisionEdit = new QLineEdit("0.0001");
        gLayout->addWidget(m_precisionEdit, 2, 1);

        gLayout->addWidget(new QLabel("Output directory:"), 3, 0);
        m_outputDir = new QLineEdit;
        gLayout->addWidget(m_outputDir, 3, 1);
        auto* browseDir = new QPushButton("Browse...");
        gLayout->addWidget(browseDir, 3, 2);

        auto* runBtn = new QPushButton("Run GeOSIRIS Model");
        gLayout->addWidget(runBtn, 4, 0, 1, 3);

        gLayout->addWidget(new QLabel("Output summary:"), 5, 0);
        m_summary = new QTextEdit;
        m_summary->setReadOnly(true);
        m_summary->setMaximumHeight(200);
        m_summary->setPlaceholderText(
            "Output images:\n"
            "- Carbon emissions with/without REDD\n"
            "- Change in emissions due to REDD\n"
            "- Deforestation with/without REDD\n"
            "- Effective opportunity cost\n"
            "- Emission factor, Emittable forest carbon\n"
            "- Site level baseline\n"
            "- Would-be deforestation/emissions\n"
            "- District decisions (opt-in/opt-out)\n"
            "- Proportional change in ag price");
        gLayout->addWidget(m_summary, 5, 1, 1, 2);

        layout->addWidget(box);
        layout->addStretch();

        connect(browseDir, &QPushButton::clicked, this, [this]() {
            QString dir = QFileDialog::getExistingDirectory(this, "Select Output Directory");
            if (!dir.isEmpty()) m_outputDir->setText(dir);
        });

        connect(runBtn, &QPushButton::clicked, this, [this]() {
            m_session->set("geo/max_iterations", m_maxIter->value());
            m_session->set("geo/precision", m_precisionEdit->text().toDouble());
            m_session->set("geo/output_dir", m_outputDir->text().trimmed());
            m_summary->setText("GeOSIRIS model configured.\n"
                "Max iterations: " + QString::number(m_maxIter->value()) +
                "\nPrecision: " + m_precisionEdit->text() +
                "\n\nReady to execute when all input data is available.");
            emit statusMessage("GeOSIRIS output parameters configured.");
        });
    }

    QString tabName() const override { return "Output Parameters"; }
    QString tabDescription() const override { return "Configure and run the GeOSIRIS iterative model."; }
    bool isComplete() const override { return true; } // final tab

private:
    QSpinBox* m_maxIter = nullptr;
    QLineEdit* m_precisionEdit = nullptr;
    QLineEdit* m_outputDir = nullptr;
    QTextEdit* m_summary = nullptr;
};

ModelerTab* createGeoOutputParamsTab(ModelerSession* s, QWidget* p) {
    return new GeoOutputParamsTab(s, p);
}

} // namespace aplaceholder

#include "GeoOutputParamsTab.moc"
