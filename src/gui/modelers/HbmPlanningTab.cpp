#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QMessageBox>

namespace aplaceholder {

class HbmPlanningTab : public ModelerTab {
    Q_OBJECT
public:
    explicit HbmPlanningTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        // --- Corridor Planning ---
        auto* corrBox = new QGroupBox("Corridor Planning");
        auto* cLayout = new QGridLayout(corrBox);

        cLayout->addWidget(new QLabel(
            "Build biological corridors between two terminal regions\n"
            "using cost-distance analysis and habitat suitability."), 0, 0, 1, 3);

        cLayout->addWidget(new QLabel("Terminal region A (Boolean):"), 1, 0);
        m_termAEdit = new QLineEdit;
        cLayout->addWidget(m_termAEdit, 1, 1);
        auto* browseA = new QPushButton("Browse...");
        cLayout->addWidget(browseA, 1, 2);

        cLayout->addWidget(new QLabel("Terminal region B (Boolean):"), 2, 0);
        m_termBEdit = new QLineEdit;
        cLayout->addWidget(m_termBEdit, 2, 1);
        auto* browseB = new QPushButton("Browse...");
        cLayout->addWidget(browseB, 2, 2);

        cLayout->addWidget(new QLabel("Habitat suitability map:"), 3, 0);
        m_suitEdit = new QLineEdit;
        cLayout->addWidget(m_suitEdit, 3, 1);
        auto* browseSuit = new QPushButton("Browse...");
        cLayout->addWidget(browseSuit, 3, 2);

        cLayout->addWidget(new QLabel("Development suitability (optional):"), 4, 0);
        m_devEdit = new QLineEdit;
        cLayout->addWidget(m_devEdit, 4, 1);
        auto* browseDev = new QPushButton("Browse...");
        cLayout->addWidget(browseDev, 4, 2);

        cLayout->addWidget(new QLabel("Conservation value (optional):"), 5, 0);
        m_consEdit = new QLineEdit;
        cLayout->addWidget(m_consEdit, 5, 1);
        auto* browseCons = new QPushButton("Browse...");
        cLayout->addWidget(browseCons, 5, 2);

        cLayout->addWidget(new QLabel("Protected lands (optional):"), 6, 0);
        m_protEdit = new QLineEdit;
        cLayout->addWidget(m_protEdit, 6, 1);
        auto* browseProt = new QPushButton("Browse...");
        cLayout->addWidget(browseProt, 6, 2);

        cLayout->addWidget(new QLabel("Ideal corridor width:"), 7, 0);
        m_corrWidth = new QDoubleSpinBox;
        m_corrWidth->setRange(1, 100000);
        m_corrWidth->setValue(1000);
        cLayout->addWidget(m_corrWidth, 7, 1);

        cLayout->addWidget(new QLabel("Number of branches:"), 8, 0);
        m_numBranches = new QSpinBox;
        m_numBranches->setRange(1, 20);
        m_numBranches->setValue(1);
        cLayout->addWidget(m_numBranches, 8, 1);

        auto* runCorridor = new QPushButton("Build Corridors");
        cLayout->addWidget(runCorridor, 9, 0, 1, 3);
        layout->addWidget(corrBox);

        layout->addStretch();

        // Browse connections
        auto browse = [this](QLineEdit* edit, const QString& title) {
            QString path = QFileDialog::getOpenFileName(this, title,
                QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
            if (!path.isEmpty()) edit->setText(path);
        };
        connect(browseA, &QPushButton::clicked, this, [=]() { browse(m_termAEdit, "Select Terminal Region A"); });
        connect(browseB, &QPushButton::clicked, this, [=]() { browse(m_termBEdit, "Select Terminal Region B"); });
        connect(browseSuit, &QPushButton::clicked, this, [=]() { browse(m_suitEdit, "Select Habitat Suitability"); });
        connect(browseDev, &QPushButton::clicked, this, [=]() { browse(m_devEdit, "Select Development Suitability"); });
        connect(browseCons, &QPushButton::clicked, this, [=]() { browse(m_consEdit, "Select Conservation Value"); });
        connect(browseProt, &QPushButton::clicked, this, [=]() { browse(m_protEdit, "Select Protected Lands"); });

        connect(runCorridor, &QPushButton::clicked, this, [this]() {
            if (m_termAEdit->text().trimmed().isEmpty() || m_termBEdit->text().trimmed().isEmpty()) {
                QMessageBox::warning(this, "Missing Input", "Please specify both terminal regions.");
                return;
            }
            if (m_suitEdit->text().trimmed().isEmpty()) {
                QMessageBox::warning(this, "Missing Input", "Please specify a habitat suitability map.");
                return;
            }
            m_session->set("hbm/terminal_a", m_termAEdit->text().trimmed());
            m_session->set("hbm/terminal_b", m_termBEdit->text().trimmed());
            m_session->set("hbm/corridor_suitability", m_suitEdit->text().trimmed());
            m_session->set("hbm/corridor_width", m_corrWidth->value());
            m_session->set("hbm/corridor_branches", m_numBranches->value());
            emit statusMessage("Corridor planning parameters saved.");
        });
    }

    QString tabName() const override { return "Planning"; }
    QString tabDescription() const override {
        return "Biological corridor planning between terminal regions.";
    }
    bool isComplete() const override { return true; } // optional

private:
    QLineEdit* m_termAEdit = nullptr;
    QLineEdit* m_termBEdit = nullptr;
    QLineEdit* m_suitEdit = nullptr;
    QLineEdit* m_devEdit = nullptr;
    QLineEdit* m_consEdit = nullptr;
    QLineEdit* m_protEdit = nullptr;
    QDoubleSpinBox* m_corrWidth = nullptr;
    QSpinBox* m_numBranches = nullptr;
};

ModelerTab* createHbmPlanningTab(ModelerSession* s, QWidget* p) {
    return new HbmPlanningTab(s, p);
}

} // namespace aplaceholder

#include "HbmPlanningTab.moc"
