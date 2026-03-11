#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QFileDialog>

namespace aplaceholder {

class GeoProjectParamsTab : public ModelerTab {
    Q_OBJECT
public:
    explicit GeoProjectParamsTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        auto* box = new QGroupBox("GeOSIRIS Project Parameters");
        auto* gLayout = new QGridLayout(box);

        gLayout->addWidget(new QLabel("Project name:"), 0, 0);
        m_projectName = new QLineEdit;
        gLayout->addWidget(m_projectName, 0, 1);

        gLayout->addWidget(new QLabel("Project file:"), 1, 0);
        m_projectFile = new QLineEdit;
        gLayout->addWidget(m_projectFile, 1, 1);
        auto* browseProj = new QPushButton("Browse...");
        gLayout->addWidget(browseProj, 1, 2);

        gLayout->addWidget(new QLabel("Vector overlay (optional):"), 2, 0);
        m_overlayEdit = new QLineEdit;
        gLayout->addWidget(m_overlayEdit, 2, 1);
        auto* browseOverlay = new QPushButton("Browse...");
        gLayout->addWidget(browseOverlay, 2, 2);

        gLayout->addWidget(new QLabel("Custom palette (optional):"), 3, 0);
        m_paletteEdit = new QLineEdit;
        gLayout->addWidget(m_paletteEdit, 3, 1);
        auto* browsePalette = new QPushButton("Browse...");
        gLayout->addWidget(browsePalette, 3, 2);

        auto* saveBtn = new QPushButton("Save Project Parameters");
        gLayout->addWidget(saveBtn, 4, 0, 1, 3);
        layout->addWidget(box);

        auto* loadBtn = new QPushButton("Load Existing Project...");
        layout->addWidget(loadBtn);

        layout->addStretch();

        auto browse = [this](QLineEdit* edit, const QString& title, const QString& filter) {
            QString path = QFileDialog::getOpenFileName(this, title, QString(), filter);
            if (!path.isEmpty()) edit->setText(path);
        };
        connect(browseProj, &QPushButton::clicked, this, [=]() {
            browse(m_projectFile, "Select Project File", "GeOSIRIS Projects (*.gos);;All Files (*)");
        });
        connect(browseOverlay, &QPushButton::clicked, this, [=]() {
            browse(m_overlayEdit, "Select Vector Overlay", "Vector Files (*.vct *.shp);;All Files (*)");
        });
        connect(browsePalette, &QPushButton::clicked, this, [=]() {
            browse(m_paletteEdit, "Select Palette", "Palette Files (*.smp *.pal);;All Files (*)");
        });

        connect(saveBtn, &QPushButton::clicked, this, [this]() {
            m_session->set("geo/project_name", m_projectName->text().trimmed());
            m_session->set("geo/project_file", m_projectFile->text().trimmed());
            m_configured = true;
            emit statusMessage("Project parameters saved.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "Project Parameters"; }
    QString tabDescription() const override { return "Create or load a GeOSIRIS project."; }
    bool isComplete() const override { return m_configured; }

private:
    QLineEdit* m_projectName = nullptr;
    QLineEdit* m_projectFile = nullptr;
    QLineEdit* m_overlayEdit = nullptr;
    QLineEdit* m_paletteEdit = nullptr;
    bool m_configured = false;
};

ModelerTab* createGeoProjectParamsTab(ModelerSession* s, QWidget* p) {
    return new GeoProjectParamsTab(s, p);
}

} // namespace aplaceholder

#include "GeoProjectParamsTab.moc"
