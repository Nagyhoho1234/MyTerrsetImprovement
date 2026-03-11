#include "ModelerTab.h"
#include "ModelerSession.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>
#include <QFileDialog>
#include <QMessageBox>

namespace aplaceholder {

class GeoInputImagesTab : public ModelerTab {
    Q_OBJECT
public:
    explicit GeoInputImagesTab(ModelerSession* session, QWidget* parent = nullptr)
        : ModelerTab(session, parent)
    {
        auto* layout = new QVBoxLayout(this);

        auto* box = new QGroupBox("Input Image Files");
        auto* gLayout = new QGridLayout(box);

        gLayout->addWidget(new QLabel("Land cover type:"), 0, 0);
        m_lcType = new QComboBox;
        m_lcType->addItems({"Boolean", "Fractional"});
        gLayout->addWidget(m_lcType, 0, 1);

        struct InputRow { const char* label; int row; };
        InputRow rows[] = {
            {"Reference/Land Area:", 1},
            {"Forest Before REDD+:", 2},
            {"Deforested Area:", 3},
            {"Biomass Carbon (tC/ha):", 4},
            {"Soil Carbon (tC/ha):", 5},
            {"Peat (Boolean presence):", 6},
            {"Administrative Levels:", 7},
            {"Agricultural Revenue ($/ha NPV):", 8},
        };

        for (auto& r : rows) {
            gLayout->addWidget(new QLabel(r.label), r.row, 0);
            auto* edit = new QLineEdit;
            gLayout->addWidget(edit, r.row, 1);
            auto* btn = new QPushButton("Browse...");
            gLayout->addWidget(btn, r.row, 2);
            m_edits.push_back(edit);
            connect(btn, &QPushButton::clicked, this, [this, edit]() {
                QString path = QFileDialog::getOpenFileName(this, "Select Input Image",
                    QString(), "Raster Files (*.tif *.tiff *.rst *.img);;All Files (*)");
                if (!path.isEmpty()) edit->setText(path);
            });
        }

        auto* applyBtn = new QPushButton("Apply Input Images");
        gLayout->addWidget(applyBtn, 9, 0, 1, 3);
        layout->addWidget(box);

        layout->addStretch();

        connect(applyBtn, &QPushButton::clicked, this, [this]() {
            QStringList keys = {"reference_area", "forest_before", "deforested_area",
                               "biomass_carbon", "soil_carbon", "peat",
                               "admin_levels", "ag_revenue"};
            for (int i = 0; i < keys.size() && i < (int)m_edits.size(); ++i) {
                m_session->set("geo/" + keys[i], m_edits[i]->text().trimmed());
            }
            m_session->set("geo/lc_type", m_lcType->currentText());
            m_applied = true;
            emit statusMessage("Input images configured.");
            emit completionChanged();
        });
    }

    QString tabName() const override { return "Input Images"; }
    QString tabDescription() const override { return "Specify input raster images for GeOSIRIS analysis."; }
    bool isComplete() const override { return m_applied; }

private:
    QComboBox* m_lcType = nullptr;
    std::vector<QLineEdit*> m_edits;
    bool m_applied = false;
};

ModelerTab* createGeoInputImagesTab(ModelerSession* s, QWidget* p) {
    return new GeoInputImagesTab(s, p);
}

} // namespace aplaceholder

#include "GeoInputImagesTab.moc"
