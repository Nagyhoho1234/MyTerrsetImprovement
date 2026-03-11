#include "ModuleDialog.h"
#include "Module.h"
#include "Parameter.h"
#include <QFormLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QFileDialog>
#include <QLabel>
#include <QGroupBox>
#include <QDialogButtonBox>
#include <QToolButton>

namespace aplaceholder {

ModuleDialog::ModuleDialog(Module* module, QWidget* parent)
    : QDialog(parent), m_module(module)
{
    setWindowTitle(module->name());
    setMinimumWidth(500);

    auto* mainLayout = new QVBoxLayout(this);

    // Description
    auto* descLabel = new QLabel(module->description());
    descLabel->setWordWrap(true);
    descLabel->setStyleSheet("color: #666; margin-bottom: 10px;");
    mainLayout->addWidget(descLabel);

    // Parameter form
    auto* form = new QFormLayout();
    form->setLabelAlignment(Qt::AlignRight);

    for (const auto& pdef : module->parameterDefs()) {
        QWidget* widget = nullptr;

        switch (pdef.type) {
        case ParameterDef::File:
        case ParameterDef::OutputFile: {
            auto* row = new QWidget();
            auto* rowLayout = new QHBoxLayout(row);
            rowLayout->setContentsMargins(0, 0, 0, 0);
            auto* edit = new QLineEdit();
            auto* btn = new QToolButton();
            btn->setText("...");
            rowLayout->addWidget(edit);
            rowLayout->addWidget(btn);

            bool isOutput = (pdef.type == ParameterDef::OutputFile);
            connect(btn, &QToolButton::clicked, [edit, isOutput, this]() {
                QString path = isOutput
                    ? QFileDialog::getSaveFileName(this, "Output File")
                    : QFileDialog::getOpenFileName(this, "Input File");
                if (!path.isEmpty())
                    edit->setText(path);
            });
            widget = row;
            m_widgets[pdef.key] = edit;
            break;
        }
        case ParameterDef::Integer: {
            auto* spin = new QSpinBox();
            spin->setRange(static_cast<int>(pdef.minValue), static_cast<int>(pdef.maxValue));
            spin->setValue(pdef.defaultValue.toInt());
            widget = spin;
            m_widgets[pdef.key] = spin;
            break;
        }
        case ParameterDef::Double: {
            auto* spin = new QDoubleSpinBox();
            spin->setRange(pdef.minValue, pdef.maxValue);
            spin->setValue(pdef.defaultValue.toDouble());
            spin->setDecimals(6);
            widget = spin;
            m_widgets[pdef.key] = spin;
            break;
        }
        case ParameterDef::String: {
            auto* edit = new QLineEdit(pdef.defaultValue.toString());
            widget = edit;
            m_widgets[pdef.key] = edit;
            break;
        }
        case ParameterDef::Combo: {
            auto* combo = new QComboBox();
            combo->addItems(pdef.options);
            combo->setCurrentIndex(pdef.defaultValue.toInt());
            widget = combo;
            m_widgets[pdef.key] = combo;
            break;
        }
        case ParameterDef::Bool: {
            auto* check = new QCheckBox();
            check->setChecked(pdef.defaultValue.toBool());
            widget = check;
            m_widgets[pdef.key] = check;
            break;
        }
        case ParameterDef::RasterLayer:
        case ParameterDef::MultiRaster:
        case ParameterDef::BandSelect: {
            auto* edit = new QLineEdit();
            edit->setPlaceholderText("Select raster layer...");
            widget = edit;
            m_widgets[pdef.key] = edit;
            break;
        }
        case ParameterDef::Separator: {
            auto* line = new QLabel(pdef.label);
            line->setStyleSheet("font-weight: bold; margin-top: 8px;");
            form->addRow(line);
            continue;
        }
        }

        if (widget) {
            if (!pdef.tooltip.isEmpty())
                widget->setToolTip(pdef.tooltip);
            form->addRow(pdef.label + ":", widget);
        }
    }

    mainLayout->addLayout(form);

    // Buttons
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttons, &QDialogButtonBox::accepted, this, &ModuleDialog::onAccept);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    mainLayout->addWidget(buttons);
}

void ModuleDialog::onAccept()
{
    // Transfer widget values to module parameters
    for (auto it = m_widgets.begin(); it != m_widgets.end(); ++it) {
        QVariant val;
        if (auto* edit = qobject_cast<QLineEdit*>(it.value()))
            val = edit->text();
        else if (auto* spin = qobject_cast<QSpinBox*>(it.value()))
            val = spin->value();
        else if (auto* dspin = qobject_cast<QDoubleSpinBox*>(it.value()))
            val = dspin->value();
        else if (auto* combo = qobject_cast<QComboBox*>(it.value()))
            val = combo->currentIndex();
        else if (auto* check = qobject_cast<QCheckBox*>(it.value()))
            val = check->isChecked();

        m_module->setParameter(it.key(), val);
    }
    accept();
}

} // namespace aplaceholder
