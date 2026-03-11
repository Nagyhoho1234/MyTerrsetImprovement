#pragma once

#include <QDialog>
#include <QFormLayout>
#include <QMap>
#include <QWidget>

namespace aplaceholder {

class Module;

// Auto-generates a dialog from a module's parameter definitions
class ModuleDialog : public QDialog {
    Q_OBJECT

public:
    ModuleDialog(Module* module, QWidget* parent = nullptr);

private slots:
    void onAccept();

private:
    Module* m_module;
    QMap<QString, QWidget*> m_widgets;
};

} // namespace aplaceholder
