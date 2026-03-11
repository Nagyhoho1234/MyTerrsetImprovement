#pragma once

#include "ModelerWizard.h"

namespace aplaceholder {

class GeOsirisWizard : public ModelerWizard {
    Q_OBJECT
public:
    explicit GeOsirisWizard(QWidget* parent = nullptr);

private:
    void setupTabs();
};

} // namespace aplaceholder
