#pragma once

#include "ModelerWizard.h"

namespace aplaceholder {

class HbmWizard : public ModelerWizard {
    Q_OBJECT
public:
    explicit HbmWizard(QWidget* parent = nullptr);

private:
    void setupTabs();
};

} // namespace aplaceholder
