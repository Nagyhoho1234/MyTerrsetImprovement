#pragma once

#include "ModelerWizard.h"

namespace aplaceholder {

class LcmWizard : public ModelerWizard {
    Q_OBJECT
public:
    explicit LcmWizard(QWidget* parent = nullptr);

private:
    void setupTabs();
};

} // namespace aplaceholder
