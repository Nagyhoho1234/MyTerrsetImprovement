#pragma once

#include "ModelerWizard.h"

namespace aplaceholder {

class EtmWizard : public ModelerWizard {
    Q_OBJECT
public:
    explicit EtmWizard(QWidget* parent = nullptr);

private:
    void setupTabs();
};

} // namespace aplaceholder
