#include "modelers/EtmWizard.h"
#include "ModelerTab.h"

namespace aplaceholder {
ModelerTab* createEtmExploreTab(ModelerSession* s, QWidget* p);
ModelerTab* createEtmAnalysisTab(ModelerSession* s, QWidget* p);
ModelerTab* createEtmPreprocessTab(ModelerSession* s, QWidget* p);
}

namespace aplaceholder {

EtmWizard::EtmWizard(QWidget* parent)
    : ModelerWizard("Earth Trends Modeler", parent)
{
    setupTabs();
}

void EtmWizard::setupTabs()
{
    addTab(createEtmExploreTab(&m_session, this));
    addTab(createEtmAnalysisTab(&m_session, this));
    addTab(createEtmPreprocessTab(&m_session, this));
}

} // namespace aplaceholder
