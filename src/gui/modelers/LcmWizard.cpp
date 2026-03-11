#include "modelers/LcmWizard.h"
#include "ModelerTab.h"

// Factory functions implemented in each tab's .cpp
namespace aplaceholder {
ModelerTab* createLcmChangeAnalysisTab(ModelerSession* s, QWidget* p);
ModelerTab* createLcmTransitionPotentialsTab(ModelerSession* s, QWidget* p);
ModelerTab* createLcmChangePredictionTab(ModelerSession* s, QWidget* p);
ModelerTab* createLcmPlanningTab(ModelerSession* s, QWidget* p);
ModelerTab* createLcmReddTab(ModelerSession* s, QWidget* p);
ModelerTab* createLcmHarmonizeTab(ModelerSession* s, QWidget* p);
}

namespace aplaceholder {

LcmWizard::LcmWizard(QWidget* parent)
    : ModelerWizard("Land Change Modeler", parent)
{
    setupTabs();
}

void LcmWizard::setupTabs()
{
    addTab(createLcmChangeAnalysisTab(&m_session, this));
    addTab(createLcmTransitionPotentialsTab(&m_session, this));
    addTab(createLcmChangePredictionTab(&m_session, this));
    addTab(createLcmPlanningTab(&m_session, this));
    addTab(createLcmReddTab(&m_session, this));
    addTab(createLcmHarmonizeTab(&m_session, this));
}

} // namespace aplaceholder
