#include "modelers/HbmWizard.h"
#include "ModelerTab.h"

namespace aplaceholder {
ModelerTab* createHbmSpeciesTab(ModelerSession* s, QWidget* p);
ModelerTab* createHbmBiodiversityTab(ModelerSession* s, QWidget* p);
ModelerTab* createHbmLandscapeTab(ModelerSession* s, QWidget* p);
ModelerTab* createHbmPlanningTab(ModelerSession* s, QWidget* p);
}

namespace aplaceholder {

HbmWizard::HbmWizard(QWidget* parent)
    : ModelerWizard("Habitat and Biodiversity Modeler", parent)
{
    setupTabs();
}

void HbmWizard::setupTabs()
{
    addTab(createHbmSpeciesTab(&m_session, this));
    addTab(createHbmBiodiversityTab(&m_session, this));
    addTab(createHbmLandscapeTab(&m_session, this));
    addTab(createHbmPlanningTab(&m_session, this));
}

} // namespace aplaceholder
