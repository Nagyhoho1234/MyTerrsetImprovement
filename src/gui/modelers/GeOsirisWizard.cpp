#include "modelers/GeOsirisWizard.h"
#include "ModelerTab.h"

namespace aplaceholder {
ModelerTab* createGeoProjectParamsTab(ModelerSession* s, QWidget* p);
ModelerTab* createGeoExternalFactorsTab(ModelerSession* s, QWidget* p);
ModelerTab* createGeoReddRulesTab(ModelerSession* s, QWidget* p);
ModelerTab* createGeoModelParamsTab(ModelerSession* s, QWidget* p);
ModelerTab* createGeoInputImagesTab(ModelerSession* s, QWidget* p);
ModelerTab* createGeoOpportunityCostTab(ModelerSession* s, QWidget* p);
ModelerTab* createGeoOutputParamsTab(ModelerSession* s, QWidget* p);
}

namespace aplaceholder {

GeOsirisWizard::GeOsirisWizard(QWidget* parent)
    : ModelerWizard("GeOSIRIS", parent)
{
    setupTabs();
}

void GeOsirisWizard::setupTabs()
{
    addTab(createGeoProjectParamsTab(&m_session, this));
    addTab(createGeoExternalFactorsTab(&m_session, this));
    addTab(createGeoReddRulesTab(&m_session, this));
    addTab(createGeoModelParamsTab(&m_session, this));
    addTab(createGeoInputImagesTab(&m_session, this));
    addTab(createGeoOpportunityCostTab(&m_session, this));
    addTab(createGeoOutputParamsTab(&m_session, this));
}

} // namespace aplaceholder
