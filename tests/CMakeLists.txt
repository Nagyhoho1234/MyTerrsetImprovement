find_package(Qt6 REQUIRED COMPONENTS Test)

add_executable(test_raster test_raster.cpp)
target_link_libraries(test_raster PRIVATE aplaceholder_core Qt6::Test)
add_test(NAME RasterTest COMMAND test_raster)

add_executable(test_io test_io.cpp)
target_link_libraries(test_io PRIVATE aplaceholder_io Qt6::Test)
add_test(NAME IoTest COMMAND test_io)

add_executable(smoke_test smoke_test.cpp)
target_link_libraries(smoke_test PRIVATE
    aplaceholder_core
    aplaceholder_io
    aplaceholder_gis
    aplaceholder_imgproc
    aplaceholder_lcm
    aplaceholder_etm
    aplaceholder_hbm
    aplaceholder_esm
    aplaceholder_ccam
    Qt6::Core
    GDAL::GDAL
)
# Force-link all module static libs so REGISTER_MODULE global constructors run
set(_module_libs aplaceholder_gis aplaceholder_imgproc aplaceholder_lcm
    aplaceholder_etm aplaceholder_hbm aplaceholder_esm aplaceholder_ccam)
foreach(_lib ${_module_libs})
    target_link_options(smoke_test PRIVATE "/WHOLEARCHIVE:$<TARGET_FILE:${_lib}>")
endforeach()
add_test(NAME SmokeTest COMMAND smoke_test)

add_executable(comprehensive_test comprehensive_test.cpp)
target_link_libraries(comprehensive_test PRIVATE
    aplaceholder_core
    aplaceholder_io
    aplaceholder_gis
    aplaceholder_imgproc
    aplaceholder_lcm
    aplaceholder_etm
    aplaceholder_hbm
    aplaceholder_esm
    aplaceholder_ccam
    Qt6::Core
    GDAL::GDAL
)
foreach(_lib ${_module_libs})
    target_link_options(comprehensive_test PRIVATE "/WHOLEARCHIVE:$<TARGET_FILE:${_lib}>")
endforeach()
add_test(NAME ComprehensiveTest COMMAND comprehensive_test)

set(_wizard_sources
    ${CMAKE_SOURCE_DIR}/src/gui/ModelerWizard.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/ModelerTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/LcmWizard.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/LcmChangeAnalysisTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/LcmTransitionPotentialsTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/LcmChangePredictionTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/LcmPlanningTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/LcmReddTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/LcmHarmonizeTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/HbmWizard.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/HbmSpeciesTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/HbmBiodiversityTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/HbmLandscapeTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/HbmPlanningTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/GeOsirisWizard.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/GeoProjectParamsTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/GeoExternalFactorsTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/GeoReddRulesTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/GeoModelParamsTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/GeoInputImagesTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/GeoOpportunityCostTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/GeoOutputParamsTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/EtmWizard.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/EtmExploreTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/EtmAnalysisTab.cpp
    ${CMAKE_SOURCE_DIR}/src/gui/modelers/EtmPreprocessTab.cpp
    ${CMAKE_SOURCE_DIR}/include/ModelerWizard.h
    ${CMAKE_SOURCE_DIR}/include/ModelerTab.h
    ${CMAKE_SOURCE_DIR}/include/ModelerSession.h
    ${CMAKE_SOURCE_DIR}/include/modelers/LcmWizard.h
    ${CMAKE_SOURCE_DIR}/include/modelers/HbmWizard.h
    ${CMAKE_SOURCE_DIR}/include/modelers/GeOsirisWizard.h
    ${CMAKE_SOURCE_DIR}/include/modelers/EtmWizard.h
)
add_executable(test_wizards test_wizards.cpp ${_wizard_sources})
target_link_libraries(test_wizards PRIVATE
    aplaceholder_core
    aplaceholder_io
    Qt6::Widgets
    Qt6::Gui
    GDAL::GDAL
)
add_test(NAME WizardTest COMMAND test_wizards)
