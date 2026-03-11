add_executable(MyTerrsetImprovement
    main.cpp
    MainWindow.cpp
    MapCanvas.cpp
    ModuleDialog.cpp
    ChartWidget.cpp
    RasterInfoPanel.cpp
    ModelerTab.cpp
    ModelerWizard.cpp
    modelers/LcmWizard.cpp
    modelers/LcmChangeAnalysisTab.cpp
    modelers/LcmTransitionPotentialsTab.cpp
    modelers/LcmChangePredictionTab.cpp
    modelers/LcmPlanningTab.cpp
    modelers/LcmReddTab.cpp
    modelers/LcmHarmonizeTab.cpp
    modelers/HbmWizard.cpp
    modelers/HbmSpeciesTab.cpp
    modelers/HbmBiodiversityTab.cpp
    modelers/HbmLandscapeTab.cpp
    modelers/HbmPlanningTab.cpp
    modelers/GeOsirisWizard.cpp
    modelers/GeoProjectParamsTab.cpp
    modelers/GeoExternalFactorsTab.cpp
    modelers/GeoReddRulesTab.cpp
    modelers/GeoModelParamsTab.cpp
    modelers/GeoInputImagesTab.cpp
    modelers/GeoOpportunityCostTab.cpp
    modelers/GeoOutputParamsTab.cpp
    modelers/EtmWizard.cpp
    modelers/EtmExploreTab.cpp
    modelers/EtmAnalysisTab.cpp
    modelers/EtmPreprocessTab.cpp
    ${CMAKE_SOURCE_DIR}/include/MainWindow.h
    ${CMAKE_SOURCE_DIR}/include/MapCanvas.h
    ${CMAKE_SOURCE_DIR}/include/ModuleDialog.h
    ${CMAKE_SOURCE_DIR}/include/RasterInfoPanel.h
    ${CMAKE_SOURCE_DIR}/include/ModelerSession.h
    ${CMAKE_SOURCE_DIR}/include/ModelerTab.h
    ${CMAKE_SOURCE_DIR}/include/ModelerWizard.h
    ${CMAKE_SOURCE_DIR}/include/modelers/LcmWizard.h
    ${CMAKE_SOURCE_DIR}/include/modelers/HbmWizard.h
    ${CMAKE_SOURCE_DIR}/include/modelers/GeOsirisWizard.h
    ${CMAKE_SOURCE_DIR}/include/ChartWidget.h
    ${CMAKE_SOURCE_DIR}/include/modelers/EtmWizard.h
)

target_link_libraries(MyTerrsetImprovement PRIVATE
    aplaceholder_core
    aplaceholder_io
    aplaceholder_gis
    aplaceholder_imgproc
    aplaceholder_lcm
    aplaceholder_etm
    aplaceholder_hbm
    aplaceholder_esm
    aplaceholder_ccam
    Qt6::Widgets
    Qt6::Gui
    Qt6::OpenGLWidgets
    GDAL::GDAL
)

# Force-link all module static libs so REGISTER_MODULE global constructors run
set(_module_libs aplaceholder_gis aplaceholder_imgproc aplaceholder_lcm
    aplaceholder_etm aplaceholder_hbm aplaceholder_esm aplaceholder_ccam)
foreach(_lib ${_module_libs})
    target_link_options(MyTerrsetImprovement PRIVATE "/WHOLEARCHIVE:$<TARGET_FILE:${_lib}>")
endforeach()
