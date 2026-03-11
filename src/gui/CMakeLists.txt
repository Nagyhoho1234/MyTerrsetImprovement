add_executable(aplaceholder
    main.cpp
    MainWindow.cpp
    MapCanvas.cpp
    ModuleDialog.cpp
    RasterInfoPanel.cpp
    ${CMAKE_SOURCE_DIR}/include/MainWindow.h
    ${CMAKE_SOURCE_DIR}/include/MapCanvas.h
    ${CMAKE_SOURCE_DIR}/include/ModuleDialog.h
    ${CMAKE_SOURCE_DIR}/include/RasterInfoPanel.h
)

target_link_libraries(aplaceholder PRIVATE
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
    target_link_options(aplaceholder PRIVATE "/WHOLEARCHIVE:$<TARGET_FILE:${_lib}>")
endforeach()
