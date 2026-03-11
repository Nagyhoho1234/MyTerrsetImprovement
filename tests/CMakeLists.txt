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
