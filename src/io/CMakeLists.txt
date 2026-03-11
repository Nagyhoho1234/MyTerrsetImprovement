add_library(aplaceholder_io STATIC
    GdalIO.cpp
)

target_include_directories(aplaceholder_io PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(aplaceholder_io PUBLIC
    aplaceholder_core
    GDAL::GDAL
)
