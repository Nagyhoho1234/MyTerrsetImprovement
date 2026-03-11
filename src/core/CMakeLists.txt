add_library(aplaceholder_core STATIC
    Raster.cpp
    RasterMetadata.cpp
    Module.cpp
    ModuleRegistry.cpp
    Pipeline.cpp
    Parameter.cpp
    Database.cpp
)

target_include_directories(aplaceholder_core PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

# Find sqlite3 (installed via vcpkg as GDAL dependency)
find_package(unofficial-sqlite3 CONFIG REQUIRED)

target_link_libraries(aplaceholder_core PUBLIC
    GDAL::GDAL
    Qt6::Core
    unofficial::sqlite3::sqlite3
)
