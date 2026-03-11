cmake_minimum_required(VERSION 3.24)
project(MyTerrsetImprovement
    VERSION 0.2.0
    DESCRIPTION "Open-source raster GIS and remote sensing analysis suite"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Dependencies
find_package(Qt6 REQUIRED COMPONENTS Core Widgets Gui OpenGLWidgets)
find_package(GDAL REQUIRED)

# Core library
add_subdirectory(src/core)

# I/O library
add_subdirectory(src/io)

# Module libraries
add_subdirectory(src/modules)

# GUI
add_subdirectory(src/gui)

# Tests
option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# =============================================================================
# Packaging / Installer (CPack + NSIS)
# =============================================================================
include(packaging/deploy_windows.cmake)

# -- Install the main executable -----------------------------------------------
install(TARGETS MyTerrsetImprovement
    RUNTIME DESTINATION bin
    COMPONENT Runtime
)

# -- Qt6 deployment (windeployqt) ----------------------------------------------
# qt_generate_deploy_app_script requires Qt 6.3+.  It produces a cmake script
# that runs windeployqt at install time, copying Qt DLLs, plugins, and
# translations into the install prefix.
qt_generate_deploy_app_script(
    TARGET MyTerrsetImprovement
    OUTPUT_SCRIPT QT_DEPLOY_SCRIPT
    NO_UNSUPPORTED_PLATFORM_ERROR
)
install(SCRIPT ${QT_DEPLOY_SCRIPT} COMPONENT Runtime)

# -- GDAL data files -----------------------------------------------------------
deploy_gdal_runtime("bin" "share")

# -- Runtime DLL collection (GDAL, SQLite, proj, etc. from vcpkg) --------------
deploy_runtime_dependencies(MyTerrsetImprovement "bin")

# -- CPack metadata ------------------------------------------------------------
set(CPACK_PACKAGE_NAME "MyTerrsetImprovement")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_VENDOR "MyTerrsetImprovement Contributors")
set(CPACK_PACKAGE_CONTACT "https://github.com/MyTerrsetImprovement/MyTerrsetImprovement")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "MyTerrsetImprovement")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")

# -- NSIS-specific settings ----------------------------------------------------
set(CPACK_GENERATOR "NSIS")
set(CPACK_NSIS_DISPLAY_NAME "MyTerrsetImprovement ${PROJECT_VERSION}")
set(CPACK_NSIS_PACKAGE_NAME "MyTerrsetImprovement")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")

# Default install directory: C:\Program Files\MyTerrsetImprovement
set(CPACK_NSIS_INSTALL_DIRECTORY "MyTerrsetImprovement")

# Start Menu shortcut
set(CPACK_NSIS_CREATE_ICONS_EXTRA
    "CreateShortCut '$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\MyTerrsetImprovement.lnk' '$INSTDIR\\\\bin\\\\MyTerrsetImprovement.exe'"
)
set(CPACK_NSIS_DELETE_ICONS_EXTRA
    "Delete '$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\MyTerrsetImprovement.lnk'"
)

# Uninstaller
set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
set(CPACK_NSIS_MODIFY_PATH OFF)

# Optional: set an icon for the installer (uncomment and set path when available)
# set(CPACK_NSIS_MUI_ICON "${CMAKE_SOURCE_DIR}/resources/icons/app.ico")
# set(CPACK_NSIS_MUI_UNIICON "${CMAKE_SOURCE_DIR}/resources/icons/app.ico")

# Environment variable so the app can find GDAL/PROJ data at runtime
set(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "
    WriteRegStr HKLM 'SYSTEM\\\\CurrentControlSet\\\\Control\\\\Session Manager\\\\Environment' \
        'GDAL_DATA' '$INSTDIR\\\\share\\\\gdal'
    WriteRegStr HKLM 'SYSTEM\\\\CurrentControlSet\\\\Control\\\\Session Manager\\\\Environment' \
        'PROJ_DATA' '$INSTDIR\\\\share\\\\proj'
")
set(CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS "
    DeleteRegValue HKLM 'SYSTEM\\\\CurrentControlSet\\\\Control\\\\Session Manager\\\\Environment' 'GDAL_DATA'
    DeleteRegValue HKLM 'SYSTEM\\\\CurrentControlSet\\\\Control\\\\Session Manager\\\\Environment' 'PROJ_DATA'
")

include(CPack)
