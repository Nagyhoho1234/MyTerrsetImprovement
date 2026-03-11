# deploy_windows.cmake
# Handles finding and copying runtime dependencies for Windows packaging.
# Designed to work with vcpkg-installed DLLs and Qt's deployment tools.
#
# This script is included from the root CMakeLists.txt and runs at install time
# via install(CODE ...) and install(SCRIPT ...) directives.

# ---------------------------------------------------------------------------
# 1. Qt6 runtime deployment
# ---------------------------------------------------------------------------
# Qt provides qt_generate_deploy_app_script() (Qt 6.3+) which produces a script
# that calls windeployqt at install time. This is the officially supported way
# to bundle Qt DLLs, plugins (platforms/, imageformats/, etc.), and translations.
#
# Usage from the root CMakeLists.txt:
#   qt_generate_deploy_app_script(
#       TARGET aplaceholder
#       OUTPUT_SCRIPT QT_DEPLOY_SCRIPT
#       NO_UNSUPPORTED_PLATFORM_ERROR
#   )
#   install(SCRIPT ${QT_DEPLOY_SCRIPT})

# ---------------------------------------------------------------------------
# 2. GDAL runtime files
# ---------------------------------------------------------------------------
# When using vcpkg, GDAL DLLs end up in <vcpkg_installed>/<triplet>/bin.
# GDAL also needs its data directory (proj.db, coordinate system definitions,
# etc.) which lives in <vcpkg_installed>/<triplet>/share/gdal.

function(deploy_gdal_runtime INSTALL_BINDIR INSTALL_DATADIR)
    # --- DLLs ---
    # Find the GDAL shared library via the imported target.
    # GDAL::GDAL is defined by find_package(GDAL); its IMPORTED_LOCATION
    # (or IMPORTED_IMPLIB on Windows) points to the .lib, but the .dll
    # lives next to it in ../bin.  We also grab transitive deps (proj, etc.)
    # via file(GET_RUNTIME_DEPENDENCIES).
    install(CODE "
        set(_gdal_target \"\$<TARGET_FILE:GDAL::GDAL>\")
    " COMPONENT Runtime)

    # We rely on CMake 3.21+ file(GET_RUNTIME_DEPENDENCIES) which is invoked
    # in the main CPack install rules (see root CMakeLists.txt).  That call
    # will pick up gdal*.dll, proj*.dll, sqlite3.dll, etc. automatically.

    # --- GDAL data files (proj.db, gcs.csv, etc.) ---
    # Try GDAL_DATA env, then common vcpkg paths.
    if(DEFINED ENV{GDAL_DATA} AND EXISTS "$ENV{GDAL_DATA}")
        set(_gdal_data_dir "$ENV{GDAL_DATA}")
    elseif(DEFINED VCPKG_INSTALLED_DIR)
        # vcpkg standard layout
        set(_gdal_data_candidate "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/gdal")
        if(EXISTS "${_gdal_data_candidate}")
            set(_gdal_data_dir "${_gdal_data_candidate}")
        endif()
    endif()

    if(DEFINED _gdal_data_dir AND EXISTS "${_gdal_data_dir}")
        install(DIRECTORY "${_gdal_data_dir}/"
            DESTINATION "${INSTALL_DATADIR}/gdal"
            COMPONENT Runtime
        )
        message(STATUS "GDAL data directory for packaging: ${_gdal_data_dir}")
    else()
        message(WARNING
            "Could not locate GDAL data directory. "
            "Set GDAL_DATA environment variable or ensure vcpkg layout is standard. "
            "The installer will still include the DLLs but may lack projection data."
        )
    endif()

    # --- PROJ data (proj.db) ---
    # proj data is often separate from gdal data in vcpkg
    if(DEFINED VCPKG_INSTALLED_DIR)
        set(_proj_data_candidate "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share/proj")
        if(EXISTS "${_proj_data_candidate}")
            install(DIRECTORY "${_proj_data_candidate}/"
                DESTINATION "${INSTALL_DATADIR}/proj"
                COMPONENT Runtime
            )
            message(STATUS "PROJ data directory for packaging: ${_proj_data_candidate}")
        endif()
    endif()
endfunction()


# ---------------------------------------------------------------------------
# 3. Generic vcpkg DLL collection via file(GET_RUNTIME_DEPENDENCIES)
# ---------------------------------------------------------------------------
# This function generates install-time code that introspects the built
# executable and copies every DLL it depends on (excluding Windows system
# DLLs).  Works for GDAL, SQLite, proj, tiff, curl, etc. -- anything
# vcpkg put in the bin directory.

function(deploy_runtime_dependencies TARGET_NAME INSTALL_BINDIR)
    # Build a list of directories CMake should search for DLLs.
    # At configure time we collect hints; at install time the generator
    # expressions are resolved.
    set(_search_dirs "")

    # vcpkg bin directory
    if(DEFINED VCPKG_INSTALLED_DIR)
        list(APPEND _search_dirs "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin")
    endif()

    # Qt bin directory
    get_target_property(_qt_core_loc Qt6::Core IMPORTED_LOCATION_RELEASE)
    if(_qt_core_loc)
        get_filename_component(_qt_bin_dir "${_qt_core_loc}" DIRECTORY)
        list(APPEND _search_dirs "${_qt_bin_dir}")
    else()
        # Fallback: try the non-config-specific property
        get_target_property(_qt_core_loc Qt6::Core IMPORTED_LOCATION)
        if(_qt_core_loc)
            get_filename_component(_qt_bin_dir "${_qt_core_loc}" DIRECTORY)
            list(APPEND _search_dirs "${_qt_bin_dir}")
        endif()
    endif()

    install(CODE "
        # Resolve runtime dependencies of the main executable.
        file(GET_RUNTIME_DEPENDENCIES
            RESOLVED_DEPENDENCIES_VAR _resolved
            UNRESOLVED_DEPENDENCIES_VAR _unresolved
            CONFLICTING_DEPENDENCIES_PREFIX _conflicts
            EXECUTABLES \"\${CMAKE_INSTALL_PREFIX}/${INSTALL_BINDIR}/$<TARGET_FILE_NAME:${TARGET_NAME}>\"
            DIRECTORIES ${_search_dirs}
            PRE_EXCLUDE_REGEXES
                [[api-ms-win-.*]]
                [[ext-ms-.*]]
            POST_EXCLUDE_REGEXES
                [[.*[/\\\\]Windows[/\\\\].*]]
                [[.*[/\\\\]system32[/\\\\].*]]
                [[.*[/\\\\]SysWOW64[/\\\\].*]]
        )

        message(STATUS \"Resolved runtime dependencies: \${_resolved}\")
        if(_unresolved)
            message(WARNING \"Unresolved runtime dependencies: \${_unresolved}\")
        endif()

        foreach(_dep IN LISTS _resolved)
            file(INSTALL
                DESTINATION \"\${CMAKE_INSTALL_PREFIX}/${INSTALL_BINDIR}\"
                TYPE SHARED_LIBRARY
                FILES \"\${_dep}\"
            )
        endforeach()
    " COMPONENT Runtime)
endfunction()
