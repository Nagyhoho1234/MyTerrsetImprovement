@echo off
echo ============================================================
echo  APlaceholder - Dependency Setup Script
echo  Installs CMake, Qt 6, GDAL, and PROJ via vcpkg
echo ============================================================
echo.

:: Check for git
where git >nul 2>nul
if errorlevel 1 (
    echo ERROR: Git is required. Install from https://git-scm.com/
    exit /b 1
)

:: Install vcpkg if not present
if not exist "C:\vcpkg" (
    echo Installing vcpkg...
    cd /d C:\
    git clone https://github.com/microsoft/vcpkg.git
    cd vcpkg
    call bootstrap-vcpkg.bat
) else (
    echo vcpkg already installed at C:\vcpkg
)

:: Set vcpkg path
set VCPKG_ROOT=C:\vcpkg
set PATH=%VCPKG_ROOT%;%PATH%

:: Install dependencies (64-bit)
echo.
echo Installing GDAL (with PROJ, GeoTIFF, etc.)...
vcpkg install gdal:x64-windows

echo.
echo Installing Qt 6...
vcpkg install qt6:x64-windows

echo.
echo ============================================================
echo  Dependencies installed!
echo.
echo  To build:
echo    cd C:\Terrset64
echo    mkdir build ^& cd build
echo    cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
echo    cmake --build . --config Release
echo ============================================================
