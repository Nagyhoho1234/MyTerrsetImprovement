# MyTerrsetImprovement

> **BETA** -- This software is in active development. Expect bugs and incomplete features.

An open-source, 64-bit raster GIS and remote sensing analysis suite inspired by and built as an improvement upon **[TerrSet / IDRISI](https://clarklabs.org/terrset/)** by Clark Labs (Clark University).

---

## Acknowledgment & Attribution

> **This project would not exist without the decades of work by Dr. J. Ronald Eastman and the Clark Labs team at Clark University.** TerrSet (formerly IDRISI) has been the gold standard for raster-based GIS and remote sensing education since the 1980s. The analytical methods, module designs, workflow patterns, and documentation that guided this implementation all originate from the TerrSet software suite.
>
> **TerrSet / IDRISI** is a product of **Clark Labs, Clark University**, Worcester, Massachusetts, USA.
> - Website: https://clarklabs.org/terrset/
> - As of 2024, TerrSet has been released as freeware under the name **liberaGIS**.
> - The TerrSet manual and tutorial documentation are licensed under **CC BY-NC 4.0**.
>
> **We strongly encourage users to explore the original TerrSet software and its extensive documentation.** This project is not affiliated with, endorsed by, or a replacement for TerrSet -- it is an independent reimplementation that aims to bring similar analytical capabilities to a modern, open-source, 64-bit platform.

---

## Why This Project?

TerrSet is an exceptional analytical tool, but it has inherent limitations due to its legacy architecture:

| Limitation in TerrSet | Improvement Here |
|----------------------|------------------|
| 32-bit only (max ~2 GB RAM) | **Native 64-bit**, no memory limits |
| Float32 maximum precision | **Float64 (double precision)** support |
| Proprietary .rst raster format | **160+ formats** via GDAL (GeoTIFF, HDF, NetCDF, .rst, ...) |
| Microsoft Access/Jet database | **SQLite** -- embedded, zero-configuration |
| Windows-only, closed source | **Open source** (MIT license), cross-platform potential |
| No multi-threading | **Multi-threaded** capable (C++20) |
| Manual dependency management | **Bundled installer** with all dependencies |

## Features

- **205 analytical modules** covering:
  - **GIS Analysis** (91 modules): terrain analysis, hydrology, cost distance, interpolation, Boolean operations, statistics, viewshed, buffering, and more
  - **Image Processing** (68 modules): supervised & unsupervised classification (MaxLike, SVM, MLP, Random Forest, KNN, ISODATA, K-Means), PCA, spectral unmixing, Fourier analysis, GLCM texture, atmospheric correction, pan-sharpening, hyperspectral analysis, and more
  - **Earth Trends Modeler** (17 modules): time series analysis, HANTS, wavelet decomposition, CCA, harmonic regression, seasonal trend analysis
  - **Land Change Modeler** (9 modules): change analysis, Markov chains, cellular automata, transition potential modeling, logistic regression
  - **Habitat & Biodiversity Modeler** (7 modules): biodiversity metrics, habitat suitability, species distribution (BIOCLIM)
  - **Ecosystem Services Modeler** (7 modules): carbon stock, pollination, coastal vulnerability, wave energy, sediment delivery
  - **Climate Change Adaptation** (3 modules): vulnerability assessment, crop suitability, scenario generation

- **Native IDRISI .rst format support** -- reads TerrSet's proprietary raster format directly via GDAL
- **Qt 6 GUI** with auto-generated parameter dialogs for every module
- **Pipeline system** for chaining modules (JSON-serializable, like TerrSet's Macro Modeler)
- **SQLite database** replacing TerrSet's Access/Jet Engine dependency

## Technology Stack

- **C++20** with MSVC 2022
- **Qt 6.8** for GUI
- **GDAL 3.12** for raster I/O (160+ formats)
- **SQLite 3** for attribute database
- **CMake** build system with **vcpkg** package management

## Building from Source

### Prerequisites

- Windows 10/11 (64-bit)
- Visual Studio 2022 Build Tools (or full VS 2022)
- CMake 3.24+
- vcpkg (for GDAL and dependencies)
- Qt 6.8+ (MSVC 2022 64-bit)

### Build Steps

```bash
# 1. Install dependencies
setup_deps.bat

# 2. Configure
cmake -B build -S . ^
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake ^
  -DCMAKE_PREFIX_PATH=C:/Qt/6.8.3/msvc2022_64

# 3. Build
cmake --build build --config Release

# 4. Deploy runtime DLLs
windeployqt build/bin/Release/aplaceholder.exe --release
copy C:\vcpkg\installed\x64-windows\bin\*.dll build\bin\Release\

# 5. Run
set GDAL_DATA=build/bin/Release/share/gdal
set PROJ_DATA=build/bin/Release/share/proj
build\bin\Release\aplaceholder.exe
```

## Module Reference

All 205 modules are auto-registered and available through the GUI or programmatically via the `ModuleRegistry`. The module structure closely follows TerrSet's organization to ease the transition for existing TerrSet users.

For detailed module documentation, refer to the **TerrSet Manual** (available from [Clark Labs](https://clarklabs.org/terrset/)) which describes the analytical methods implemented here.

## Project Status

**Current: BETA v0.1.0**

- [x] Core raster engine (Float64, multi-band, NoData)
- [x] GDAL I/O (160+ formats including IDRISI .rst)
- [x] 205 analytical modules implemented
- [x] Qt 6 GUI with auto-generated dialogs
- [x] SQLite database layer
- [x] Pipeline/macro system
- [x] 50-test comprehensive test suite passing
- [ ] Comprehensive per-module validation against TerrSet outputs
- [ ] Cross-platform build (Linux, macOS)
- [ ] Plugin system for user-defined modules
- [ ] Full documentation
- [ ] NSIS installer packaging

## License

This project is licensed under the **MIT License** -- see [LICENSE](LICENSE).

The analytical methods and module designs are inspired by TerrSet/IDRISI. The TerrSet manual and tutorial materials are copyright Clark Labs and licensed under CC BY-NC 4.0. This project contains no TerrSet source code.

## References

- Eastman, J.R. (2020). *TerrSet Manual*. Clark Labs, Clark University.
- Eastman, J.R. (2020). *TerrSet Tutorial*. Clark Labs, Clark University.
- GDAL/OGR contributors (2024). *GDAL/OGR Geospatial Data Abstraction software Library*. https://gdal.org
- Qt Project (2024). *Qt 6 Framework*. https://www.qt.io
