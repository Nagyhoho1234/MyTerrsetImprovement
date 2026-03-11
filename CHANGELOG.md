# Changelog

## v0.2.0 (2026-03-11)

### New Features
- **4 Modeler Wizards**: LCM (6 tabs), HBM (4 tabs), GeOSIRIS (7 tabs), ETM (3 tabs)
  - All accessible as direct menu-bar actions (click to open, no submenu)
  - Full tabbed interfaces with session state, tab gating, progress bars
- **Chart Visualization System**: Custom QPainter-based chart rendering
  - 5 chart types: Histogram, Scatter, Line, Bar, Pie
  - Export to PNG/BMP/JPG image and CSV data
  - Auto-displays after chart-producing modules run
- **103 New Modules** (total now 308):
  - ~50 format converters (GeoTIFF, JPEG, BMP, KML, ENVI, ERDAS, Landsat, Sentinel, MODIS, etc.)
  - ~25 GIS utilities (ASPECT, CURVATURE, BREAKOUT, CELLATOM, GEOMOD, DISAGGREGATE, etc.)
  - ~15 image processing modules (CANCOMP, CANCOR, TFA, THERMAL, SNOWINDEX, FILTFQ, etc.)
  - ~14 UI/system modules (DISPLAY, EDIT, MACROMODELER, etc.)
- **16 Modules with Chart Output**: HISTO, SCATTER, PROFILE, CORRELATE, CROSSTAB, AREA, ROC, SENSITIVITY, ERRMAT, PCA (scree plot), NDVI, STRETCH, CLUSTER, FREQDIST, SIGCOMP, CALIBRATE

### Improvements
- All menus fully populated matching original TerrSet layout
  - File menu: Display, Import (50+ formats), Export, Reformat, Data Entry
  - GIS Analysis: 9 submenus with 120+ items
  - Image Processing: 10 submenus with 100+ items
- Unregistered modules shown grayed out instead of hidden
- Wizard test suite added (test_wizards.exe)

### Technical
- 169 files changed, +14,309 lines
- Build: zero errors (MSVC, C++20)
- Tests: smoke 7/7, comprehensive 50/50, wizard 12/12

---

## v0.1.1-beta (2026-03-11)
- Added Zenodo metadata (.zenodo.json) and CITATION.cff for DOI minting

## v0.1.0-beta (2026-03-11)
- Initial release: 205 modules, Qt6 GUI with TerrSet Explorer panel
- Core libraries: core, io, gis, imgproc, lcm, etm, hbm, esm, ccam
- GDAL raster I/O, IDRISI .rst format support
- Module auto-registration system
- Smoke test and comprehensive test suites
