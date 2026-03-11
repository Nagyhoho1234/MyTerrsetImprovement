#include "MainWindow.h"
#include "MapCanvas.h"
#include "RasterInfoPanel.h"
#include "Module.h"
#include "ModuleRegistry.h"
#include "ModuleDialog.h"
#include "ChartWidget.h"
#include "GdalIO.h"
#include "modelers/LcmWizard.h"
#include "modelers/HbmWizard.h"
#include "modelers/GeOsirisWizard.h"
#include "modelers/EtmWizard.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QSplitter>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QApplication>
#include <QPushButton>
#include <QHeaderView>
#include <QCheckBox>
#include <QDir>
#include <QFileInfo>
#include <QGroupBox>
#include <QStandardPaths>
#include <QStyle>
#include <QPainter>
#include <QPixmap>

namespace aplaceholder {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    setWindowTitle(QString("MyTerrsetImprovement v%1  Geospatial Monitoring and Modeling System")
                       .arg(QApplication::applicationVersion()));
    resize(1400, 900);

    // Central widget — map canvas
    m_canvas = new MapCanvas(this);
    setCentralWidget(m_canvas);

    createMenus();
    createToolBars();
    createExplorerPanel();
    populateModuleCombo();

    // Default working directory
    m_workingDir = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    populateFileList(m_workingDir);

    statusBar()->showMessage("Ready");
}

MainWindow::~MainWindow() = default;

// ---------------------------------------------------------------------------
// Menus — matching TerrSet's menu bar layout
// ---------------------------------------------------------------------------

void MainWindow::addMod(QMenu* menu, const QString& moduleId, const QString& label)
{
    QString displayName = label.isEmpty() ? moduleId : label;
    bool registered = ModuleRegistry::instance().moduleNames().contains(moduleId);
    auto* action = menu->addAction(displayName, this, [this, moduleId]() { runModule(moduleId); });
    if (!registered)
        action->setEnabled(false);
}

void MainWindow::createMenus()
{
    // Apply dark maroon style to the menu bar
    menuBar()->setStyleSheet(
        "QMenuBar { background-color: #5B1A18; color: white; }"
        "QMenuBar::item { padding: 4px 10px; }"
        "QMenuBar::item:selected { background-color: #7B3A38; }"
        "QMenuBar::item:pressed { background-color: #9B5A58; }"
    );

    // --- File ---
    auto* fileMenu = menuBar()->addMenu("&File");
    buildFileMenu(fileMenu);

    // --- IDRISI GIS Analysis ---
    auto* gisMenu = menuBar()->addMenu("IDRISI &GIS Analysis");
    buildGisAnalysisMenu(gisMenu);

    // --- IDRISI Image Processing ---
    auto* imgMenu = menuBar()->addMenu("IDRISI Image &Processing");
    buildImageProcessingMenu(imgMenu);

    // --- Modeler direct actions (no submenus) ---
    menuBar()->addAction("&Land Change Modeler", this, [this]() {
        LcmWizard wizard(this);
        wizard.exec();
    });
    menuBar()->addAction("&Habitat and Biodiversity Modeler", this, [this]() {
        HbmWizard wizard(this);
        wizard.exec();
    });
    menuBar()->addAction("Ge&OSIRIS", this, [this]() {
        GeOsirisWizard wizard(this);
        wizard.exec();
    });
    menuBar()->addAction("&Earth Trends Modeler", this, [this]() {
        EtmWizard wizard(this);
        wizard.exec();
    });

    // Spacer — push Window and Help to the right
    menuBar()->setCornerWidget(new QWidget(), Qt::TopRightCorner);

    // --- Window ---
    auto* windowMenu = menuBar()->addMenu("&Window");
    windowMenu->addAction("Tile", this, []() {});
    windowMenu->addAction("Cascade", this, []() {});

    // --- Help ---
    auto* helpMenu = menuBar()->addMenu("Hel&p");
    helpMenu->addAction("Contents", this, []() {});
    helpMenu->addAction("Using Help", this, []() {});
    helpMenu->addSeparator();
    helpMenu->addAction("Quick Start", this, []() {});
    helpMenu->addSeparator();
    helpMenu->addAction("What's New", this, []() {});
    helpMenu->addAction("Manual", this, []() {});
    helpMenu->addAction("Tutorial", this, []() {});
    helpMenu->addSeparator();
    helpMenu->addAction("&About MyTerrsetImprovement", this, &MainWindow::about);
}

void MainWindow::buildFileMenu(QMenu* menu)
{
    // --- Display submenu ---
    auto* displayMenu = menu->addMenu("Display");
    displayMenu->addAction("DISPLAY Launcher", this, [this]() { runModule("DISPLAY"); });
    addMod(displayMenu, "PYRAMID");
    addMod(displayMenu, "ORTHO");
    displayMenu->addAction("VFIELD", this, [this]() { runModule("VFIELD"); });
    displayMenu->addAction("Fly Through", this, [this]() { runModule("FLYTHROUGH"); });
    displayMenu->addAction("Media Viewer", this, [this]() { runModule("MEDIAVIEWER"); });
    displayMenu->addAction("Create Photo Layer", this, [this]() { runModule("PHOTOLAYER"); });
    displayMenu->addSeparator();
    displayMenu->addAction("Symbol Workshop", this, [this]() { runModule("SYMBOLWORKSHOP"); });
    displayMenu->addSeparator();
    addMod(displayMenu, "COMPOSITE");
    addMod(displayMenu, "SEPARATE");
    addMod(displayMenu, "ILLUMINATE");
    displayMenu->addSeparator();
    addMod(displayMenu, "HISTO");
    addMod(displayMenu, "STRETCH");

    // --- Import submenu ---
    auto* importMenu = menu->addMenu("Import");
    importMenu->addAction("GDAL Raster Conversion Utility", this, [this]() { runModule("CONVERT"); });

    auto* impGeneral = importMenu->addMenu("General Conversion Tools");
    impGeneral->addAction("CRLF (Unix/Intel ASCII)", this, [this]() { runModule("CRLF"); });
    impGeneral->addAction("CSV2AVL (Comma Separated to AVL)", this, [this]() { runModule("CSV2AVL"); });
    impGeneral->addAction("GenericRaster (BIL, BIP, and BSQ)", this, [this]() { runModule("GENERICRASTER"); });
    impGeneral->addAction("NanFix (Replace NAN or INF values)", this, [this]() { runModule("NANFIX"); });
    impGeneral->addAction("SSTIDRIS (ASCII Grid/Spreadsheet)", this, [this]() { runModule("SSTIDRIS"); });
    impGeneral->addAction("VAR2FIX (Variable/Fixed Length ASCII)", this, [this]() { runModule("VAR2FIX"); });
    impGeneral->addAction("XYZIDRIS (ASCII XYZ)", this, [this]() { runModule("XYZIDRIS"); });

    auto* impGovt = importMenu->addMenu("Government / Data Provider Formats");
    impGovt->addAction("ASDIDRISI", this, [this]() { runModule("ASDIDRISI"); });
    impGovt->addAction("DigitalGlobe", this, [this]() { runModule("DIGITALGLOBE"); });
    impGovt->addAction("GACPIDRISI", this, [this]() { runModule("GACPIDRISI"); });
    impGovt->addAction("GeoEye", this, [this]() { runModule("GEOEYE"); });
    impGovt->addAction("GEOTIFF/TIFF", this, [this]() { runModule("GEOTIFF"); });
    impGovt->addAction("GPCIDRISI", this, [this]() { runModule("GPCIDRISI"); });
    impGovt->addAction("HDFEOS (HDF 4 or HDF-EOS 4)", this, [this]() { runModule("HDFEOS"); });
    impGovt->addAction("IKONOS", this, [this]() { runModule("IKONOS"); });
    impGovt->addAction("Landsat data archive", this, [this]() { runModule("LANDSAT"); });
    impGovt->addAction("MODISCONV", this, [this]() { runModule("MODISCONV"); });
    impGovt->addAction("MODISQC", this, [this]() { runModule("MODISQC"); });
    impGovt->addAction("NDVI3g", this, [this]() { runModule("NDVI3G"); });
    impGovt->addAction("NetCDF", this, [this]() { runModule("NETCDF"); });
    impGovt->addAction("PSDIDRISI", this, [this]() { runModule("PSDIDRISI"); });
    impGovt->addAction("QuickBird", this, [this]() { runModule("QUICKBIRD"); });
    impGovt->addAction("RADARSAT", this, [this]() { runModule("RADARSAT"); });
    impGovt->addAction("SACIDRIS", this, [this]() { runModule("SACIDRIS"); });
    impGovt->addAction("Sentinel", this, [this]() { runModule("SENTINEL"); });
    impGovt->addAction("SPOT (GeoTIFF/SPOT Scene/GeoSPOT-SPOTView)", this, [this]() { runModule("SPOT"); });
    impGovt->addAction("WorldView", this, [this]() { runModule("WORLDVIEW"); });
    impGovt->addAction("XYZMONTHLY", this, [this]() { runModule("XYZMONTHLY"); });
    impGovt->addSeparator();
    impGovt->addAction("DEMIDRIS (USGS)", this, [this]() { runModule("DEMIDRIS"); });
    impGovt->addAction("DLG (USGS)", this, [this]() { runModule("DLG"); });
    impGovt->addAction("SDTS", this, [this]() { runModule("SDTS"); });

    auto* impDesktop = importMenu->addMenu("Desktop Publishing Formats");
    impDesktop->addAction("BMPIDRIS", this, [this]() { runModule("BMPIDRIS"); });
    impDesktop->addAction("DXFIDRIS", this, [this]() { runModule("DXFIDRIS"); });
    impDesktop->addAction("GEOTIFF/TIFF", this, [this]() { runModule("GEOTIFF"); });
    impDesktop->addAction("JPGIDRIS", this, [this]() { runModule("JPGIDRIS"); });
    impDesktop->addAction("KMLIDRIS", this, [this]() { runModule("KMLIDRIS"); });

    auto* impSoftware = importMenu->addMenu("Software-Specific Formats");
    auto* impEsri = impSoftware->addMenu("ESRI Formats");
    impEsri->addAction("SHAPEIDR", this, [this]() { runModule("SHAPEIDR"); });
    impEsri->addAction("ARCRASTER", this, [this]() { runModule("ARCRASTER"); });
    impEsri->addAction("ARCIDRIS (GEN Format)", this, [this]() { runModule("ARCIDRIS"); });
    impSoftware->addAction("ECWIDRIS (ECW Format)", this, [this]() { runModule("ECWIDRIS"); });
    impSoftware->addAction("ENVIIDRIS (ENVI)", this, [this]() { runModule("ENVIIDRIS"); });
    impSoftware->addAction("ERDIDRIS (ERDAS)", this, [this]() { runModule("ERDIDRIS"); });
    impSoftware->addAction("ERMIDRIS (ERMapper)", this, [this]() { runModule("ERMIDRIS"); });
    impSoftware->addAction("GRASSIDR (GRASS)", this, [this]() { runModule("GRASSIDR"); });
    impSoftware->addAction("MIFIDRIS (MapInfo)", this, [this]() { runModule("MIFIDRIS"); });
    impSoftware->addAction("SHAPEIDR (Esri Shape)", this, [this]() { runModule("SHAPEIDR"); });
    impSoftware->addAction("SPLUSIDRIS (S-Plus)", this, [this]() { runModule("SPLUSIDRIS"); });
    impSoftware->addAction("SRFIDRIS (Surfer)", this, [this]() { runModule("SRFIDRIS"); });
    impSoftware->addAction("STATIDRIS (Statistica)", this, [this]() { runModule("STATIDRIS"); });
    impSoftware->addSeparator();
    impSoftware->addAction("IDRISI Vector Export (VXP)", this, [this]() { runModule("VXP"); });
    impSoftware->addAction("IDRISI File Conversion (16/32)", this, [this]() { runModule("IDRISICONV"); });

    // --- Export submenu ---
    auto* exportMenu = menu->addMenu("Export");
    exportMenu->addAction("GDAL Raster Conversion Utility", this, [this]() { runModule("CONVERT"); });

    auto* expGeneral = exportMenu->addMenu("General Conversions Tools");
    expGeneral->addAction("CRLF (Unix/Intel ASCII)", this, [this]() { runModule("CRLF"); });
    expGeneral->addAction("XYZIDRIS (ASCII XYZ)", this, [this]() { runModule("XYZIDRIS"); });
    expGeneral->addAction("AVL2CSV (Values to Comma Separated)", this, [this]() { runModule("AVL2CSV"); });

    auto* expDesktop = exportMenu->addMenu("Desktop Publishing Formats");
    expDesktop->addAction("BMPIDRIS", this, [this]() { runModule("BMPIDRIS"); });
    expDesktop->addAction("DXFIDRIS", this, [this]() { runModule("DXFIDRIS"); });
    expDesktop->addAction("GEOTIFF/TIFF", this, [this]() { runModule("GEOTIFF"); });
    expDesktop->addAction("JPGIDRIS", this, [this]() { runModule("JPGIDRIS"); });
    expDesktop->addAction("KMLIDRIS", this, [this]() { runModule("KMLIDRIS"); });

    auto* expSoftware = exportMenu->addMenu("Software-Specific Formats");
    auto* expEsri = expSoftware->addMenu("ESRI Formats");
    expEsri->addAction("SHAPEIDR", this, [this]() { runModule("SHAPEIDR"); });
    expEsri->addAction("ARCRASTER", this, [this]() { runModule("ARCRASTER"); });
    expEsri->addAction("ARCIDRIS (GEN Format)", this, [this]() { runModule("ARCIDRIS"); });
    expSoftware->addAction("ECWIDRIS (ECW Format)", this, [this]() { runModule("ECWIDRIS"); });
    expSoftware->addAction("ENVIIDRIS (ENVI)", this, [this]() { runModule("ENVIIDRIS"); });
    expSoftware->addAction("ERDIDRIS (ERDAS)", this, [this]() { runModule("ERDIDRIS"); });
    expSoftware->addAction("ERMIDRIS (ERMapper)", this, [this]() { runModule("ERMIDRIS"); });
    expSoftware->addAction("GEOTIFF/TIFF", this, [this]() { runModule("GEOTIFF"); });
    expSoftware->addAction("GRASSIDR (GRASS)", this, [this]() { runModule("GRASSIDR"); });
    expSoftware->addAction("MIFIDRIS (MapInfo)", this, [this]() { runModule("MIFIDRIS"); });
    expSoftware->addAction("SHAPEIDR (Esri Shape)", this, [this]() { runModule("SHAPEIDR"); });
    expSoftware->addAction("SPLUSIDRIS (S-Plus)", this, [this]() { runModule("SPLUSIDRIS"); });
    expSoftware->addAction("SRFIDRIS (Surfer)", this, [this]() { runModule("SRFIDRIS"); });
    expSoftware->addAction("STATIDRIS (Statistica)", this, [this]() { runModule("STATIDRIS"); });
    expSoftware->addSeparator();
    expSoftware->addAction("IDRISI Vector Export (VXP)", this, [this]() { runModule("VXP"); });
    expSoftware->addAction("IDRISI File Conversion (16/32)", this, [this]() { runModule("IDRISICONV"); });

    // --- Reformat submenu ---
    auto* reformatMenu = menu->addMenu("Reformat");
    addMod(reformatMenu, "CONVERT");
    reformatMenu->addSeparator();
    addMod(reformatMenu, "PROJECT");
    addMod(reformatMenu, "RESAMPLE");
    reformatMenu->addSeparator();
    addMod(reformatMenu, "WINDOW");
    addMod(reformatMenu, "EXPAND");
    addMod(reformatMenu, "CONTRACT");
    addMod(reformatMenu, "CONCAT");
    addMod(reformatMenu, "FLIP", "TRANSPOSE");
    reformatMenu->addSeparator();
    addMod(reformatMenu, "METAUPDATE");
    reformatMenu->addAction("RASTERVECTOR", this, [this]() { runModule("RASTERVECTOR"); });
    reformatMenu->addSeparator();
    addMod(reformatMenu, "GENERALIZE", "GENERALIZATION");
    reformatMenu->addAction("LINTOPNT", this, [this]() { runModule("LINTOPNT"); });
    reformatMenu->addAction("POLY2LINE", this, [this]() { runModule("POLY2LINE"); });

    // --- Data Entry submenu ---
    auto* dataEntryMenu = menu->addMenu("Data Entry");
    dataEntryMenu->addAction("Edit", this, [this]() { runModule("EDIT"); });
    addMod(dataEntryMenu, "ASSIGN");
    dataEntryMenu->addSeparator();
    dataEntryMenu->addAction("INITIAL", this, [this]() { runModule("INITIAL"); });
    dataEntryMenu->addAction("UPDATE", this, [this]() { runModule("UPDATE"); });
    dataEntryMenu->addAction("UTMRef", this, [this]() { runModule("UTMREF"); });
    dataEntryMenu->addSeparator();
    auto* interpMenu = dataEntryMenu->addMenu("Surface Interpolation");
    addMod(interpMenu, "INTERPOL", "IDW (INTERPOL)");
    addMod(interpMenu, "INTERCON");
    auto* tinMenu = interpMenu->addMenu("TIN Interpolation");
    addMod(tinMenu, "TIN");
    addMod(tinMenu, "TINSURF");
    addMod(tinMenu, "GENERALIZE", "GENERALIZATION");
    addMod(tinMenu, "LINTOPNT");
    addMod(tinMenu, "TINPREP");
    addMod(interpMenu, "THIESSEN");
    addMod(interpMenu, "TREND");
    dataEntryMenu->addSeparator();
    dataEntryMenu->addAction("Database Workshop", this, [this]() { runModule("DBWORKSHOP"); });

    menu->addSeparator();
    menu->addAction("Collection Editor", this, [this]() { runModule("COLLECTION_EDITOR"); });
    menu->addAction("Create TSF", this, [this]() { runModule("CREATE_TSF"); });

    menu->addSeparator();
    menu->addAction("User Preferences", this, []() {});

    menu->addSeparator();
    auto* helpMenu = menu->addMenu("Help");
    helpMenu->addAction("Contents", this, [this]() { about(); });
    helpMenu->addAction("Using Help", this, []() {});

    menu->addSeparator();
    menu->addAction("E&xit", this, &QWidget::close, QKeySequence::Quit);
}

void MainWindow::buildGisAnalysisMenu(QMenu* menu)
{
    // Database Query (screenshot 17)
    auto* dbQuery = menu->addMenu("Database Query");
    addMod(dbQuery, "RECLASS");
    addMod(dbQuery, "OVERLAY");
    dbQuery->addAction("LOGIC", this, [this]() { runModule("LOGIC"); });
    addMod(dbQuery, "CROSSTAB");
    dbQuery->addSeparator();
    dbQuery->addAction("Edit", this, [this]() { runModule("EDIT"); });
    addMod(dbQuery, "ASSIGN");
    addMod(dbQuery, "EXTRACT");
    dbQuery->addAction("BREAKOUT", this, [this]() { runModule("BREAKOUT"); });
    addMod(dbQuery, "CONSENSUS");
    dbQuery->addSeparator();
    addMod(dbQuery, "HISTO");
    addMod(dbQuery, "AREA");
    addMod(dbQuery, "PERIM");
    addMod(dbQuery, "PROFILE");
    addMod(dbQuery, "QUERY");
    addMod(dbQuery, "PCLASS");
    dbQuery->addSeparator();
    dbQuery->addAction("Database Workshop", this, [this]() { runModule("DBWORKSHOP"); });
    dbQuery->addAction("Image Calculator", this, [this]() { runModule("OVERLAY"); });

    // Mathematical Operators (screenshot 18)
    auto* mathOps = menu->addMenu("Mathematical Operators");
    addMod(mathOps, "OVERLAY");
    addMod(mathOps, "SCALAR");
    addMod(mathOps, "TRANSFORM");
    mathOps->addAction("LOGIC", this, [this]() { runModule("LOGIC"); });
    mathOps->addSeparator();
    mathOps->addAction("Image Calculator", this, [this]() { runModule("OVERLAY"); });

    // Distance Operators (screenshot 19)
    auto* distOps = menu->addMenu("Distance Operators");
    addMod(distOps, "DISTANCE");
    addMod(distOps, "SPDIST");
    addMod(distOps, "COST");
    addMod(distOps, "BUFFER");
    distOps->addSeparator();
    addMod(distOps, "VARCOST");
    addMod(distOps, "DISPERSE");
    addMod(distOps, "RESULTANT");
    addMod(distOps, "DECOMP");
    distOps->addSeparator();
    addMod(distOps, "PATHWAY");
    addMod(distOps, "ALLOCATE");
    addMod(distOps, "RELOCATE");
    addMod(distOps, "THIESSEN");

    // Context Operators (screenshot 20)
    auto* ctxOps = menu->addMenu("Context Operators");
    addMod(ctxOps, "SURFACE");
    addMod(ctxOps, "FILTER");
    ctxOps->addAction("PATTERN", this, [this]() { runModule("PATTERN"); });
    addMod(ctxOps, "TEXTURE");
    addMod(ctxOps, "GROUP");
    ctxOps->addAction("AREAFILTER", this, [this]() { runModule("AREAFILTER"); });
    addMod(ctxOps, "VIEWSHED");
    addMod(ctxOps, "WATERSHED");
    ctxOps->addAction("HINTERLAND", this, [this]() { runModule("HINTERLAND"); });
    ctxOps->addAction("PIXEL LOCATION", this, [this]() { runModule("PIXELLOCATION"); });

    // Statistics (screenshot 21)
    auto* stats = menu->addMenu("Statistics");
    addMod(stats, "HISTO");
    addMod(stats, "EXTRACT");
    stats->addAction("PATTERN", this, [this]() { runModule("PATTERN"); });
    stats->addAction("COUNT", this, [this]() { runModule("COUNT"); });
    stats->addSeparator();
    addMod(stats, "REGRESS");
    stats->addAction("MULTIREG", this, [this]() { runModule("MULTI_REG"); });
    addMod(stats, "LOGISTIC_REG", "LOGISTICREG");
    addMod(stats, "MULTILOGISTICREG");
    stats->addAction("KNNREGRESS", this, [this]() { runModule("KNNREGRESS"); });
    addMod(stats, "TREND");
    stats->addSeparator();
    addMod(stats, "AUTOCORR");
    addMod(stats, "DURBINWATSON", "DURBIN WATSON");
    addMod(stats, "QUADRAT");
    stats->addAction("CENTER", this, [this]() { runModule("CENTER"); });
    stats->addAction("CRATIO", this, [this]() { runModule("CRATIO"); });
    stats->addSeparator();
    addMod(stats, "CROSSTAB");
    addMod(stats, "JACCARD", "JACCARD (IoU)");
    addMod(stats, "VALIDATE");
    addMod(stats, "ROC");
    stats->addSeparator();
    addMod(stats, "SAMPLE");
    addMod(stats, "SUBSAMPLE");
    addMod(stats, "RANDOM");
    addMod(stats, "STANDARD");
    stats->addSeparator();
    stats->addAction("SPLUSIDRIS (S-Plus)", this, [this]() { runModule("SPLUSIDRIS"); });
    stats->addAction("STATIDRIS (Statistica)", this, [this]() { runModule("STATIDRIS"); });

    // Decision Support (screenshot 22)
    auto* decision = menu->addMenu("Decision Support");
    addMod(decision, "WEIGHT");
    addMod(decision, "MCE");
    addMod(decision, "RANK");
    addMod(decision, "TOPRANK");
    addMod(decision, "MOLA");
    decision->addSeparator();
    addMod(decision, "STANDARD");
    addMod(decision, "FUZZY");
    decision->addSeparator();
    decision->addAction("COUNT", this, [this]() { runModule("COUNT"); });
    decision->addAction("MDCHOICE", this, [this]() { runModule("MDCHOICE"); });
    decision->addSeparator();
    addMod(decision, "PCLASS");
    addMod(decision, "BAYES");
    decision->addAction("PCORBINARY", this, [this]() { runModule("PCORBINARY"); });
    addMod(decision, "BELIEF");
    addMod(decision, "RANDOM");
    decision->addSeparator();
    addMod(decision, "SAMPLE");
    addMod(decision, "ERRMAT");

    // Change / Time Series (screenshot 23)
    auto* changeSeries = menu->addMenu("Change / Time Series");
    addMod(changeSeries, "IMAGEDIFF");
    addMod(changeSeries, "IMAGERATIO");
    addMod(changeSeries, "CVA");
    addMod(changeSeries, "CALIBRATE");
    addMod(changeSeries, "CROSSTAB");
    changeSeries->addSeparator();
    addMod(changeSeries, "PROFILE");
    changeSeries->addAction("TFA", this, [this]() { runModule("TFA"); });
    addMod(changeSeries, "CORRELATE");
    addMod(changeSeries, "KENDALL");
    changeSeries->addAction("KENDALL TAU", this, [this]() { runModule("KENDALLTAU"); });
    addMod(changeSeries, "TIMESERIES_STATS", "TSTATS");
    addMod(changeSeries, "TCOR");
    changeSeries->addAction("Media Viewer", this, [this]() { runModule("MEDIAVIEWER"); });
    changeSeries->addSeparator();
    addMod(changeSeries, "MARKOV_CHAIN", "MARKOV");
    changeSeries->addAction("STCHOICE", this, [this]() { runModule("STCHOICE"); });
    changeSeries->addAction("DISAGGREGATE", this, [this]() { runModule("DISAGGREGATE"); });
    addMod(changeSeries, "NORMALIZE");
    addMod(changeSeries, "LOGISTIC_REG", "LOGISTICREG");
    changeSeries->addSeparator();
    changeSeries->addAction("CELLATOM", this, [this]() { runModule("CELLATOM"); });
    addMod(changeSeries, "CELLULAR_AUTOMATA", "CA_MARKOV");
    changeSeries->addAction("GEOMOD", this, [this]() { runModule("GEOMOD"); });
    changeSeries->addSeparator();
    addMod(changeSeries, "VALIDATE");
    addMod(changeSeries, "ROC");

    // Surface Analysis (screenshots 24-26)
    auto* surfAnalysis = menu->addMenu("Surface Analysis");

    auto* saInterp = surfAnalysis->addMenu("Interpolation");
    addMod(saInterp, "INTERPOL", "IDW (INTERPOL)");
    addMod(saInterp, "INTERCON");
    auto* saTin = saInterp->addMenu("TIN Interpolation");
    addMod(saTin, "TIN");
    addMod(saTin, "TINSURF");
    saTin->addSeparator();
    addMod(saTin, "GENERALIZE", "GENERALIZATION");
    saTin->addAction("LINTOPNT", this, [this]() { runModule("LINTOPNT"); });
    addMod(saTin, "TINPREP");
    addMod(saInterp, "THIESSEN");
    addMod(saInterp, "TREND");

    auto* saTopo = surfAnalysis->addMenu("Topographic Variables");
    addMod(saTopo, "SURFACE", "SLOPE");
    saTopo->addAction("ASPECT", this, [this]() { runModule("ASPECT"); });
    addMod(saTopo, "HILLSHADE");
    saTopo->addAction("CURVATURE", this, [this]() { runModule("CURVATURE"); });
    saTopo->addAction("FRACTAL", this, [this]() { runModule("FRACTAL"); });

    auto* saFeature = surfAnalysis->addMenu("Feature Extraction");
    addMod(saFeature, "CONTOUR");
    saFeature->addAction("TOPOSHAPE", this, [this]() { runModule("TOPOSHAPE"); });
    addMod(saFeature, "PITREMOVAL", "PIT REMOVAL");
    addMod(saFeature, "RUNOFF");
    addMod(saFeature, "FLOW");
    addMod(saFeature, "RUSLE");
    addMod(saFeature, "WATERSHED");
    addMod(saFeature, "SLOPELENGTH");
    addMod(saFeature, "FACET");
    saFeature->addAction("SEDIMENTATION", this, [this]() { runModule("SEDIMENTATION"); });

    // Model Deployment Tools (screenshot 27)
    auto* modelTools = menu->addMenu("Model Deployment Tools");
    modelTools->addAction("Image Calculator", this, [this]() { runModule("OVERLAY"); });
    modelTools->addAction("Macro Modeler", this, [this]() { runModule("MACROMODELER"); });
    modelTools->addAction("Spatial Decision Modeler", this, [this]() { runModule("SDM"); });
    modelTools->addSeparator();
    modelTools->addAction("Python", this, [this]() { runModule("PYTHON"); });
    modelTools->addAction("Run Macro", this, [this]() { runModule("RUNMACRO"); });
    modelTools->addAction("Edit", this, [this]() { runModule("EDIT"); });
}

void MainWindow::buildImageProcessingMenu(QMenu* menu)
{
    // Restoration (screenshot 28)
    auto* restoration = menu->addMenu("Restoration");
    addMod(restoration, "RESAMPLE");
    restoration->addAction("LANDSAT_UPSAMPLE", this, [this]() { runModule("LANDSAT_UPSAMPLE"); });
    restoration->addAction("LOCALAFFINE", this, [this]() { runModule("LOCALAFFINE"); });
    addMod(restoration, "MOSAIC");
    addMod(restoration, "DESTRIPE");
    addMod(restoration, "RADIANCE");
    addMod(restoration, "ATMOSC");
    addMod(restoration, "NDVICOMP");
    addMod(restoration, "SCREEN");

    // Enhancement (screenshot 29)
    auto* enhancement = menu->addMenu("Enhancement");
    addMod(enhancement, "STRETCH");
    addMod(enhancement, "COMPOSITE");
    addMod(enhancement, "FILTER");
    addMod(enhancement, "PANSHARPEN");
    enhancement->addAction("AREAFILTER", this, [this]() { runModule("AREAFILTER"); });

    // Transformation (screenshot 30)
    auto* transformation = menu->addMenu("Transformation");
    addMod(transformation, "PCA");
    transformation->addAction("CANCOMP", this, [this]() { runModule("CANCOMP"); });
    transformation->addAction("CANCOR", this, [this]() { runModule("CANCOR"); });
    addMod(transformation, "MNF");
    transformation->addAction("TFA", this, [this]() { runModule("TFA"); });
    addMod(transformation, "COLSPACE");
    addMod(transformation, "TEXTURE");
    transformation->addAction("THERMAL", this, [this]() { runModule("THERMAL"); });
    addMod(transformation, "NDVI", "VEGINDEX");
    transformation->addAction("SNOWINDEX", this, [this]() { runModule("SNOWINDEX"); });
    addMod(transformation, "WATERINDEX");
    transformation->addAction("SOILSALINITY", this, [this]() { runModule("SOILSALINITY"); });
    addMod(transformation, "TASSCAP");

    // Fourier Analysis (screenshot 31)
    auto* fourier = menu->addMenu("Fourier Analysis");
    addMod(fourier, "FOURIER");
    addMod(fourier, "ZEROPAD");
    addMod(fourier, "FILTFQ", "FILTERFQ");
    addMod(fourier, "FREQDIST");
    addMod(fourier, "DRAWFILT");

    // Signature Development (screenshot 32)
    auto* sigDev = menu->addMenu("Signature Development");
    addMod(sigDev, "MAKESIG");
    addMod(sigDev, "ENDSIG", "EndSig");
    addMod(sigDev, "FUZSIG");
    sigDev->addAction("PURIFY", this, [this]() { runModule("PURIFY"); });
    addMod(sigDev, "SUBSAMPLE");
    sigDev->addSeparator();
    addMod(sigDev, "HYPERSIG");
    addMod(sigDev, "HYPERAUTOSIG");
    sigDev->addSeparator();
    addMod(sigDev, "SIGCOMP");
    addMod(sigDev, "SEPSIG");
    addMod(sigDev, "SCATTER");

    // Hard Classifiers (screenshot 33)
    auto* hardClass = menu->addMenu("Hard Classifiers");
    addMod(hardClass, "PIPED");
    addMod(hardClass, "MINDIST");
    addMod(hardClass, "MAXLIKE");
    addMod(hardClass, "MULTILOGISTICREG");
    addMod(hardClass, "FISHER", "FISHER (LDA)");
    addMod(hardClass, "KNN");
    addMod(hardClass, "SEGCLASS");
    hardClass->addSeparator();
    addMod(hardClass, "CLUSTER");
    addMod(hardClass, "ISOCLUST");
    addMod(hardClass, "ISODATA");
    addMod(hardClass, "KMEANS");
    addMod(hardClass, "MAXSET");
    hardClass->addAction("CHAINCLUSTER", this, [this]() { runModule("CHAINCLUSTER"); });
    hardClass->addSeparator();
    addMod(hardClass, "MLP");
    addMod(hardClass, "SOM");
    addMod(hardClass, "FUZZYARTMAP", "Fuzzy ARTMAP");
    addMod(hardClass, "RBFNN");
    addMod(hardClass, "CTA");
    addMod(hardClass, "DECISIONFOREST");
    addMod(hardClass, "SVM");
    hardClass->addSeparator();
    addMod(hardClass, "CONSENSUS");

    // Soft Classifiers / Mixture Analysis (screenshot 34)
    auto* softClass = menu->addMenu("Soft Classifiers / Mixture Analysis");
    addMod(softClass, "BAYCLASS");
    addMod(softClass, "MAHALCLASS");
    addMod(softClass, "BELCLASS");
    addMod(softClass, "FUZCLASS");
    addMod(softClass, "MULTILOGISTICREG");
    addMod(softClass, "KNN");
    addMod(softClass, "MLP");
    addMod(softClass, "SOM");
    addMod(softClass, "DECISIONFOREST");
    addMod(softClass, "SVM");
    softClass->addSeparator();
    addMod(softClass, "UNMIX");
    addMod(softClass, "HYPERUSP");
    addMod(softClass, "HYPEROSP");
    addMod(softClass, "HYPERUNMIX");
    addMod(softClass, "HYPERABSORB");
    softClass->addSeparator();
    addMod(softClass, "BELCALC");
    addMod(softClass, "BELIEF");
    softClass->addSeparator();
    addMod(softClass, "HARDEN");
    softClass->addAction("PRIORCORRECT", this, [this]() { runModule("PRIORCORRECT"); });

    // Segmentation Classifiers (screenshot 35)
    auto* segClass = menu->addMenu("Segmentation Classifiers");
    addMod(segClass, "SEGMENT", "SEGMENTATION");
    addMod(segClass, "SEGTRAIN");
    addMod(segClass, "SEGCLASS");

    // Hyperspectral Image Analysis (screenshot 36)
    auto* hyperSpec = menu->addMenu("Hyperspectral Image Analysis");
    addMod(hyperSpec, "HYPERSIG");
    hyperSpec->addAction("ASDIDRISI", this, [this]() { runModule("ASDIDRISI"); });
    addMod(hyperSpec, "HYPERAUTOSIG");
    addMod(hyperSpec, "SCREEN");
    hyperSpec->addSeparator();
    addMod(hyperSpec, "HYPERSAM");
    addMod(hyperSpec, "HYPERMIN");
    hyperSpec->addSeparator();
    addMod(hyperSpec, "HYPERUSP");
    addMod(hyperSpec, "HYPEROSP");
    addMod(hyperSpec, "HYPERUNMIX");
    addMod(hyperSpec, "HYPERABSORB");

    // Accuracy Assessment
    auto* accuracy = menu->addMenu("Accuracy Assessment");
    addMod(accuracy, "SAMPLE");
    addMod(accuracy, "ERRMAT");
}

// ---------------------------------------------------------------------------
// Toolbars — matching TerrSet's icon bar + module combo
// ---------------------------------------------------------------------------

// Helper to create a simple colored icon
static QIcon makeIcon(const QColor& bg, const QString& letter, const QColor& fg = Qt::white)
{
    QPixmap px(20, 20);
    px.fill(Qt::transparent);
    QPainter p(&px);
    p.setRenderHint(QPainter::Antialiasing);
    p.setBrush(bg);
    p.setPen(Qt::NoPen);
    p.drawRoundedRect(1, 1, 18, 18, 3, 3);
    p.setPen(fg);
    p.setFont(QFont("Arial", 9, QFont::Bold));
    p.drawText(QRect(0, 0, 20, 20), Qt::AlignCenter, letter);
    p.end();
    return QIcon(px);
}

void MainWindow::createToolBars()
{
    auto* s = style();
    const int iconSz = 20;

    auto* mainBar = addToolBar("Main");
    mainBar->setMovable(false);
    mainBar->setIconSize(QSize(iconSz, iconSz));

    // File operations
    auto* actOpen = mainBar->addAction(s->standardIcon(QStyle::SP_DialogOpenButton), "Open");
    actOpen->setToolTip("Open Raster");
    connect(actOpen, &QAction::triggered, this, &MainWindow::openRaster);

    auto* actSave = mainBar->addAction(s->standardIcon(QStyle::SP_DialogSaveButton), "Save");
    actSave->setToolTip("Save Raster");
    connect(actSave, &QAction::triggered, this, &MainWindow::saveRaster);

    mainBar->addAction(s->standardIcon(QStyle::SP_FileIcon), "")->setToolTip("New");
    mainBar->addAction(s->standardIcon(QStyle::SP_DirOpenIcon), "")->setToolTip("Open Project");

    mainBar->addSeparator();

    // Edit / clipboard
    mainBar->addAction(makeIcon(QColor("#4A6FA5"), "C"), "")->setToolTip("Copy");
    mainBar->addAction(makeIcon(QColor("#4A6FA5"), "P"), "")->setToolTip("Paste");
    mainBar->addAction(s->standardIcon(QStyle::SP_BrowserReload), "")->setToolTip("Refresh");

    mainBar->addSeparator();

    // Display & navigation
    mainBar->addAction(makeIcon(QColor("#2E7D32"), "D"), "")->setToolTip("Display");
    mainBar->addAction(s->standardIcon(QStyle::SP_FileDialogInfoView), "")->setToolTip("Info");

    mainBar->addSeparator();

    // Map tools
    auto* actZoomIn = mainBar->addAction(makeIcon(QColor("#1565C0"), "+"), "");
    actZoomIn->setToolTip("Zoom In");
    connect(actZoomIn, &QAction::triggered, m_canvas, [this]() { m_canvas->zoomIn(); });

    auto* actZoomOut = mainBar->addAction(makeIcon(QColor("#1565C0"), "\u2212"), "");
    actZoomOut->setToolTip("Zoom Out");
    connect(actZoomOut, &QAction::triggered, m_canvas, [this]() { m_canvas->zoomOut(); });

    auto* actFit = mainBar->addAction(makeIcon(QColor("#1565C0"), "\u2B1C"), "");
    actFit->setToolTip("Zoom to Fit");
    connect(actFit, &QAction::triggered, m_canvas, [this]() { m_canvas->zoomToFit(); });

    mainBar->addAction(makeIcon(QColor("#1565C0"), "\u2316"), "")->setToolTip("Pan");

    mainBar->addSeparator();

    // GIS tool icons (matching TerrSet's toolbar icons)
    mainBar->addAction(makeIcon(QColor("#C62828"), "\u2261"), "")->setToolTip("OVERLAY");
    mainBar->addAction(makeIcon(QColor("#6A1B9A"), "\u2234"), "")->setToolTip("DISTANCE");
    mainBar->addAction(makeIcon(QColor("#AD1457"), "\u2248"), "")->setToolTip("SURFACE");
    mainBar->addAction(makeIcon(QColor("#00695C"), "T"), "")->setToolTip("TRANSFORM");

    mainBar->addSeparator();

    // Classification / Image Processing icons
    mainBar->addAction(makeIcon(QColor("#E65100"), "Cl"), "")->setToolTip("CLUSTER");
    mainBar->addAction(makeIcon(QColor("#BF360C"), "Ml"), "")->setToolTip("MAXLIKE");
    mainBar->addAction(makeIcon(QColor("#33691E"), "PC"), "")->setToolTip("PCA");

    mainBar->addSeparator();

    // Composite / Stretch / Filter
    mainBar->addAction(makeIcon(QColor("#0D47A1"), "Co"), "")->setToolTip("COMPOSITE");
    mainBar->addAction(makeIcon(QColor("#4E342E"), "St"), "")->setToolTip("STRETCH");
    mainBar->addAction(makeIcon(QColor("#37474F"), "Fi"), "")->setToolTip("FILTER");

    mainBar->addSeparator();

    // --- Module launcher combo (right side, like TerrSet) ---
    // Window & Help labels (matching TerrSet's toolbar layout)
    mainBar->addAction("Window", this, []() {});
    mainBar->addAction("Help", this, []() {});

    auto* spacer = new QWidget;
    spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    mainBar->addWidget(spacer);

    m_moduleCombo = new QComboBox;
    m_moduleCombo->setEditable(true);
    m_moduleCombo->setInsertPolicy(QComboBox::NoInsert);
    m_moduleCombo->setMinimumWidth(180);
    m_moduleCombo->setToolTip("Type a module name or select from the dropdown");
    mainBar->addWidget(m_moduleCombo);

    // Green play button
    auto* runBtn = new QPushButton;
    QPixmap playPx(20, 20);
    playPx.fill(Qt::transparent);
    QPainter pp(&playPx);
    pp.setRenderHint(QPainter::Antialiasing);
    pp.setBrush(QColor("#2E7D32"));
    pp.setPen(Qt::NoPen);
    QPolygon tri;
    tri << QPoint(4, 2) << QPoint(18, 10) << QPoint(4, 18);
    pp.drawPolygon(tri);
    pp.end();
    runBtn->setIcon(QIcon(playPx));
    runBtn->setIconSize(QSize(20, 20));
    runBtn->setFixedSize(26, 26);
    runBtn->setToolTip("Run selected module");
    runBtn->setFlat(true);
    mainBar->addWidget(runBtn);

    connect(runBtn, &QPushButton::clicked, this, &MainWindow::runSelectedModule);
    connect(m_moduleCombo, &QComboBox::activated, this, &MainWindow::runSelectedModule);
}

void MainWindow::populateModuleCombo()
{
    auto& reg = ModuleRegistry::instance();
    QStringList names = reg.moduleNames();
    names.sort(Qt::CaseInsensitive);
    m_moduleCombo->clear();
    m_moduleCombo->addItems(names);
}

void MainWindow::runSelectedModule()
{
    QString name = m_moduleCombo->currentText().trimmed();
    if (!name.isEmpty())
        runModule(name);
}

// ---------------------------------------------------------------------------
// TerrSet Explorer — left panel with Projects/Files/Filters + Metadata
// ---------------------------------------------------------------------------

void MainWindow::createExplorerPanel()
{
    auto* explorerDock = new QDockWidget("TerrSet Explorer", this);
    explorerDock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);

    auto* container = new QWidget;
    auto* vLayout = new QVBoxLayout(container);
    vLayout->setContentsMargins(0, 0, 0, 0);
    vLayout->setSpacing(0);

    // --- Tab widget: Projects | Files | Filters ---
    m_explorerTabs = new QTabWidget;

    // Projects tab
    m_projectTree = new QTreeWidget;
    m_projectTree->setHeaderHidden(true);
    auto* projectRoot = new QTreeWidgetItem(m_projectTree, {"Projects"});
    auto* defaultProject = new QTreeWidgetItem(projectRoot, {"Default"});
    Q_UNUSED(defaultProject);
    projectRoot->setExpanded(true);
    m_explorerTabs->addTab(m_projectTree, "Projects");

    // Files tab
    auto* filesPage = new QWidget;
    auto* filesLayout = new QVBoxLayout(filesPage);
    filesLayout->setContentsMargins(0, 0, 0, 0);

    m_fileTree = new QTreeWidget;
    m_fileTree->setHeaderHidden(true);
    m_fileTree->setRootIsDecorated(false);
    filesLayout->addWidget(m_fileTree);

    connect(m_fileTree, &QTreeWidget::itemClicked,
            this, &MainWindow::onFileSelected);
    connect(m_fileTree, &QTreeWidget::itemDoubleClicked,
            this, [this](QTreeWidgetItem* item, int) {
        QString path = m_workingDir + "/" + item->text(0);
        if (QFileInfo(path).isFile()) {
            // Try to open as raster
            auto raster = GdalIO::read(path);
            if (raster) {
                m_infoPanel->showRasterInfo(*raster, path);
                m_canvas->displayRaster(std::move(raster));
                statusBar()->showMessage("Loaded: " + path);
            }
        }
    });

    m_explorerTabs->addTab(filesPage, "Files");

    // Filters tab
    m_filtersPage = new QWidget;
    auto* filtersLayout = new QVBoxLayout(m_filtersPage);
    filtersLayout->setAlignment(Qt::AlignTop);

    auto* filterLabel = new QLabel("File Type Filters");
    filterLabel->setStyleSheet("font-weight: bold;");
    filtersLayout->addWidget(filterLabel);

    // File type checkboxes matching TerrSet
    const QStringList filterTypes = {
        "Raster Image (*.tif, *.tiff, *.rst)",
        "GeoTIFF (*.tif)",
        "NetCDF (*.nc)",
        "HDF (*.hdf)",
        "PNG Image (*.png)",
        "JPEG Image (*.jpg, *.jpeg)",
        "IDRISI Raster (*.rst)",
        "All Files (*.*)"
    };

    for (const auto& ft : filterTypes) {
        auto* cb = new QCheckBox(ft);
        cb->setChecked(ft.contains("*.tif") || ft.contains("*.rst") || ft.contains("All"));
        filtersLayout->addWidget(cb);
    }

    m_explorerTabs->addTab(m_filtersPage, "Filters");

    vLayout->addWidget(m_explorerTabs, 1);

    // --- Metadata panel (below explorer tabs) ---
    auto* metaLabel = new QLabel("  Metadata");
    metaLabel->setStyleSheet("background-color: #4682B4; color: white; "
                              "font-weight: bold; padding: 2px;");
    vLayout->addWidget(metaLabel);

    m_metadataTable = new QTableWidget;
    m_metadataTable->setColumnCount(2);
    m_metadataTable->setHorizontalHeaderLabels({"Property", "Value"});
    m_metadataTable->horizontalHeader()->setStretchLastSection(true);
    m_metadataTable->horizontalHeader()->setVisible(false);
    m_metadataTable->verticalHeader()->setVisible(false);
    m_metadataTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_metadataTable->setSelectionMode(QAbstractItemView::NoSelection);
    m_metadataTable->setShowGrid(true);
    m_metadataTable->setMaximumHeight(300);
    vLayout->addWidget(m_metadataTable);

    // Raster info panel (used internally for metadata display)
    m_infoPanel = new RasterInfoPanel;
    m_infoPanel->setVisible(false); // hidden, we use metadata table instead

    explorerDock->setWidget(container);
    addDockWidget(Qt::LeftDockWidgetArea, explorerDock);
}

void MainWindow::populateFileList(const QString& directory)
{
    m_fileTree->clear();

    QDir dir(directory);
    if (!dir.exists()) return;

    // Show path in the tree header-like first item
    auto* headerItem = new QTreeWidgetItem(m_fileTree, {dir.absolutePath() + " :"});
    QFont boldFont = headerItem->font(0);
    boldFont.setBold(true);
    headerItem->setFont(0, boldFont);

    QStringList filters = {"*.tif", "*.tiff", "*.rst", "*.img", "*.nc", "*.hdf",
                           "*.png", "*.jpg", "*.jpeg", "*.vct", "*.rgf", "*.tsf"};
    QFileInfoList files = dir.entryInfoList(filters, QDir::Files, QDir::Name);

    for (const auto& fi : files) {
        new QTreeWidgetItem(m_fileTree, {fi.fileName()});
    }
}

void MainWindow::onFileSelected(QTreeWidgetItem* item, int)
{
    QString filename = item->text(0);
    if (filename.endsWith(":")) return; // header item

    QString path = m_workingDir + "/" + filename;
    QFileInfo fi(path);

    if (!fi.exists()) return;

    // Show metadata in the table
    m_metadataTable->setRowCount(0);

    auto addRow = [this](const QString& prop, const QString& val) {
        int row = m_metadataTable->rowCount();
        m_metadataTable->insertRow(row);
        m_metadataTable->setItem(row, 0, new QTableWidgetItem(prop));
        m_metadataTable->setItem(row, 1, new QTableWidgetItem(val));
    };

    addRow("Name", fi.baseName());
    addRow("File format", fi.suffix().toUpper());
    addRow("File size", QString::number(fi.size()) + " bytes");

    // Try to read raster metadata via GDAL
    auto raster = GdalIO::read(path);
    if (raster) {
        addRow("Columns", QString::number(raster->cols()));
        addRow("Rows", QString::number(raster->rows()));
        addRow("Bands", QString::number(raster->bands()));
        addRow("Ref. system", raster->projection().left(60));

        auto gt = raster->geoTransform();
        addRow("Min. X", QString::number(gt.originX, 'f', 6));
        addRow("Max. X", QString::number(gt.originX + gt.pixelWidth * raster->cols(), 'f', 6));
        addRow("Min. Y", QString::number(gt.originY + gt.pixelHeight * raster->rows(), 'f', 6));
        addRow("Max. Y", QString::number(gt.originY, 'f', 6));
        addRow("Resolution", QString("%1 x %2").arg(std::abs(gt.pixelWidth)).arg(std::abs(gt.pixelHeight)));

        auto stats = raster->computeStats(0);
        addRow("Min. value", QString::number(stats.min, 'f', 4));
        addRow("Max. value", QString::number(stats.max, 'f', 4));
        addRow("Mean", QString::number(stats.mean, 'f', 4));
        addRow("Std. dev.", QString::number(stats.stddev, 'f', 4));
    }

    m_metadataTable->resizeColumnsToContents();
}

void MainWindow::refreshFileList()
{
    populateFileList(m_workingDir);
}

// ---------------------------------------------------------------------------
// Actions
// ---------------------------------------------------------------------------

void MainWindow::openRaster()
{
    QString path = QFileDialog::getOpenFileName(this, "Open Raster",
        m_workingDir,
        "All Raster Files (*.tif *.tiff *.rst *.img *.nc *.hdf *.png *.jpg);;All Files (*)");

    if (path.isEmpty())
        return;

    // Update working directory
    m_workingDir = QFileInfo(path).absolutePath();
    populateFileList(m_workingDir);

    auto raster = GdalIO::read(path);
    if (!raster) {
        QMessageBox::critical(this, "Error", "Failed to open: " + path);
        return;
    }

    m_canvas->displayRaster(std::move(raster));
    statusBar()->showMessage("Loaded: " + path);
}

void MainWindow::saveRaster()
{
    // TODO: implement save
    statusBar()->showMessage("Save not yet implemented");
}

void MainWindow::runModule(const QString& moduleName)
{
    auto module = ModuleRegistry::instance().create(moduleName);
    if (!module) {
        QMessageBox::information(this, moduleName,
            QString("The module '%1' is not yet available in this version.\n\n"
                    "It is planned for a future release. You can use the GDAL\n"
                    "Raster Conversion Utility or equivalent tools in the meantime.")
            .arg(moduleName));
        return;
    }

    ModuleDialog dialog(module.get(), this);
    if (dialog.exec() == QDialog::Accepted) {
        if (!module->execute()) {
            QMessageBox::critical(this, "Module Error", module->lastError());
        } else if (module->chartResult().hasData()) {
            // Show chart output
            ChartDialog chartDlg(module->chartResult(), this);
            chartDlg.exec();
        }
    }
}

void MainWindow::about()
{
    QMessageBox::about(this, "About MyTerrsetImprovement",
        "<h2>MyTerrsetImprovement v0.2.0</h2>"
        "<p>Open-source 64-bit raster GIS and remote sensing analysis suite</p>"
        "<p>Inspired by and built as an improvement upon TerrSet/IDRISI by Clark Labs (Clark University).</p>"
        "<hr>"
        "<p><b>Author:</b> Prof. Zsolt Zolt&aacute;n Feh&eacute;r Dr.<br>"
        "University of Debrecen, Hungary<br>"
        "ORCID: 0009-0007-6659-4197</p>"
        "<hr>"
        "<p>Built with Qt, GDAL, and PROJ.<br>"
        "Licensed under open-source terms.</p>");
}

} // namespace aplaceholder
