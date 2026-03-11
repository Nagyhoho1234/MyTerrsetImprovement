#include "MainWindow.h"
#include "MapCanvas.h"
#include "RasterInfoPanel.h"
#include "Module.h"
#include "ModuleRegistry.h"
#include "ModuleDialog.h"
#include "GdalIO.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QSplitter>
#include <QVBoxLayout>
#include <QLabel>
#include <QApplication>

namespace aplaceholder {

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    setWindowTitle(QString("APlaceholder v%1").arg(QApplication::applicationVersion()));
    resize(1400, 900);

    // Central widget — map canvas
    m_canvas = new MapCanvas(this);
    setCentralWidget(m_canvas);

    createMenus();
    createToolBars();
    createDockWidgets();
    populateModuleTree();

    statusBar()->showMessage("Ready");
}

MainWindow::~MainWindow() = default;

void MainWindow::createMenus()
{
    // File menu
    auto* fileMenu = menuBar()->addMenu("&File");
    fileMenu->addAction("&Open Raster...", this, &MainWindow::openRaster, QKeySequence::Open);
    fileMenu->addAction("&Save Raster...", this, &MainWindow::saveRaster, QKeySequence::Save);
    fileMenu->addSeparator();
    fileMenu->addAction("E&xit", this, &QWidget::close, QKeySequence::Quit);

    // GIS Analysis menu
    auto* gisMenu = menuBar()->addMenu("&GIS Analysis");
    auto addModuleAction = [&](QMenu* menu, const QString& name) {
        menu->addAction(name, this, [this, name]() { runModule(name); });
    };
    addModuleAction(gisMenu, "OVERLAY");
    addModuleAction(gisMenu, "DISTANCE");
    addModuleAction(gisMenu, "COST");
    addModuleAction(gisMenu, "SURFACE");
    addModuleAction(gisMenu, "RECLASS");
    addModuleAction(gisMenu, "CROSSTAB");
    addModuleAction(gisMenu, "FILTER");
    addModuleAction(gisMenu, "AUTOCORR");
    gisMenu->addSeparator();
    addModuleAction(gisMenu, "MCE");
    addModuleAction(gisMenu, "MOLA");

    // Image Processing menu
    auto* imgMenu = menuBar()->addMenu("&Image Processing");
    auto* classMenu = imgMenu->addMenu("Classification");
    addModuleAction(classMenu, "MAXLIKE");
    addModuleAction(classMenu, "MINDIST");
    addModuleAction(classMenu, "CLUSTER");
    auto* transMenu = imgMenu->addMenu("Transformation");
    addModuleAction(transMenu, "PCA");
    addModuleAction(transMenu, "NDVI");
    addModuleAction(imgMenu, "STRETCH");
    addModuleAction(imgMenu, "COMPOSITE");
    addModuleAction(imgMenu, "SEGMENT");

    // Modelers menu
    auto* modelMenu = menuBar()->addMenu("&Modelers");
    auto* lcmMenu = modelMenu->addMenu("Land Change Modeler");
    addModuleAction(lcmMenu, "Change Analysis");
    addModuleAction(lcmMenu, "Transition Potential");
    addModuleAction(lcmMenu, "Markov Chain");
    addModuleAction(lcmMenu, "Prediction");
    auto* etmMenu = modelMenu->addMenu("Earth Trends Modeler");
    addModuleAction(etmMenu, "Trend Analysis");
    addModuleAction(etmMenu, "Seasonal Decomposition");
    auto* hbmMenu = modelMenu->addMenu("Habitat && Biodiversity");
    addModuleAction(hbmMenu, "Species Distribution");
    addModuleAction(hbmMenu, "Biodiversity Metrics");

    // Help menu
    auto* helpMenu = menuBar()->addMenu("&Help");
    helpMenu->addAction("&About", this, &MainWindow::about);
}

void MainWindow::createToolBars()
{
    auto* fileBar = addToolBar("File");
    fileBar->addAction("Open", this, &MainWindow::openRaster);
    fileBar->addAction("Save", this, &MainWindow::saveRaster);

    auto* navBar = addToolBar("Navigation");
    navBar->addAction("Zoom In", m_canvas, [this]() { m_canvas->zoomIn(); });
    navBar->addAction("Zoom Out", m_canvas, [this]() { m_canvas->zoomOut(); });
    navBar->addAction("Fit", m_canvas, [this]() { m_canvas->zoomToFit(); });
}

void MainWindow::createDockWidgets()
{
    // Module Explorer (left dock)
    auto* moduleDock = new QDockWidget("Module Explorer", this);
    m_moduleTree = new QTreeWidget();
    m_moduleTree->setHeaderLabel("Modules");
    moduleDock->setWidget(m_moduleTree);
    addDockWidget(Qt::LeftDockWidgetArea, moduleDock);

    connect(m_moduleTree, &QTreeWidget::itemDoubleClicked,
            this, [this](QTreeWidgetItem* item) {
        if (item->childCount() == 0) {
            runModule(item->text(0));
        }
    });

    // Layer list (left dock, tabbed with modules)
    auto* layerDock = new QDockWidget("Layers", this);
    m_layerTree = new QTreeWidget();
    m_layerTree->setHeaderLabels({"Layer", "Type"});
    layerDock->setWidget(m_layerTree);
    addDockWidget(Qt::LeftDockWidgetArea, layerDock);
    tabifyDockWidget(moduleDock, layerDock);

    // Info panel (bottom dock)
    auto* infoDock = new QDockWidget("Raster Info", this);
    m_infoPanel = new RasterInfoPanel();
    infoDock->setWidget(m_infoPanel);
    addDockWidget(Qt::BottomDockWidgetArea, infoDock);
}

void MainWindow::populateModuleTree()
{
    auto& reg = ModuleRegistry::instance();
    auto categories = reg.categories();

    for (const auto& cat : categories) {
        auto* catItem = new QTreeWidgetItem(m_moduleTree, {cat});
        catItem->setExpanded(true);
        for (const auto& mod : reg.modulesInCategory(cat)) {
            new QTreeWidgetItem(catItem, {mod});
        }
    }

    // If no modules registered yet, show placeholder structure
    if (categories.isEmpty()) {
        auto* gis = new QTreeWidgetItem(m_moduleTree, {"GIS Analysis"});
        for (auto& m : {"OVERLAY","DISTANCE","COST","SURFACE","RECLASS","CROSSTAB",
                        "FILTER","AUTOCORR","MCE","MOLA","AREA"})
            new QTreeWidgetItem(gis, {m});

        auto* img = new QTreeWidgetItem(m_moduleTree, {"Image Processing"});
        for (auto& m : {"MAXLIKE","MINDIST","CLUSTER","PCA","STRETCH","NDVI",
                        "COMPOSITE","SEGMENT"})
            new QTreeWidgetItem(img, {m});

        auto* lcm = new QTreeWidgetItem(m_moduleTree, {"Land Change Modeler"});
        for (auto& m : {"Change Analysis","Transition Potential","Markov Chain",
                        "Cellular Automata","Prediction"})
            new QTreeWidgetItem(lcm, {m});

        auto* etm = new QTreeWidgetItem(m_moduleTree, {"Earth Trends Modeler"});
        for (auto& m : {"Trend Analysis","Seasonal Decomposition","Time Series Stats"})
            new QTreeWidgetItem(etm, {m});

        auto* hbm = new QTreeWidgetItem(m_moduleTree, {"Habitat & Biodiversity"});
        for (auto& m : {"Species Distribution","Biodiversity Metrics","Habitat Assessment"})
            new QTreeWidgetItem(hbm, {m});
    }
}

void MainWindow::openRaster()
{
    QString path = QFileDialog::getOpenFileName(this, "Open Raster",
        QString(), "All Raster Files (*.tif *.tiff *.rst *.img *.nc *.hdf *.png *.jpg);;All Files (*)");

    if (path.isEmpty())
        return;

    auto raster = GdalIO::read(path);
    if (!raster) {
        QMessageBox::critical(this, "Error", "Failed to open: " + path);
        return;
    }

    auto stats = raster->computeStats(0);
    m_infoPanel->showRasterInfo(*raster, path);
    m_canvas->displayRaster(std::move(raster));

    // Add to layer tree
    QFileInfo fi(path);
    auto* item = new QTreeWidgetItem(m_layerTree, {fi.fileName(), "Raster"});
    Q_UNUSED(item);

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
            QString("Module '%1' is not yet implemented.\n\n"
                    "This is a placeholder in the module explorer.")
            .arg(moduleName));
        return;
    }

    ModuleDialog dialog(module.get(), this);
    if (dialog.exec() == QDialog::Accepted) {
        if (!module->execute()) {
            QMessageBox::critical(this, "Module Error", module->lastError());
        }
    }
}

void MainWindow::about()
{
    QMessageBox::about(this, "About APlaceholder",
        "<h2>APlaceholder v2026.1</h2>"
        "<p>Open-source raster GIS and remote sensing analysis suite</p>"
        "<p>A modern, 64-bit alternative inspired by TerrSet/IDRISI.</p>"
        "<hr>"
        "<p><b>Author:</b> Prof. Zsolt Zolt&aacute;n Feh&eacute;r Dr.<br>"
        "University of Debrecen, Hungary<br>"
        "ORCID: 0009-0007-6659-4197</p>"
        "<hr>"
        "<p>Built with Qt, GDAL, and PROJ.<br>"
        "Licensed under open-source terms.</p>");
}

} // namespace aplaceholder
