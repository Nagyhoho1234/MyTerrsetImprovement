#pragma once

#include <QMainWindow>
#include <QTreeWidget>
#include <QDockWidget>
#include <QStatusBar>
#include <QMenuBar>
#include <QToolBar>

namespace aplaceholder {

class MapCanvas;
class RasterInfoPanel;

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;

private slots:
    void openRaster();
    void saveRaster();
    void runModule(const QString& moduleName);
    void about();

private:
    void createMenus();
    void createToolBars();
    void createDockWidgets();
    void populateModuleTree();

    MapCanvas* m_canvas = nullptr;
    QTreeWidget* m_moduleTree = nullptr;
    QTreeWidget* m_layerTree = nullptr;
    RasterInfoPanel* m_infoPanel = nullptr;
};

} // namespace aplaceholder
