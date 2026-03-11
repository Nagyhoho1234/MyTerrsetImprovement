#pragma once

#include <QMainWindow>
#include <QTreeWidget>
#include <QDockWidget>
#include <QStatusBar>
#include <QMenuBar>
#include <QToolBar>
#include <QComboBox>
#include <QTableWidget>
#include <QTabWidget>
#include <QLineEdit>

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
    void runSelectedModule();
    void about();
    void onFileSelected(QTreeWidgetItem* item, int column);
    void refreshFileList();

private:
    void createMenus();
    void createToolBars();
    void createExplorerPanel();
    void populateModuleCombo();
    void populateFileList(const QString& directory);

    // Menu builders
    void buildFileMenu(QMenu* menu);
    void buildGisAnalysisMenu(QMenu* menu);
    void buildImageProcessingMenu(QMenu* menu);

    // Helper to add a module action to a menu
    void addMod(QMenu* menu, const QString& moduleId, const QString& label = {});

    MapCanvas*       m_canvas       = nullptr;
    RasterInfoPanel* m_infoPanel    = nullptr;

    // TerrSet Explorer
    QTabWidget*   m_explorerTabs    = nullptr;
    QTreeWidget*  m_projectTree     = nullptr;
    QTreeWidget*  m_fileTree        = nullptr;
    QWidget*      m_filtersPage     = nullptr;
    QTableWidget* m_metadataTable   = nullptr;

    // Module launcher combo
    QComboBox*    m_moduleCombo     = nullptr;

    // Working directory
    QString       m_workingDir;
};

} // namespace aplaceholder
