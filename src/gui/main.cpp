#include <QApplication>
#include <QDir>
#include <QSettings>
#include <QStandardPaths>
#include <cstdlib>
#include "MainWindow.h"
#include "GdalIO.h"

int main(int argc, char* argv[])
{
    // Use INI files for settings — no Windows Registry, no admin needed
    QSettings::setDefaultFormat(QSettings::IniFormat);

    QApplication app(argc, argv);
    app.setApplicationName("MyTerrsetImprovement");
    app.setApplicationVersion("0.1.1-beta");
    app.setOrganizationName("University of Debrecen");

    // Auto-detect GDAL/PROJ data paths relative to the executable
    QString appDir = QCoreApplication::applicationDirPath();
    QString gdalData = appDir + "/share/gdal";
    QString projData = appDir + "/share/proj";
    if (QDir(gdalData).exists())
        qputenv("GDAL_DATA", gdalData.toUtf8());
    if (QDir(projData).exists())
        qputenv("PROJ_DATA", projData.toUtf8());

    // Set GDAL temp directory to user-writable location
    QString tempDir = QStandardPaths::writableLocation(QStandardPaths::TempLocation);
    qputenv("CPL_TMPDIR", tempDir.toUtf8());

    // Ensure user data directory exists (for project databases, recent files, etc.)
    QString userDataDir = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);
    QDir().mkpath(userDataDir);

    // Store settings INI next to the exe if writable (portable), otherwise in AppData
    QString portableIni = appDir + "/MyTerrsetImprovement.ini";
    QFile testWrite(portableIni);
    if (testWrite.open(QIODevice::WriteOnly | QIODevice::Append)) {
        testWrite.close();
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, appDir);
    }
    // else: falls back to %APPDATA% automatically

    // Initialize GDAL
    aplaceholder::GdalIO::initialize();

    aplaceholder::MainWindow window;
    window.show();

    int result = app.exec();

    aplaceholder::GdalIO::cleanup();
    return result;
}
