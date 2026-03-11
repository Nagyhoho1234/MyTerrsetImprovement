#include <QApplication>
#include "MainWindow.h"
#include "GdalIO.h"

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);
    app.setApplicationName("APlaceholder");
    app.setApplicationVersion("2026.1");
    app.setOrganizationName("University of Debrecen");

    // Initialize GDAL
    aplaceholder::GdalIO::initialize();

    aplaceholder::MainWindow window;
    window.show();

    int result = app.exec();

    aplaceholder::GdalIO::cleanup();
    return result;
}
