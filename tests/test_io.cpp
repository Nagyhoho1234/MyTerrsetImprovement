#include <QtTest>
#include "GdalIO.h"
#include "Raster.h"
#include <QTemporaryDir>

using namespace aplaceholder;

class TestIO : public QObject {
    Q_OBJECT

private slots:
    void initTestCase() {
        GdalIO::initialize();
    }

    void cleanupTestCase() {
        GdalIO::cleanup();
    }

    void testWriteReadGeoTiff() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        QString path = dir.path() + "/test.tif";

        // Create and write
        Raster r(50, 40, 1, DataType::Float64);
        GeoTransform gt;
        gt.originX = 500000;
        gt.originY = 5200000;
        gt.pixelWidth = 30;
        gt.pixelHeight = -30;
        r.setGeoTransform(gt);
        r.setNoDataValue(-9999);
        r.setValue(10, 20, 42.5);

        QVERIFY(GdalIO::write(r, path));

        // Read back
        auto r2 = GdalIO::read(path);
        QVERIFY(r2 != nullptr);
        QCOMPARE(r2->cols(), 50);
        QCOMPARE(r2->rows(), 40);
        QCOMPARE(r2->value(10, 20), 42.5);
    }

    void testReadMetadata() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        QString path = dir.path() + "/meta.tif";

        Raster r(200, 300, 3, DataType::Float32);
        QVERIFY(GdalIO::write(r, path));

        auto meta = GdalIO::readMetadata(path);
        QVERIFY(meta != nullptr);
        QCOMPARE(meta->cols(), 200);
        QCOMPARE(meta->rows(), 300);
        QCOMPARE(meta->bands(), 3);
    }

    void testSupportedFormats() {
        auto formats = GdalIO::supportedReadFormats();
        QVERIFY(formats.contains("GTiff"));
    }

    void testDriverDetection() {
        QCOMPARE(GdalIO::detectDriver("test.tif"), QString("GTiff"));
        QCOMPARE(GdalIO::detectDriver("test.rst"), QString("RST"));
        QCOMPARE(GdalIO::detectDriver("test.img"), QString("HFA"));
    }
};

QTEST_MAIN(TestIO)
#include "test_io.moc"
