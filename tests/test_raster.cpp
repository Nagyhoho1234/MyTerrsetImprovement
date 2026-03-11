#include <QtTest>
#include "Raster.h"

using namespace aplaceholder;

class TestRaster : public QObject {
    Q_OBJECT

private slots:
    void testCreate() {
        Raster r(100, 200, 1, DataType::Float64);
        QCOMPARE(r.cols(), 100);
        QCOMPARE(r.rows(), 200);
        QCOMPARE(r.bands(), 1);
        QCOMPARE(r.cellCount(), (int64_t)20000);
    }

    void testSetGetValue() {
        Raster r(10, 10);
        r.setValue(5, 3, 42.0);
        QCOMPARE(r.value(5, 3), 42.0);
        QCOMPARE(r.value(0, 0), 0.0);
    }

    void testBoundsCheck() {
        Raster r(10, 10);
        r.setNoDataValue(-9999);
        QCOMPARE(r.value(-1, 0), -9999.0);
        QCOMPARE(r.value(10, 0), -9999.0);
    }

    void testStats() {
        Raster r(3, 3);
        for (int i = 0; i < 9; ++i)
            r.data(0)[i] = i + 1.0;

        auto stats = r.computeStats(0);
        QCOMPARE(stats.min, 1.0);
        QCOMPARE(stats.max, 9.0);
        QCOMPARE(stats.mean, 5.0);
        QCOMPARE(stats.validCount, (int64_t)9);
    }

    void testGeoTransform() {
        Raster r(100, 100);
        GeoTransform gt;
        gt.originX = 500000;
        gt.originY = 5200000;
        gt.pixelWidth = 30;
        gt.pixelHeight = -30;
        r.setGeoTransform(gt);

        double x, y;
        r.colRowToXY(0, 0, x, y);
        QCOMPARE(x, 500000.0);
        QCOMPARE(y, 5200000.0);

        r.colRowToXY(10, 5, x, y);
        QCOMPARE(x, 500300.0);
        QCOMPARE(y, 5199850.0);
    }

    void testMultiBand() {
        Raster r(10, 10, 3);
        r.setValue(0, 0, 100, 0);
        r.setValue(0, 0, 200, 1);
        r.setValue(0, 0, 50, 2);
        QCOMPARE(r.value(0, 0, 0), 100.0);
        QCOMPARE(r.value(0, 0, 1), 200.0);
        QCOMPARE(r.value(0, 0, 2), 50.0);
    }

    void testNoData() {
        Raster r(5, 5);
        r.setNoDataValue(-9999);
        r.setValue(2, 2, -9999);
        r.setValue(0, 0, 10);
        r.setValue(1, 1, 20);

        auto stats = r.computeStats(0);
        // Only 2 non-zero, non-nodata values but there are 24 zeros too
        QVERIFY(stats.validCount == 24); // 25 - 1 nodata = 24
    }
};

QTEST_MAIN(TestRaster)
#include "test_raster.moc"
