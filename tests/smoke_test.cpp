#include "GdalIO.h"
#include "Module.h"
#include "ModuleRegistry.h"
#include "Raster.h"
#include <QCoreApplication>
#include <iostream>

using namespace aplaceholder;

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    GdalIO::initialize();

    int passed = 0, failed = 0;

    // Test 1: Read a GeoTIFF
    std::cout << "TEST 1: Read GeoTIFF... ";
    auto raster = GdalIO::read("C:/Terrset64/test_data/test_dem.tif");
    if (raster && raster->cols() == 100 && raster->rows() == 100) {
        std::cout << "PASS (" << raster->cols() << "x" << raster->rows()
                  << ", bands=" << raster->bands() << ")\n";
        passed++;
    } else {
        std::cout << "FAIL\n";
        failed++;
    }

    // Test 2: Check pixel values
    std::cout << "TEST 2: Pixel values... ";
    if (raster) {
        double center = raster->value(50, 50, 0);
        double corner = raster->value(0, 0, 0);
        std::cout << "PASS (center=" << center << ", corner=" << corner << ")\n";
        passed++;
    } else {
        std::cout << "FAIL (no raster)\n";
        failed++;
    }

    // Test 3: List registered modules
    std::cout << "TEST 3: Module registry... ";
    auto names = ModuleRegistry::instance().moduleNames();
    std::cout << "PASS (" << names.size() << " modules registered: ";
    for (int i = 0; i < std::min((int)names.size(), 5); ++i)
        std::cout << names[i].toStdString() << " ";
    if (names.size() > 5) std::cout << "...";
    std::cout << ")\n";
    passed++;

    // Test 4: Run SURFACE module (slope)
    std::cout << "TEST 4: Run SURFACE (slope)... ";
    auto surfMod = ModuleRegistry::instance().create("SURFACE");
    if (surfMod) {
        surfMod->setParameter("input", "C:/Terrset64/test_data/test_dem.tif");
        surfMod->setParameter("output", "C:/Terrset64/test_data/out_slope.tif");
        surfMod->setParameter("operation", "slope_degrees");
        if (surfMod->execute()) {
            auto slope = GdalIO::read("C:/Terrset64/test_data/out_slope.tif");
            if (slope && slope->cols() == 100) {
                std::cout << "PASS (slope center=" << slope->value(50, 50, 0) << ")\n";
                passed++;
            } else {
                std::cout << "FAIL (can't read output)\n";
                failed++;
            }
        } else {
            std::cout << "FAIL (" << surfMod->lastError().toStdString() << ")\n";
            failed++;
        }
    } else {
        std::cout << "FAIL (module not found)\n";
        failed++;
    }

    // Test 5: Run OVERLAY module (add)
    std::cout << "TEST 5: Run OVERLAY (add)... ";
    auto ovMod = ModuleRegistry::instance().create("OVERLAY");
    if (ovMod) {
        ovMod->setParameter("input1", "C:/Terrset64/test_data/test_dem.tif");
        ovMod->setParameter("input2", "C:/Terrset64/test_data/test_friction.tif");
        ovMod->setParameter("output", "C:/Terrset64/test_data/out_overlay.tif");
        ovMod->setParameter("operation", "add");
        if (ovMod->execute()) {
            auto result = GdalIO::read("C:/Terrset64/test_data/out_overlay.tif");
            if (result && result->cols() == 100) {
                std::cout << "PASS\n";
                passed++;
            } else {
                std::cout << "FAIL (can't read output)\n";
                failed++;
            }
        } else {
            std::cout << "FAIL (" << ovMod->lastError().toStdString() << ")\n";
            failed++;
        }
    } else {
        std::cout << "FAIL (module not found)\n";
        failed++;
    }

    // Test 6: Run DISTANCE module
    std::cout << "TEST 6: Run DISTANCE... ";
    auto distMod = ModuleRegistry::instance().create("DISTANCE");
    if (distMod) {
        distMod->setParameter("input", "C:/Terrset64/test_data/test_landcover.tif");
        distMod->setParameter("output", "C:/Terrset64/test_data/out_distance.tif");
        if (distMod->execute()) {
            std::cout << "PASS\n";
            passed++;
        } else {
            std::cout << "FAIL (" << distMod->lastError().toStdString() << ")\n";
            failed++;
        }
    } else {
        std::cout << "FAIL (module not found)\n";
        failed++;
    }

    // Test 7: Read IDRISI .rst format
    std::cout << "TEST 7: Read IDRISI .rst... ";
    auto rst = GdalIO::read("C:/Program Files (x86)/TerrSet_liberaGIS/Georef/Naduslat.rst");
    if (rst) {
        std::cout << "PASS (" << rst->cols() << "x" << rst->rows()
                  << ", bands=" << rst->bands() << ")\n";
        passed++;
    } else {
        std::cout << "FAIL (could not read IDRISI format)\n";
        failed++;
    }

    // Summary
    std::cout << "\n=== SMOKE TEST RESULTS ===\n";
    std::cout << "Passed: " << passed << "/" << (passed + failed) << "\n";
    if (failed > 0)
        std::cout << "FAILED: " << failed << " test(s)\n";
    else
        std::cout << "ALL TESTS PASSED\n";

    GdalIO::cleanup();
    return failed;
}
