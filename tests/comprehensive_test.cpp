#include "GdalIO.h"
#include "Module.h"
#include "ModuleRegistry.h"
#include "Raster.h"
#include <QCoreApplication>
#include <iostream>
#include <fstream>

using namespace aplaceholder;

static int passed = 0, failed = 0;

// Helper: run a module with params, check if it succeeds
bool runModule(const QString& name, std::initializer_list<std::pair<QString,QVariant>> params)
{
    auto mod = ModuleRegistry::instance().create(name);
    if (!mod) {
        std::cout << "FAIL (module '" << name.toStdString() << "' not found)\n";
        failed++;
        return false;
    }
    for (const auto& p : params)
        mod->setParameter(p.first, p.second);
    if (mod->execute()) {
        std::cout << "PASS\n";
        passed++;
        return true;
    } else {
        std::cout << "FAIL (" << mod->lastError().toStdString() << ")\n";
        failed++;
        return false;
    }
}

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    GdalIO::initialize();

    const QString D = "C:/Terrset64/test_data/";
    const QString dem = D + "test_dem.tif";
    const QString friction = D + "test_friction.tif";
    const QString lc = D + "test_landcover.tif";

    // Create a simple points CSV for interpolation tests
    {
        std::ofstream f((D + "test_points.csv").toStdString());
        f << "x,y,value\n";
        f << "500100,4499900,100\n500300,4499700,200\n500500,4499500,500\n";
        f << "500700,4499300,300\n500900,4499100,150\n500200,4499800,250\n";
        f << "500600,4499400,400\n500800,4499200,350\n500400,4499600,450\n";
        f << "500150,4499650,175\n";
    }

    // Create a reclass file for RECLASS test
    {
        std::ofstream f((D + "reclass.txt").toStdString());
        f << "1 1 10\n2 2 20\n3 3 30\n4 4 40\n5 5 50\n";
    }

    // Create a CSV value table for ecosystem services
    {
        std::ofstream f((D + "eco_values.csv").toStdString());
        f << "class_id,class_name,value_per_ha\n";
        f << "1,Forest,5000\n2,Grassland,2500\n3,Cropland,1500\n4,Urban,0\n5,Water,8000\n";
    }

    std::cout << "=== COMPREHENSIVE MODULE TEST ===\n";
    std::cout << "Testing modules across all categories...\n\n";

    // ===== GIS ANALYSIS =====
    std::cout << "--- GIS ANALYSIS ---\n";

    std::cout << " 1. SURFACE (slope)... ";
    runModule("SURFACE", {{"input", dem}, {"output", D+"out_slope.tif"}, {"operation", "slope_degrees"}});

    std::cout << " 2. HILLSHADE... ";
    runModule("HILLSHADE", {{"dem", dem}, {"output", D+"out_hillshade.tif"},
        {"sun_azimuth", 315.0}, {"sun_elevation", 45.0}});

    std::cout << " 3. OVERLAY (multiply)... ";
    runModule("OVERLAY", {{"input1", dem}, {"input2", friction}, {"output", D+"out_ov.tif"}, {"operation", "multiply"}});

    std::cout << " 4. SCALAR (sqrt)... ";
    runModule("SCALAR", {{"input", dem}, {"output", D+"out_scalar.tif"}, {"operation", 7}, {"value", 1.0}});

    std::cout << " 5. DISTANCE... ";
    runModule("DISTANCE", {{"input", lc}, {"output", D+"out_dist.tif"}});

    std::cout << " 6. COST... ";
    runModule("COST", {{"input_features", lc}, {"friction_surface", friction}, {"output", D+"out_cost.tif"}});

    std::cout << " 7. RECLASS... ";
    runModule("RECLASS", {{"input", lc}, {"output", D+"out_reclass.tif"},
        {"method", 1}, {"reclass_file", D+"reclass.txt"}});

    std::cout << " 8. FILTER (mean 3x3)... ";
    runModule("FILTER", {{"input", dem}, {"output", D+"out_filter.tif"}, {"filter_type", 0}, {"kernel_size", 3}});

    std::cout << " 9. NORMALIZE... ";
    runModule("NORMALIZE", {{"input", dem}, {"output", D+"out_norm.tif"}});

    std::cout << "10. STANDARD (z-score)... ";
    runModule("STANDARD", {{"input", dem}, {"output", D+"out_std.tif"}});

    std::cout << "11. BUFFER... ";
    runModule("BUFFER", {{"input", lc}, {"output", D+"out_buf.tif"}, {"buffer_distance", 50.0}});

    std::cout << "12. CONTOUR... ";
    runModule("CONTOUR", {{"input", dem}, {"output", D+"out_contour.tif"}, {"contour_interval", 50.0}});

    std::cout << "13. FLOW (D8)... ";
    runModule("FLOW", {{"dem", dem}, {"output_direction", D+"out_flowdir.tif"},
        {"output_accumulation", D+"out_flowacc.tif"}});

    std::cout << "14. PITREMOVAL... ";
    runModule("PITREMOVAL", {{"input_dem", dem}, {"output", D+"out_pitfill.tif"}});

    std::cout << "15. VIEWSHED... ";
    runModule("VIEWSHED", {{"dem", dem}, {"output", D+"out_view.tif"},
        {"observer_x", 500500.0}, {"observer_y", 4499500.0}});

    std::cout << "16. THIESSEN... ";
    runModule("THIESSEN", {{"points_file", D+"test_points.csv"}, {"reference_raster", dem},
        {"output", D+"out_thies.tif"}});

    std::cout << "17. INTERPOL (IDW)... ";
    runModule("INTERPOL", {{"points_file", D+"test_points.csv"}, {"reference_raster", dem},
        {"output", D+"out_idw.tif"}, {"power", 2.0}});

    std::cout << "18. TREND (1st order)... ";
    runModule("TREND", {{"points_file", D+"test_points.csv"}, {"reference_raster", dem},
        {"output", D+"out_trend.tif"}, {"order", 1}});

    std::cout << "19. CROSSTAB... ";
    runModule("CROSSTAB", {{"input1", lc}, {"input2", lc}, {"output", D+"out_cross.tif"},
        {"output_report", D+"out_cross.txt"}});

    std::cout << "20. HISTO... ";
    runModule("HISTO", {{"input", dem}, {"output_file", D+"out_histo.txt"}, {"num_bins", 50}});

    std::cout << "21. AREA... ";
    runModule("AREA", {{"input", lc}, {"output", D+"out_area.txt"}});

    std::cout << "22. LOGICAND... ";
    runModule("LOGICAND", {{"input1", lc}, {"input2", lc}, {"output", D+"out_and.tif"}});

    std::cout << "23. FLIP (vertical)... ";
    runModule("FLIP", {{"input", dem}, {"output", D+"out_flip.tif"}, {"direction", 1}});

    std::cout << "24. WINDOW (mean 5x5)... ";
    runModule("WINDOW", {{"input", dem}, {"output", D+"out_win.tif"}, {"statistic", 0}, {"window_size", 5}});

    std::cout << "25. RUSLE... ";
    runModule("RUSLE", {{"r_factor", friction}, {"k_factor", friction}, {"ls_factor", friction},
        {"c_factor", friction}, {"p_factor", friction}, {"output", D+"out_rusle.tif"}});

    std::cout << "26. CORRELATE... ";
    runModule("CORRELATE", {{"input1", dem}, {"input2", friction}, {"output_report", D+"out_corr.txt"}});

    std::cout << "27. VALIDATE... ";
    runModule("VALIDATE", {{"predicted", dem}, {"observed", friction}, {"output_report", D+"out_valid.txt"}});

    std::cout << "28. ERRMAT... ";
    runModule("ERRMAT", {{"classified", lc}, {"reference", lc}, {"output", D+"out_errmat.txt"}});

    // ===== IMAGE PROCESSING =====
    std::cout << "\n--- IMAGE PROCESSING ---\n";

    std::cout << "29. NDVI... ";
    runModule("NDVI", {{"nir_band", dem}, {"red_band", friction}, {"output", D+"out_ndvi.tif"}});

    std::cout << "30. STRETCH (linear)... ";
    runModule("STRETCH", {{"input", dem}, {"output", D+"out_stretch.tif"}, {"method", 0}});

    std::cout << "31. CLUSTER (K=3)... ";
    runModule("CLUSTER", {{"input", dem}, {"output", D+"out_cluster.tif"}, {"num_classes", 3}});

    std::cout << "32. PCA... ";
    runModule("PCA", {{"input", dem}, {"output", D+"out_pca.tif"}, {"num_components", 1}});

    std::cout << "33. SEGMENT... ";
    runModule("SEGMENT", {{"input", dem}, {"output", D+"out_seg.tif"}});

    std::cout << "34. TEXTURE (contrast)... ";
    runModule("TEXTURE", {{"input", dem}, {"output", D+"out_tex.tif"},
        {"measure", 0}, {"window_size", 5}, {"direction", 0}, {"num_levels", 16}});

    std::cout << "35. DENOISE (Lee)... ";
    runModule("DENOISE", {{"input", dem}, {"output", D+"out_denoise.tif"}, {"window_size", 5}});

    std::cout << "36. FOURIER (forward)... ";
    runModule("FOURIER", {{"input", dem}, {"output_magnitude", D+"out_fft_mag.tif"},
        {"output_phase", D+"out_fft_pha.tif"}, {"direction", 0}});

    std::cout << "37. ZEROPAD... ";
    runModule("ZEROPAD", {{"input", dem}, {"output", D+"out_zpad.tif"}});

    std::cout << "38. HARDEN... ";
    runModule("HARDEN", {{"input_rasters", dem+","+friction}, {"output", D+"out_hard.tif"}});

    std::cout << "39. MOSAIC... ";
    runModule("MOSAIC", {{"input_rasters", dem+","+dem}, {"output", D+"out_mosaic.tif"}, {"method", 0}});

    std::cout << "40. VEGINDEX (SAVI)... ";
    runModule("VEGINDEX", {{"nir_band", dem}, {"red_band", friction}, {"output", D+"out_savi.tif"},
        {"index_type", 0}});

    // ===== LAND CHANGE MODELER =====
    std::cout << "\n--- LAND CHANGE MODELER ---\n";

    std::cout << "41. CHANGE_ANALYSIS... ";
    runModule("CHANGE_ANALYSIS", {{"earlier_map", lc}, {"later_map", lc},
        {"output", D+"out_change.tif"}, {"output_report", D+"out_change.txt"}});

    std::cout << "42. MARKOV_CHAIN... ";
    runModule("MARKOV_CHAIN", {{"earlier_image", lc}, {"later_image", lc},
        {"output_matrix", D+"out_markov.txt"}});

    // ===== EARTH TRENDS MODELER =====
    std::cout << "\n--- EARTH TRENDS MODELER ---\n";

    std::cout << "43. TREND_ANALYSIS... ";
    runModule("TREND_ANALYSIS", {{"input_series", dem+","+friction+","+dem},
        {"output_trend", D+"out_tslope.tif"}, {"output_significance", D+"out_tpval.tif"}});

    std::cout << "44. TIMESERIES_STATS... ";
    runModule("TIMESERIES_STATS", {{"input_series", dem+","+friction+","+dem},
        {"output", D+"out_tsstats.tif"}, {"statistic", 0}});

    // ===== HBM =====
    std::cout << "\n--- HABITAT & BIODIVERSITY ---\n";

    std::cout << "45. BIODIVERSITY_METRICS (Shannon)... ";
    runModule("BIODIVERSITY_METRICS", {{"input", lc}, {"output", D+"out_shannon.tif"},
        {"metric", 0}, {"window_size", 5}});

    // ===== ESM =====
    std::cout << "\n--- ECOSYSTEM SERVICES ---\n";

    std::cout << "46. ECOSYSTEM_SERVICES... ";
    runModule("ECOSYSTEM_SERVICES", {{"land_cover", lc},
        {"valuation_table", D+"eco_values.csv"}, {"output", D+"out_ecoserv.tif"}});

    // ===== CCAM =====
    std::cout << "\n--- CLIMATE ADAPTATION ---\n";

    std::cout << "47. CLIMATE_ADAPTATION... ";
    runModule("CLIMATE_ADAPTATION", {{"exposure", dem}, {"sensitivity", friction},
        {"adaptive_capacity", dem}, {"output", D+"out_vuln.tif"}});

    // ===== RASTER UTILITIES =====
    std::cout << "\n--- RASTER UTILITIES ---\n";

    std::cout << "48. RESAMPLE... ";
    runModule("RESAMPLE", {{"input", dem}, {"output", D+"out_resamp.tif"},
        {"cell_size", 20.0}, {"method", 1}});

    std::cout << "49. GROUP... ";
    runModule("GROUP", {{"input", lc}, {"output", D+"out_group.tif"}, {"connectivity", 1}});

    std::cout << "50. CONVERT (to float32)... ";
    runModule("CONVERT", {{"input", lc}, {"output", D+"out_conv.tif"}, {"data_type", 3}});

    // ===== SUMMARY =====
    std::cout << "\n=== COMPREHENSIVE TEST RESULTS ===\n";
    std::cout << "Passed: " << passed << "/" << (passed + failed) << "\n";
    if (failed > 0)
        std::cout << "FAILED: " << failed << " test(s)\n";
    else
        std::cout << "ALL TESTS PASSED\n";

    GdalIO::cleanup();
    return failed;
}
