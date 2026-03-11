#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <gdal_priv.h>
#include <gdalwarper.h>
#include <ogr_spatialref.h>
#include <cpl_string.h>

namespace aplaceholder {

class ProjectModule : public Module {
public:
    QString name() const override { return "PROJECT"; }
    QString description() const override {
        return "Reproject raster to a different coordinate reference system using GDAL warp.";
    }
    QString category() const override { return "GIS Analysis"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input raster"),
            ParameterDef::output("output", "Output raster"),
            ParameterDef::integer("target_epsg", "Target EPSG code", 4326, 1, 99999,
                "EPSG code for the target projection"),
            ParameterDef({
                "cell_size", "Output cell size", ParameterDef::Double, 0.0,
                {}, 0.0, 1e10, "Output cell size (0 for auto)", false
            }),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString outputPath = parameter("output").toString();
        int targetEpsg = parameter("target_epsg").toInt();
        double cellSize = parameter("cell_size").toDouble();

        reportProgress(0.1, "Opening source dataset...");

        GDALDataset* srcDs = static_cast<GDALDataset*>(
            GDALOpen(inputPath.toUtf8().constData(), GA_ReadOnly));
        if (!srcDs) {
            setError("Failed to open input raster with GDAL");
            return false;
        }

        // Set up target SRS
        OGRSpatialReference targetSrs;
        targetSrs.importFromEPSG(targetEpsg);
        char* targetWkt = nullptr;
        targetSrs.exportToWkt(&targetWkt);

        // Get source WKT
        const char* srcWkt = GDALGetProjectionRef(srcDs);

        reportProgress(0.2, "Computing output bounds...");

        // Create warp options
        GDALWarpOptions* warpOpts = GDALCreateWarpOptions();
        warpOpts->hSrcDS = srcDs;
        warpOpts->nBandCount = GDALGetRasterCount(srcDs);
        warpOpts->panSrcBands = static_cast<int*>(CPLMalloc(sizeof(int) * warpOpts->nBandCount));
        warpOpts->panDstBands = static_cast<int*>(CPLMalloc(sizeof(int) * warpOpts->nBandCount));
        for (int i = 0; i < warpOpts->nBandCount; ++i) {
            warpOpts->panSrcBands[i] = i + 1;
            warpOpts->panDstBands[i] = i + 1;
        }

        // Create coordinate transformation
        warpOpts->pTransformerArg = GDALCreateGenImgProjTransformer(
            srcDs, srcWkt, nullptr, targetWkt, FALSE, 0.0, 1);
        if (!warpOpts->pTransformerArg) {
            setError("Failed to create coordinate transformation");
            CPLFree(targetWkt);
            GDALClose(srcDs);
            GDALDestroyWarpOptions(warpOpts);
            return false;
        }
        warpOpts->pfnTransformer = GDALGenImgProjTransform;

        // Compute output extent and dimensions
        double geoTransform[6];
        int outCols = 0, outRows = 0;
        CPLErr err = GDALSuggestedWarpOutput(srcDs,
            warpOpts->pfnTransformer, warpOpts->pTransformerArg,
            geoTransform, &outCols, &outRows);

        if (err != CE_None) {
            setError("Failed to compute output extent");
            GDALDestroyGenImgProjTransformer(warpOpts->pTransformerArg);
            CPLFree(targetWkt);
            GDALClose(srcDs);
            GDALDestroyWarpOptions(warpOpts);
            return false;
        }

        // Apply custom cell size if specified
        if (cellSize > 0.0) {
            double width = geoTransform[1] * outCols;
            double height = std::abs(geoTransform[5]) * outRows;
            outCols = static_cast<int>(std::ceil(width / cellSize));
            outRows = static_cast<int>(std::ceil(height / cellSize));
            geoTransform[1] = cellSize;
            geoTransform[5] = -cellSize;
        }

        reportProgress(0.4, "Creating output dataset...");

        // Create output dataset
        GDALDriverH hDriver = GDALGetDriverByName("GTiff");
        GDALDataType eDT = GDALGetRasterDataType(GDALGetRasterBand(srcDs, 1));
        GDALDatasetH dstDs = GDALCreate(hDriver,
            outputPath.toUtf8().constData(),
            outCols, outRows, warpOpts->nBandCount,
            eDT, nullptr);

        if (!dstDs) {
            setError("Failed to create output dataset");
            GDALDestroyGenImgProjTransformer(warpOpts->pTransformerArg);
            CPLFree(targetWkt);
            GDALClose(srcDs);
            GDALDestroyWarpOptions(warpOpts);
            return false;
        }

        GDALSetGeoTransform(dstDs, geoTransform);
        GDALSetProjection(dstDs, targetWkt);

        // Set nodata on output bands
        for (int i = 1; i <= warpOpts->nBandCount; ++i) {
            int hasNoData = 0;
            double noData = GDALGetRasterNoDataValue(GDALGetRasterBand(srcDs, i), &hasNoData);
            if (hasNoData)
                GDALSetRasterNoDataValue(GDALGetRasterBand(dstDs, i), noData);
        }

        // Recreate transformer with the actual destination dataset
        GDALDestroyGenImgProjTransformer(warpOpts->pTransformerArg);
        warpOpts->hDstDS = dstDs;
        warpOpts->pTransformerArg = GDALCreateGenImgProjTransformer(
            srcDs, srcWkt, dstDs, targetWkt, FALSE, 0.0, 1);
        warpOpts->pfnTransformer = GDALGenImgProjTransform;

        reportProgress(0.5, "Warping...");

        // Execute warp
        GDALWarpOperation warpOp;
        warpOp.Initialize(warpOpts);
        err = warpOp.ChunkAndWarpImage(0, 0, outCols, outRows);

        // Cleanup
        GDALDestroyGenImgProjTransformer(warpOpts->pTransformerArg);
        warpOpts->pTransformerArg = nullptr;
        GDALDestroyWarpOptions(warpOpts);
        CPLFree(targetWkt);
        GDALClose(dstDs);
        GDALClose(srcDs);

        if (err != CE_None) {
            setError("Warp operation failed");
            return false;
        }

        reportProgress(1.0, "Done");
        return true;
    }
};

REGISTER_MODULE(ProjectModule)

} // namespace aplaceholder
