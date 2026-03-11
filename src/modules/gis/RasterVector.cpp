#include "Module.h"
#include "ModuleRegistry.h"
#include "GdalIO.h"
#include <cmath>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <gdal_alg.h>

namespace aplaceholder {

class RasterVectorModule : public Module {
public:
    QString name() const override { return "RASTERVECTOR"; }
    QString description() const override {
        return "Convert between raster and vector formats using GDAL.";
    }
    QString category() const override { return "File"; }

    std::vector<ParameterDef> parameterDefs() const override {
        return {
            ParameterDef::file("input", "Input file (raster or vector)"),
            ParameterDef::output("output", "Output file"),
            ParameterDef::combo("direction", "Conversion direction",
                {"Raster to Vector", "Vector to Raster"}, 0,
                "Direction of conversion"),
        };
    }

    bool execute() override {
        QString inputPath = parameter("input").toString();
        QString outputPath = parameter("output").toString();
        int direction = parameter("direction").toInt();

        if (direction == 0) {
            return rasterToVector(inputPath, outputPath);
        } else {
            return vectorToRaster(inputPath, outputPath);
        }
    }

private:
    bool rasterToVector(const QString& inputPath, const QString& outputPath) {
        GDALDataset* srcDs = static_cast<GDALDataset*>(
            GDALOpen(inputPath.toStdString().c_str(), GA_ReadOnly));
        if (!srcDs) {
            setError("Failed to open input raster: " + inputPath);
            return false;
        }

        GDALRasterBand* band = srcDs->GetRasterBand(1);
        if (!band) {
            GDALClose(srcDs);
            setError("Failed to read raster band from input");
            return false;
        }

        reportProgress(0.2, "Creating vector output...");

        GDALDriver* vecDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
        if (!vecDriver) {
            GDALClose(srcDs);
            setError("ESRI Shapefile driver not available");
            return false;
        }

        GDALDataset* dstDs = vecDriver->Create(
            outputPath.toStdString().c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        if (!dstDs) {
            GDALClose(srcDs);
            setError("Failed to create output vector file: " + outputPath);
            return false;
        }

        OGRSpatialReference* srs = nullptr;
        const char* projWkt = srcDs->GetProjectionRef();
        if (projWkt && projWkt[0] != '\0') {
            srs = new OGRSpatialReference(projWkt);
        }

        OGRLayer* layer = dstDs->CreateLayer("polygonized", srs, wkbPolygon, nullptr);
        if (!layer) {
            delete srs;
            GDALClose(dstDs);
            GDALClose(srcDs);
            setError("Failed to create output layer");
            return false;
        }

        OGRFieldDefn field("DN", OFTInteger);
        layer->CreateField(&field);

        reportProgress(0.4, "Polygonizing raster...");

        CPLErr err = GDALPolygonize(band, nullptr, reinterpret_cast<OGRLayerH>(layer), 0, nullptr, nullptr, nullptr);

        delete srs;
        GDALClose(dstDs);
        GDALClose(srcDs);

        if (err != CE_None) {
            setError("GDAL polygonize failed");
            return false;
        }

        reportProgress(1.0, "Raster to vector conversion complete.");
        return true;
    }

    bool vectorToRaster(const QString& inputPath, const QString& outputPath) {
        GDALDataset* srcDs = static_cast<GDALDataset*>(
            GDALOpenEx(inputPath.toStdString().c_str(), GDAL_OF_VECTOR,
                       nullptr, nullptr, nullptr));
        if (!srcDs) {
            setError("Failed to open input vector: " + inputPath);
            return false;
        }

        OGRLayer* layer = srcDs->GetLayer(0);
        if (!layer) {
            GDALClose(srcDs);
            setError("Failed to read vector layer from input");
            return false;
        }

        reportProgress(0.2, "Computing extent...");

        OGREnvelope envelope;
        layer->GetExtent(&envelope);

        // Default resolution: 1000 cells across the longest dimension
        double xRes = (envelope.MaxX - envelope.MinX) / 1000.0;
        double yRes = (envelope.MaxY - envelope.MinY) / 1000.0;
        double cellSize = std::max(xRes, yRes);
        if (cellSize <= 0.0)
            cellSize = 1.0;

        int cols = static_cast<int>(std::ceil((envelope.MaxX - envelope.MinX) / cellSize));
        int rows = static_cast<int>(std::ceil((envelope.MaxY - envelope.MinY) / cellSize));
        if (cols < 1) cols = 1;
        if (rows < 1) rows = 1;

        reportProgress(0.4, "Creating raster output...");

        GDALDriver* rasDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
        if (!rasDriver) {
            GDALClose(srcDs);
            setError("GTiff driver not available");
            return false;
        }

        GDALDataset* dstDs = rasDriver->Create(
            outputPath.toStdString().c_str(), cols, rows, 1, GDT_Float32, nullptr);
        if (!dstDs) {
            GDALClose(srcDs);
            setError("Failed to create output raster: " + outputPath);
            return false;
        }

        double geoTransform[6] = {
            envelope.MinX, cellSize, 0.0,
            envelope.MaxY, 0.0, -cellSize
        };
        dstDs->SetGeoTransform(geoTransform);

        const OGRSpatialReference* srs = layer->GetSpatialRef();
        if (srs) {
            char* wkt = nullptr;
            const_cast<OGRSpatialReference*>(srs)->exportToWkt(&wkt);
            if (wkt) {
                dstDs->SetProjection(wkt);
                CPLFree(wkt);
            }
        }

        // Initialize raster to nodata
        GDALRasterBand* band = dstDs->GetRasterBand(1);
        band->SetNoDataValue(0.0);
        double fillVal = 0.0;
        band->Fill(fillVal);

        reportProgress(0.6, "Rasterizing vector...");

        int bandList[] = {1};
        double burnValues[] = {1.0};
        OGRLayerH layerH = reinterpret_cast<OGRLayerH>(layer);

        CPLErr err = GDALRasterizeLayers(
            dstDs, 1, bandList, 1, &layerH,
            nullptr, nullptr, burnValues, nullptr, nullptr, nullptr);

        GDALClose(dstDs);
        GDALClose(srcDs);

        if (err != CE_None) {
            setError("GDAL rasterize failed");
            return false;
        }

        reportProgress(1.0, "Vector to raster conversion complete.");
        return true;
    }
};

REGISTER_MODULE(RasterVectorModule)

} // namespace aplaceholder
