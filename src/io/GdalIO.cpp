#include "GdalIO.h"
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <QFileInfo>

namespace aplaceholder {

void GdalIO::initialize()
{
    GDALAllRegister();
}

void GdalIO::cleanup()
{
    GDALDestroyDriverManager();
}

std::unique_ptr<Raster> GdalIO::read(const QString& path)
{
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpen(path.toUtf8().constData(), GA_ReadOnly));
    if (!ds)
        return nullptr;

    int cols = ds->GetRasterXSize();
    int rows = ds->GetRasterYSize();
    int bands = ds->GetRasterCount();

    auto raster = std::make_unique<Raster>(cols, rows, bands, DataType::Float64);

    // Geo-transform
    double gt[6];
    if (ds->GetGeoTransform(gt) == CE_None) {
        GeoTransform geoTr;
        geoTr.originX = gt[0];
        geoTr.pixelWidth = gt[1];
        geoTr.rotationX = gt[2];
        geoTr.originY = gt[3];
        geoTr.rotationY = gt[4];
        geoTr.pixelHeight = gt[5];
        raster->setGeoTransform(geoTr);
    }

    // Projection
    const char* wkt = ds->GetProjectionRef();
    if (wkt && wkt[0])
        raster->setProjection(QString::fromUtf8(wkt));

    // Read each band
    for (int b = 0; b < bands; ++b) {
        GDALRasterBand* band = ds->GetRasterBand(b + 1);

        // NoData
        int hasNoData = 0;
        double noData = band->GetNoDataValue(&hasNoData);
        if (hasNoData && b == 0) {
            raster->setNoDataValue(noData);
        }

        // Read pixels as Float64
        auto& data = raster->data(b);
        CPLErr err = band->RasterIO(GF_Read, 0, 0, cols, rows,
                                     data.data(), cols, rows,
                                     GDT_Float64, 0, 0);
        if (err != CE_None) {
            GDALClose(ds);
            return nullptr;
        }
    }

    GDALClose(ds);
    return raster;
}

std::unique_ptr<Raster> GdalIO::readMetadata(const QString& path)
{
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpen(path.toUtf8().constData(), GA_ReadOnly));
    if (!ds)
        return nullptr;

    int cols = ds->GetRasterXSize();
    int rows = ds->GetRasterYSize();
    int bands = ds->GetRasterCount();

    // Determine data type from first band
    DataType dt = DataType::Float64;
    if (bands > 0) {
        GDALDataType gdt = ds->GetRasterBand(1)->GetRasterDataType();
        switch (gdt) {
            case GDT_Byte:    dt = DataType::Byte; break;
            case GDT_Int16:   dt = DataType::Int16; break;
            case GDT_Int32:   dt = DataType::Int32; break;
            case GDT_UInt16:  dt = DataType::UInt16; break;
            case GDT_UInt32:  dt = DataType::UInt32; break;
            case GDT_Float32: dt = DataType::Float32; break;
            default:          dt = DataType::Float64; break;
        }
    }

    auto raster = std::make_unique<Raster>();
    // Set dimensions without allocating pixel data
    // We use a minimal constructor approach
    *raster = Raster(cols, rows, bands, dt);
    // Clear data to save memory — metadata only
    // (The Raster constructor allocates, so we'd optimize this later)

    double gt[6];
    if (ds->GetGeoTransform(gt) == CE_None) {
        GeoTransform geoTr;
        geoTr.originX = gt[0];
        geoTr.pixelWidth = gt[1];
        geoTr.rotationX = gt[2];
        geoTr.originY = gt[3];
        geoTr.rotationY = gt[4];
        geoTr.pixelHeight = gt[5];
        raster->setGeoTransform(geoTr);
    }

    const char* wkt = ds->GetProjectionRef();
    if (wkt && wkt[0])
        raster->setProjection(QString::fromUtf8(wkt));

    int hasNoData = 0;
    if (bands > 0) {
        double noData = ds->GetRasterBand(1)->GetNoDataValue(&hasNoData);
        if (hasNoData)
            raster->setNoDataValue(noData);
    }

    GDALClose(ds);
    return raster;
}

std::unique_ptr<Raster> GdalIO::readBand(const QString& path, int band)
{
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpen(path.toUtf8().constData(), GA_ReadOnly));
    if (!ds)
        return nullptr;

    int cols = ds->GetRasterXSize();
    int rows = ds->GetRasterYSize();

    if (band < 1 || band > ds->GetRasterCount()) {
        GDALClose(ds);
        return nullptr;
    }

    auto raster = std::make_unique<Raster>(cols, rows, 1, DataType::Float64);

    double gt[6];
    if (ds->GetGeoTransform(gt) == CE_None) {
        GeoTransform geoTr{gt[0], gt[3], gt[1], gt[5], gt[2], gt[4]};
        raster->setGeoTransform(geoTr);
    }

    const char* wkt = ds->GetProjectionRef();
    if (wkt && wkt[0])
        raster->setProjection(QString::fromUtf8(wkt));

    GDALRasterBand* gdalBand = ds->GetRasterBand(band);
    int hasNoData = 0;
    double noData = gdalBand->GetNoDataValue(&hasNoData);
    if (hasNoData)
        raster->setNoDataValue(noData);

    auto& data = raster->data(0);
    gdalBand->RasterIO(GF_Read, 0, 0, cols, rows,
                        data.data(), cols, rows, GDT_Float64, 0, 0);

    GDALClose(ds);
    return raster;
}

bool GdalIO::write(const Raster& raster, const QString& path, const QString& driver)
{
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName(driver.toUtf8().constData());
    if (!drv)
        return false;

    GDALDataType gdt;
    switch (raster.dataType()) {
        case DataType::Byte:    gdt = GDT_Byte; break;
        case DataType::Int16:   gdt = GDT_Int16; break;
        case DataType::Int32:   gdt = GDT_Int32; break;
        case DataType::UInt16:  gdt = GDT_UInt16; break;
        case DataType::UInt32:  gdt = GDT_UInt32; break;
        case DataType::Float32: gdt = GDT_Float32; break;
        default:                gdt = GDT_Float64; break;
    }

    GDALDataset* ds = drv->Create(path.toUtf8().constData(),
                                   raster.cols(), raster.rows(), raster.bands(),
                                   gdt, nullptr);
    if (!ds)
        return false;

    // Geo-transform
    const auto& gt = raster.geoTransform();
    double adfGT[6] = {gt.originX, gt.pixelWidth, gt.rotationX,
                        gt.originY, gt.rotationY, gt.pixelHeight};
    ds->SetGeoTransform(adfGT);

    // Projection
    if (!raster.projection().isEmpty())
        ds->SetProjection(raster.projection().toUtf8().constData());

    // Write bands
    for (int b = 0; b < raster.bands(); ++b) {
        GDALRasterBand* band = ds->GetRasterBand(b + 1);

        if (raster.hasNoData())
            band->SetNoDataValue(raster.noDataValue());

        const auto& data = raster.data(b);
        band->RasterIO(GF_Write, 0, 0, raster.cols(), raster.rows(),
                        const_cast<double*>(data.data()),
                        raster.cols(), raster.rows(), GDT_Float64, 0, 0);
    }

    GDALClose(ds);
    return true;
}

QStringList GdalIO::supportedReadFormats()
{
    QStringList formats;
    int count = GDALGetDriverCount();
    for (int i = 0; i < count; ++i) {
        GDALDriverH drv = GDALGetDriver(i);
        const char* meta = GDALGetMetadataItem(drv, GDAL_DCAP_RASTER, nullptr);
        if (meta) {
            formats.append(QString::fromUtf8(GDALGetDriverShortName(drv)));
        }
    }
    return formats;
}

QStringList GdalIO::supportedWriteFormats()
{
    QStringList formats;
    int count = GDALGetDriverCount();
    for (int i = 0; i < count; ++i) {
        GDALDriverH drv = GDALGetDriver(i);
        const char* raster = GDALGetMetadataItem(drv, GDAL_DCAP_RASTER, nullptr);
        const char* create = GDALGetMetadataItem(drv, GDAL_DCAP_CREATE, nullptr);
        if (raster && create) {
            formats.append(QString::fromUtf8(GDALGetDriverShortName(drv)));
        }
    }
    return formats;
}

QString GdalIO::detectDriver(const QString& path)
{
    QFileInfo fi(path);
    QString ext = fi.suffix().toLower();

    if (ext == "tif" || ext == "tiff") return "GTiff";
    if (ext == "rst") return "RST";
    if (ext == "img") return "HFA";
    if (ext == "nc")  return "netCDF";
    if (ext == "hdf") return "HDF4";
    if (ext == "png") return "PNG";
    if (ext == "jpg" || ext == "jpeg") return "JPEG";

    return "GTiff"; // default
}

} // namespace aplaceholder
