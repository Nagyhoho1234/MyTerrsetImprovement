#pragma once

#include <QString>
#include <memory>
#include "Raster.h"

namespace aplaceholder {

// GDAL-based raster I/O — supports 160+ formats including IDRISI .rst
class GdalIO {
public:
    // Initialize GDAL (call once at startup)
    static void initialize();
    static void cleanup();

    // Read raster from any GDAL-supported format
    static std::unique_ptr<Raster> read(const QString& path);

    // Read only metadata (no pixel data loaded)
    static std::unique_ptr<Raster> readMetadata(const QString& path);

    // Read a single band
    static std::unique_ptr<Raster> readBand(const QString& path, int band);

    // Write raster to file
    // Supported drivers: "GTiff", "RST" (IDRISI), "HFA", "ENVI", "netCDF", etc.
    static bool write(const Raster& raster, const QString& path,
                      const QString& driver = "GTiff");

    // Get list of supported formats
    static QStringList supportedReadFormats();
    static QStringList supportedWriteFormats();

    // File format detection
    static QString detectDriver(const QString& path);
};

} // namespace aplaceholder
