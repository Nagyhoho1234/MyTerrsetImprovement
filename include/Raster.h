#pragma once

#include <QString>
#include <QSize>
#include <vector>
#include <memory>
#include <cstdint>

namespace aplaceholder {

enum class DataType {
    Byte,       // uint8
    Int16,      // int16
    Int32,      // int32
    Float32,    // float
    Float64,    // double  <-- TerrSet doesn't support this!
    UInt16,
    UInt32
};

struct GeoTransform {
    double originX = 0.0;
    double originY = 0.0;
    double pixelWidth = 1.0;
    double pixelHeight = 1.0;
    double rotationX = 0.0;
    double rotationY = 0.0;
};

class Raster {
public:
    Raster();
    Raster(int cols, int rows, int bands = 1, DataType type = DataType::Float64);
    ~Raster();

    // Dimensions
    int cols() const { return m_cols; }
    int rows() const { return m_rows; }
    int bands() const { return m_bands; }
    int64_t cellCount() const { return static_cast<int64_t>(m_cols) * m_rows; }
    DataType dataType() const { return m_dataType; }

    // Geo-referencing
    const GeoTransform& geoTransform() const { return m_geoTransform; }
    void setGeoTransform(const GeoTransform& gt) { m_geoTransform = gt; }
    const QString& projection() const { return m_projection; }
    void setProjection(const QString& proj) { m_projection = proj; }

    // NoData
    double noDataValue() const { return m_noData; }
    void setNoDataValue(double val) { m_noData = val; m_hasNoData = true; }
    bool hasNoData() const { return m_hasNoData; }

    // Pixel access (band index is 0-based)
    double value(int col, int row, int band = 0) const;
    void setValue(int col, int row, double val, int band = 0);

    // Bulk data access
    const std::vector<double>& data(int band = 0) const;
    std::vector<double>& data(int band = 0);

    // Allocate memory
    void allocate();
    bool isAllocated() const { return !m_data.empty(); }

    // Statistics
    struct Stats {
        double min = 0, max = 0, mean = 0, stddev = 0;
        int64_t validCount = 0;
    };
    Stats computeStats(int band = 0) const;

    // Coordinate transforms
    void colRowToXY(int col, int row, double& x, double& y) const;
    void xyToColRow(double x, double y, int& col, int& row) const;

private:
    int m_cols = 0;
    int m_rows = 0;
    int m_bands = 1;
    DataType m_dataType = DataType::Float64;
    GeoTransform m_geoTransform;
    QString m_projection;
    double m_noData = -9999.0;
    bool m_hasNoData = false;
    std::vector<std::vector<double>> m_data;  // m_data[band][row * cols + col]
};

} // namespace aplaceholder
