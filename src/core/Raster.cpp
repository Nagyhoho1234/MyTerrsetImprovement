#include "Raster.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace aplaceholder {

Raster::Raster() = default;

Raster::Raster(int cols, int rows, int bands, DataType type)
    : m_cols(cols), m_rows(rows), m_bands(bands), m_dataType(type)
{
    allocate();
}

Raster::~Raster() = default;

void Raster::allocate()
{
    m_data.resize(m_bands);
    for (int b = 0; b < m_bands; ++b) {
        m_data[b].resize(static_cast<size_t>(m_cols) * m_rows, 0.0);
    }
}

double Raster::value(int col, int row, int band) const
{
    if (band < 0 || band >= m_bands || col < 0 || col >= m_cols || row < 0 || row >= m_rows)
        return m_noData;
    return m_data[band][static_cast<size_t>(row) * m_cols + col];
}

void Raster::setValue(int col, int row, double val, int band)
{
    if (band < 0 || band >= m_bands || col < 0 || col >= m_cols || row < 0 || row >= m_rows)
        return;
    m_data[band][static_cast<size_t>(row) * m_cols + col] = val;
}

const std::vector<double>& Raster::data(int band) const
{
    return m_data.at(band);
}

std::vector<double>& Raster::data(int band)
{
    return m_data.at(band);
}

Raster::Stats Raster::computeStats(int band) const
{
    Stats s;
    if (band < 0 || band >= m_bands || m_data.empty())
        return s;

    const auto& d = m_data[band];
    double sum = 0, sumSq = 0;
    s.min = std::numeric_limits<double>::max();
    s.max = std::numeric_limits<double>::lowest();
    s.validCount = 0;

    for (size_t i = 0; i < d.size(); ++i) {
        double v = d[i];
        if (m_hasNoData && v == m_noData)
            continue;
        if (std::isnan(v))
            continue;
        s.validCount++;
        sum += v;
        sumSq += v * v;
        if (v < s.min) s.min = v;
        if (v > s.max) s.max = v;
    }

    if (s.validCount > 0) {
        s.mean = sum / s.validCount;
        double variance = (sumSq / s.validCount) - (s.mean * s.mean);
        s.stddev = std::sqrt(std::max(0.0, variance));
    }
    return s;
}

void Raster::colRowToXY(int col, int row, double& x, double& y) const
{
    x = m_geoTransform.originX + col * m_geoTransform.pixelWidth
        + row * m_geoTransform.rotationX;
    y = m_geoTransform.originY + col * m_geoTransform.rotationY
        + row * m_geoTransform.pixelHeight;
}

void Raster::xyToColRow(double x, double y, int& col, int& row) const
{
    // Assumes no rotation
    col = static_cast<int>((x - m_geoTransform.originX) / m_geoTransform.pixelWidth);
    row = static_cast<int>((y - m_geoTransform.originY) / m_geoTransform.pixelHeight);
}

} // namespace aplaceholder
