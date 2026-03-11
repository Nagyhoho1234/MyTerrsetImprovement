#include "RasterInfoPanel.h"
#include <QFileInfo>

namespace aplaceholder {

RasterInfoPanel::RasterInfoPanel(QWidget* parent)
    : QWidget(parent)
{
    auto* layout = new QFormLayout(this);
    layout->setLabelAlignment(Qt::AlignRight);

    m_path = new QLabel("-");
    m_dims = new QLabel("-");
    m_bands = new QLabel("-");
    m_dataType = new QLabel("-");
    m_projection = new QLabel("-");
    m_pixelSize = new QLabel("-");
    m_extent = new QLabel("-");
    m_stats = new QLabel("-");

    m_projection->setWordWrap(true);

    layout->addRow("File:", m_path);
    layout->addRow("Dimensions:", m_dims);
    layout->addRow("Bands:", m_bands);
    layout->addRow("Data Type:", m_dataType);
    layout->addRow("Pixel Size:", m_pixelSize);
    layout->addRow("Extent:", m_extent);
    layout->addRow("Projection:", m_projection);
    layout->addRow("Statistics:", m_stats);
}

void RasterInfoPanel::showRasterInfo(const Raster& raster, const QString& path)
{
    QFileInfo fi(path);
    m_path->setText(fi.fileName());
    m_dims->setText(QString("%1 x %2 (%L3 cells)")
        .arg(raster.cols()).arg(raster.rows()).arg(raster.cellCount()));
    m_bands->setText(QString::number(raster.bands()));

    QString typeStr;
    switch (raster.dataType()) {
        case DataType::Byte:    typeStr = "Byte (uint8)"; break;
        case DataType::Int16:   typeStr = "Int16"; break;
        case DataType::Int32:   typeStr = "Int32"; break;
        case DataType::Float32: typeStr = "Float32"; break;
        case DataType::Float64: typeStr = "Float64 (double)"; break;
        case DataType::UInt16:  typeStr = "UInt16"; break;
        case DataType::UInt32:  typeStr = "UInt32"; break;
    }
    m_dataType->setText(typeStr);

    const auto& gt = raster.geoTransform();
    m_pixelSize->setText(QString("%1 x %2").arg(gt.pixelWidth).arg(std::abs(gt.pixelHeight)));

    double xmin = gt.originX;
    double ymax = gt.originY;
    double xmax = xmin + raster.cols() * gt.pixelWidth;
    double ymin = ymax + raster.rows() * gt.pixelHeight;
    m_extent->setText(QString("X: %1 to %2, Y: %3 to %4")
        .arg(xmin, 0, 'f', 2).arg(xmax, 0, 'f', 2)
        .arg(ymin, 0, 'f', 2).arg(ymax, 0, 'f', 2));

    QString proj = raster.projection();
    if (proj.length() > 80)
        proj = proj.left(80) + "...";
    m_projection->setText(proj.isEmpty() ? "Unknown" : proj);

    auto stats = raster.computeStats(0);
    m_stats->setText(QString("Min: %1  Max: %2  Mean: %3  StdDev: %4  (Band 1)")
        .arg(stats.min, 0, 'f', 4)
        .arg(stats.max, 0, 'f', 4)
        .arg(stats.mean, 0, 'f', 4)
        .arg(stats.stddev, 0, 'f', 4));
}

} // namespace aplaceholder
