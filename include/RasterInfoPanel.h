#pragma once

#include <QWidget>
#include <QLabel>
#include <QFormLayout>
#include "Raster.h"

namespace aplaceholder {

class RasterInfoPanel : public QWidget {
    Q_OBJECT

public:
    explicit RasterInfoPanel(QWidget* parent = nullptr);

    void showRasterInfo(const Raster& raster, const QString& path);

private:
    QLabel* m_path;
    QLabel* m_dims;
    QLabel* m_bands;
    QLabel* m_dataType;
    QLabel* m_projection;
    QLabel* m_pixelSize;
    QLabel* m_extent;
    QLabel* m_stats;
};

} // namespace aplaceholder
