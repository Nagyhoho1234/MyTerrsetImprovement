#pragma once

#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <memory>
#include "Raster.h"

namespace aplaceholder {

class MapCanvas : public QGraphicsView {
    Q_OBJECT

public:
    explicit MapCanvas(QWidget* parent = nullptr);
    ~MapCanvas() override;

    void displayRaster(std::unique_ptr<Raster> raster);
    void zoomIn();
    void zoomOut();
    void zoomToFit();

protected:
    void wheelEvent(QWheelEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;

private:
    QImage rasterToImage(const Raster& raster, int band = 0);

    QGraphicsScene* m_scene = nullptr;
    QGraphicsPixmapItem* m_pixmapItem = nullptr;
    std::unique_ptr<Raster> m_raster;
    double m_zoomFactor = 1.15;
};

} // namespace aplaceholder
