#include "MapCanvas.h"
#include <QWheelEvent>
#include <QMouseEvent>
#include <cmath>
#include <algorithm>

namespace aplaceholder {

MapCanvas::MapCanvas(QWidget* parent)
    : QGraphicsView(parent)
{
    m_scene = new QGraphicsScene(this);
    setScene(m_scene);

    setRenderHint(QPainter::SmoothPixmapTransform);
    setDragMode(QGraphicsView::ScrollHandDrag);
    setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    setResizeAnchor(QGraphicsView::AnchorViewCenter);
    setViewportUpdateMode(QGraphicsView::SmartViewportUpdate);

    setBackgroundBrush(QBrush(QColor(40, 40, 40)));
}

MapCanvas::~MapCanvas() = default;

void MapCanvas::displayRaster(std::unique_ptr<Raster> raster)
{
    m_raster = std::move(raster);
    if (!m_raster)
        return;

    QImage img = rasterToImage(*m_raster);

    m_scene->clear();
    m_pixmapItem = m_scene->addPixmap(QPixmap::fromImage(img));
    m_scene->setSceneRect(m_pixmapItem->boundingRect());

    zoomToFit();
}

QImage MapCanvas::rasterToImage(const Raster& raster, int band)
{
    int w = raster.cols();
    int h = raster.rows();

    // If 3+ bands, create RGB composite
    if (raster.bands() >= 3) {
        QImage img(w, h, QImage::Format_RGB32);
        const auto& r = raster.data(0);
        const auto& g = raster.data(1);
        const auto& b = raster.data(2);

        // Compute min/max for each band for stretching
        auto stats_r = raster.computeStats(0);
        auto stats_g = raster.computeStats(1);
        auto stats_b = raster.computeStats(2);

        for (int y = 0; y < h; ++y) {
            auto* scanline = reinterpret_cast<QRgb*>(img.scanLine(y));
            for (int x = 0; x < w; ++x) {
                size_t idx = static_cast<size_t>(y) * w + x;

                auto stretch = [](double val, double mn, double mx) -> int {
                    if (mx <= mn) return 128;
                    double norm = (val - mn) / (mx - mn);
                    return std::clamp(static_cast<int>(norm * 255), 0, 255);
                };

                int rv = stretch(r[idx], stats_r.min, stats_r.max);
                int gv = stretch(g[idx], stats_g.min, stats_g.max);
                int bv = stretch(b[idx], stats_b.min, stats_b.max);
                scanline[x] = qRgb(rv, gv, bv);
            }
        }
        return img;
    }

    // Single band — grayscale with min-max stretch
    QImage img(w, h, QImage::Format_Grayscale8);
    auto stats = raster.computeStats(band);
    const auto& data = raster.data(band);

    for (int y = 0; y < h; ++y) {
        auto* scanline = img.scanLine(y);
        for (int x = 0; x < w; ++x) {
            size_t idx = static_cast<size_t>(y) * w + x;
            double val = data[idx];

            if (raster.hasNoData() && val == raster.noDataValue()) {
                scanline[x] = 0;
                continue;
            }

            double norm = (stats.max > stats.min)
                ? (val - stats.min) / (stats.max - stats.min)
                : 0.5;
            scanline[x] = static_cast<uint8_t>(std::clamp(norm * 255.0, 0.0, 255.0));
        }
    }
    return img;
}

void MapCanvas::zoomIn()
{
    scale(m_zoomFactor, m_zoomFactor);
}

void MapCanvas::zoomOut()
{
    scale(1.0 / m_zoomFactor, 1.0 / m_zoomFactor);
}

void MapCanvas::zoomToFit()
{
    if (m_pixmapItem)
        fitInView(m_pixmapItem, Qt::KeepAspectRatio);
}

void MapCanvas::wheelEvent(QWheelEvent* event)
{
    if (event->angleDelta().y() > 0)
        zoomIn();
    else
        zoomOut();
}

void MapCanvas::mouseMoveEvent(QMouseEvent* event)
{
    QGraphicsView::mouseMoveEvent(event);
    // TODO: show coordinates in status bar
}

} // namespace aplaceholder
