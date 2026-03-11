#pragma once

#include "Module.h"
#include <QWidget>
#include <QDialog>

namespace aplaceholder {

/// Custom chart renderer using QPainter — no external chart library needed.
/// Supports histogram, scatter, line, bar, and pie chart types.
class ChartWidget : public QWidget {
    Q_OBJECT
public:
    explicit ChartWidget(QWidget* parent = nullptr);

    void setChartResult(const ChartResult& result);
    const ChartResult& chartResult() const { return m_result; }

    QSize sizeHint() const override { return QSize(700, 500); }

protected:
    void paintEvent(QPaintEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;

private:
    void drawHistogram(QPainter& p, const QRect& area);
    void drawScatter(QPainter& p, const QRect& area);
    void drawLine(QPainter& p, const QRect& area);
    void drawBar(QPainter& p, const QRect& area);
    void drawPie(QPainter& p, const QRect& area);
    void drawAxes(QPainter& p, const QRect& area,
                  double xMin, double xMax, double yMin, double yMax);
    void drawGrid(QPainter& p, const QRect& area,
                  double xMin, double xMax, double yMin, double yMax);

    QPointF dataToPixel(double x, double y, const QRect& area,
                        double xMin, double xMax, double yMin, double yMax) const;

    ChartResult m_result;
    QPoint m_mousePos;
};

/// Convenience dialog that wraps a ChartWidget with export button
class ChartDialog : public QDialog {
    Q_OBJECT
public:
    ChartDialog(const ChartResult& result, QWidget* parent = nullptr);

private slots:
    void exportImage();
    void exportCsv();

private:
    ChartWidget* m_chart;
};

} // namespace aplaceholder
