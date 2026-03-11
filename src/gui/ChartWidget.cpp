#include "ChartWidget.h"
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox>
#include <QToolTip>
#include <QFile>
#include <QTextStream>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace aplaceholder {

// Convert ChartColor to QColor
static QColor toQColor(const ChartColor& c) { return QColor(c.r, c.g, c.b); }

// ─── Default palette for multiple series ───
static const QColor kPalette[] = {
    QColor(31, 119, 180),   // blue
    QColor(255, 127, 14),   // orange
    QColor(44, 160, 44),    // green
    QColor(214, 39, 40),    // red
    QColor(148, 103, 189),  // purple
    QColor(140, 86, 75),    // brown
    QColor(227, 119, 194),  // pink
    QColor(127, 127, 127),  // gray
    QColor(188, 189, 34),   // olive
    QColor(23, 190, 207),   // cyan
};
static const int kPaletteSize = 10;

static QColor seriesColor(const ChartSeries& s, size_t index) {
    QColor c = toQColor(s.color);
    // If it's the default blue, assign from palette by index
    if (c == kPalette[0] && index > 0)
        return kPalette[index % kPaletteSize];
    return c;
}

// ─── ChartWidget ───

ChartWidget::ChartWidget(QWidget* parent)
    : QWidget(parent)
{
    setMouseTracking(true);
    setMinimumSize(400, 300);
}

void ChartWidget::setChartResult(const ChartResult& result)
{
    m_result = result;
    update();
}

void ChartWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    // Background
    p.fillRect(rect(), Qt::white);

    if (!m_result.hasData()) {
        p.drawText(rect(), Qt::AlignCenter, "No chart data");
        return;
    }

    // Layout: margins for title, labels, axes
    int leftMargin = 70;
    int rightMargin = 20;
    int topMargin = 40;
    int bottomMargin = 50;

    // Extra right margin for legend if multiple series
    if (m_result.series.size() > 1)
        rightMargin = 150;

    QRect plotArea(leftMargin, topMargin,
                   width() - leftMargin - rightMargin,
                   height() - topMargin - bottomMargin);

    // Title
    p.setPen(Qt::black);
    QFont titleFont = font();
    titleFont.setPointSize(12);
    titleFont.setBold(true);
    p.setFont(titleFont);
    p.drawText(QRect(0, 4, width(), topMargin - 4), Qt::AlignCenter, m_result.title);

    // Reset font
    QFont normalFont = font();
    normalFont.setPointSize(9);
    p.setFont(normalFont);

    // Draw based on type
    switch (m_result.type) {
    case ChartResult::Histogram: drawHistogram(p, plotArea); break;
    case ChartResult::Scatter:   drawScatter(p, plotArea);   break;
    case ChartResult::Line:      drawLine(p, plotArea);      break;
    case ChartResult::Bar:       drawBar(p, plotArea);       break;
    case ChartResult::Pie:       drawPie(p, plotArea);       break;
    default: break;
    }

    // Legend (if multiple series)
    if (m_result.series.size() > 1) {
        int legendX = plotArea.right() + 15;
        int legendY = plotArea.top() + 10;
        for (size_t i = 0; i < m_result.series.size(); i++) {
            QColor c = seriesColor(m_result.series[i], i);
            QRect colorRect(legendX, legendY + (int)i * 20, 12, 12);
            p.fillRect(colorRect, c);
            p.setPen(Qt::black);
            p.drawRect(colorRect);
            p.drawText(legendX + 18, legendY + (int)i * 20 + 11,
                       m_result.series[i].label);
        }
    }
}

void ChartWidget::mouseMoveEvent(QMouseEvent* event)
{
    m_mousePos = event->pos();
    QWidget::mouseMoveEvent(event);
}

// ─── Data-to-pixel mapping ───

QPointF ChartWidget::dataToPixel(double x, double y, const QRect& area,
                                  double xMin, double xMax, double yMin, double yMax) const
{
    double px = area.left() + (x - xMin) / (xMax - xMin) * area.width();
    double py = area.bottom() - (y - yMin) / (yMax - yMin) * area.height();
    return QPointF(px, py);
}

// ─── Axes and grid ───

static void niceRange(double& lo, double& hi, int& ticks)
{
    double range = hi - lo;
    if (range < 1e-12) { lo -= 1; hi += 1; range = 2; }
    double rough = range / ticks;
    double mag = std::pow(10.0, std::floor(std::log10(rough)));
    double frac = rough / mag;
    double nice = (frac <= 1.5) ? 1 : (frac <= 3) ? 2 : (frac <= 7) ? 5 : 10;
    double step = nice * mag;
    lo = std::floor(lo / step) * step;
    hi = std::ceil(hi / step) * step;
    ticks = std::max(1, (int)std::round((hi - lo) / step));
}

void ChartWidget::drawGrid(QPainter& p, const QRect& area,
                            double xMin, double xMax, double yMin, double yMax)
{
    QPen gridPen(QColor(230, 230, 230), 1, Qt::DotLine);
    p.setPen(gridPen);

    int xTicks = 8, yTicks = 6;
    niceRange(xMin, xMax, xTicks);
    niceRange(yMin, yMax, yTicks);

    double xStep = (xMax - xMin) / xTicks;
    double yStep = (yMax - yMin) / yTicks;

    for (int i = 1; i < xTicks; i++) {
        double x = xMin + i * xStep;
        QPointF pt = dataToPixel(x, yMin, area, xMin, xMax, yMin, yMax);
        p.drawLine(QPointF(pt.x(), area.top()), QPointF(pt.x(), area.bottom()));
    }
    for (int i = 1; i < yTicks; i++) {
        double y = yMin + i * yStep;
        QPointF pt = dataToPixel(xMin, y, area, xMin, xMax, yMin, yMax);
        p.drawLine(QPointF(area.left(), pt.y()), QPointF(area.right(), pt.y()));
    }
}

void ChartWidget::drawAxes(QPainter& p, const QRect& area,
                            double xMin, double xMax, double yMin, double yMax)
{
    drawGrid(p, area, xMin, xMax, yMin, yMax);

    // Axis lines
    QPen axisPen(Qt::black, 1.5);
    p.setPen(axisPen);
    p.drawLine(area.bottomLeft(), area.bottomRight());
    p.drawLine(area.topLeft(), area.bottomLeft());

    // Tick labels
    QFont tickFont = font();
    tickFont.setPointSize(8);
    p.setFont(tickFont);
    p.setPen(Qt::black);

    int xTicks = 8, yTicks = 6;
    niceRange(xMin, xMax, xTicks);
    niceRange(yMin, yMax, yTicks);
    double xStep = (xMax - xMin) / xTicks;
    double yStep = (yMax - yMin) / yTicks;

    // X axis ticks
    for (int i = 0; i <= xTicks; i++) {
        double x = xMin + i * xStep;
        QPointF pt = dataToPixel(x, yMin, area, xMin, xMax, yMin, yMax);
        p.drawLine(QPointF(pt.x(), area.bottom()), QPointF(pt.x(), area.bottom() + 4));
        QString label = (std::abs(x) < 0.01 && x != 0) ? QString::number(x, 'e', 1)
                       : QString::number(x, 'g', 5);
        p.drawText(QRectF(pt.x() - 30, area.bottom() + 5, 60, 15),
                   Qt::AlignCenter, label);
    }

    // Y axis ticks
    for (int i = 0; i <= yTicks; i++) {
        double y = yMin + i * yStep;
        QPointF pt = dataToPixel(xMin, y, area, xMin, xMax, yMin, yMax);
        p.drawLine(QPointF(area.left() - 4, pt.y()), QPointF(area.left(), pt.y()));
        QString label = (std::abs(y) < 0.01 && y != 0) ? QString::number(y, 'e', 1)
                       : QString::number(y, 'g', 5);
        p.drawText(QRectF(area.left() - 65, pt.y() - 8, 60, 16),
                   Qt::AlignRight | Qt::AlignVCenter, label);
    }

    // Axis labels
    QFont labelFont = font();
    labelFont.setPointSize(9);
    labelFont.setBold(true);
    p.setFont(labelFont);

    // X label
    p.drawText(QRect(area.left(), area.bottom() + 25, area.width(), 20),
               Qt::AlignCenter, m_result.xLabel);

    // Y label (rotated)
    p.save();
    p.translate(12, area.center().y());
    p.rotate(-90);
    p.drawText(QRect(-area.height()/2, 0, area.height(), 20),
               Qt::AlignCenter, m_result.yLabel);
    p.restore();
}

// ─── Histogram ───

void ChartWidget::drawHistogram(QPainter& p, const QRect& area)
{
    if (m_result.series.empty() || m_result.series[0].x.empty()) return;
    const auto& s = m_result.series[0];
    size_t n = s.x.size();

    double xMin = s.x.front(), xMax = s.x.back();
    double yMin = 0, yMax = *std::max_element(s.y.begin(), s.y.end());
    if (yMax <= 0) yMax = 1;
    yMax *= 1.05;

    double barW = (n > 1) ? (s.x[1] - s.x[0]) : 1.0;
    xMin -= barW * 0.5;
    xMax += barW * 0.5;

    int yT = 6;
    double axXMin = xMin, axXMax = xMax, axYMin = yMin, axYMax = yMax;
    niceRange(axYMin, axYMax, yT);

    drawAxes(p, area, axXMin, axXMax, axYMin, axYMax);

    QColor fillColor = seriesColor(s, 0);
    fillColor.setAlpha(180);
    p.setPen(QPen(fillColor.darker(130), 1));
    p.setBrush(fillColor);

    for (size_t i = 0; i < n; i++) {
        double left = s.x[i] - barW * 0.45;
        double right = s.x[i] + barW * 0.45;
        QPointF topLeft = dataToPixel(left, s.y[i], area, axXMin, axXMax, axYMin, axYMax);
        QPointF botRight = dataToPixel(right, 0, area, axXMin, axXMax, axYMin, axYMax);
        p.drawRect(QRectF(topLeft, botRight));
    }
}

// ─── Scatter ───

void ChartWidget::drawScatter(QPainter& p, const QRect& area)
{
    double xMin = 1e30, xMax = -1e30, yMin = 1e30, yMax = -1e30;
    for (const auto& s : m_result.series) {
        for (double v : s.x) { xMin = std::min(xMin, v); xMax = std::max(xMax, v); }
        for (double v : s.y) { yMin = std::min(yMin, v); yMax = std::max(yMax, v); }
    }
    double xPad = (xMax - xMin) * 0.05, yPad = (yMax - yMin) * 0.05;
    xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;

    int xT = 8, yT = 6;
    niceRange(xMin, xMax, xT);
    niceRange(yMin, yMax, yT);

    drawAxes(p, area, xMin, xMax, yMin, yMax);

    for (size_t si = 0; si < m_result.series.size(); si++) {
        const auto& s = m_result.series[si];
        QColor dotColor = seriesColor(s, si);
        dotColor.setAlpha(120);
        p.setPen(Qt::NoPen);
        p.setBrush(dotColor);
        size_t n = std::min(s.x.size(), s.y.size());
        for (size_t i = 0; i < n; i++) {
            QPointF pt = dataToPixel(s.x[i], s.y[i], area, xMin, xMax, yMin, yMax);
            p.drawEllipse(pt, 3, 3);
        }
    }
}

// ─── Line ───

void ChartWidget::drawLine(QPainter& p, const QRect& area)
{
    double xMin = 1e30, xMax = -1e30, yMin = 1e30, yMax = -1e30;
    for (const auto& s : m_result.series) {
        for (double v : s.x) { xMin = std::min(xMin, v); xMax = std::max(xMax, v); }
        for (double v : s.y) { yMin = std::min(yMin, v); yMax = std::max(yMax, v); }
    }
    double xPad = (xMax - xMin) * 0.05, yPad = (yMax - yMin) * 0.05;
    xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;

    int xT = 8, yT = 6;
    niceRange(xMin, xMax, xT);
    niceRange(yMin, yMax, yT);

    drawAxes(p, area, xMin, xMax, yMin, yMax);

    for (size_t si = 0; si < m_result.series.size(); si++) {
        const auto& s = m_result.series[si];
        QColor c = seriesColor(s, si);
        QPen linePen(c, 2);
        p.setPen(linePen);
        p.setBrush(Qt::NoBrush);

        size_t n = std::min(s.x.size(), s.y.size());
        if (n < 2) continue;

        QPainterPath path;
        QPointF first = dataToPixel(s.x[0], s.y[0], area, xMin, xMax, yMin, yMax);
        path.moveTo(first);
        for (size_t i = 1; i < n; i++) {
            QPointF pt = dataToPixel(s.x[i], s.y[i], area, xMin, xMax, yMin, yMax);
            path.lineTo(pt);
        }
        p.drawPath(path);

        // Draw data points
        p.setPen(Qt::NoPen);
        p.setBrush(c);
        for (size_t i = 0; i < n; i++) {
            QPointF pt = dataToPixel(s.x[i], s.y[i], area, xMin, xMax, yMin, yMax);
            p.drawEllipse(pt, 3, 3);
        }
    }
}

// ─── Bar chart ───

void ChartWidget::drawBar(QPainter& p, const QRect& area)
{
    if (m_result.series.empty()) return;
    const auto& s = m_result.series[0];
    size_t n = s.y.size();
    if (n == 0) return;

    double yMin = 0;
    double yMax = *std::max_element(s.y.begin(), s.y.end()) * 1.1;
    if (yMax <= 0) yMax = 1;

    int yT = 6;
    double axYMin = yMin, axYMax = yMax;
    niceRange(axYMin, axYMax, yT);

    QPen axisPen(Qt::black, 1.5);
    p.setPen(axisPen);
    p.drawLine(area.bottomLeft(), area.bottomRight());
    p.drawLine(area.topLeft(), area.bottomLeft());

    // Y ticks
    QFont tickFont = font();
    tickFont.setPointSize(8);
    p.setFont(tickFont);
    double yStep = (axYMax - axYMin) / yT;
    for (int i = 0; i <= yT; i++) {
        double y = axYMin + i * yStep;
        double py = area.bottom() - (y - axYMin) / (axYMax - axYMin) * area.height();
        p.setPen(Qt::black);
        p.drawLine(QPointF(area.left() - 4, py), QPointF(area.left(), py));
        p.drawText(QRectF(area.left() - 65, py - 8, 60, 16),
                   Qt::AlignRight | Qt::AlignVCenter, QString::number(y, 'g', 5));
        // Grid line
        QPen gp(QColor(230, 230, 230), 1, Qt::DotLine);
        p.setPen(gp);
        p.drawLine(QPointF(area.left(), py), QPointF(area.right(), py));
    }

    // Bars
    double barW = area.width() / (double)(n + 1);
    for (size_t i = 0; i < n; i++) {
        double x = area.left() + (i + 0.5) * barW;
        double barH = (s.y[i] / axYMax) * area.height();
        QRectF barRect(x, area.bottom() - barH, barW * 0.8, barH);

        QColor c = kPalette[i % kPaletteSize];
        c.setAlpha(200);
        p.setPen(QPen(c.darker(130), 1));
        p.setBrush(c);
        p.drawRect(barRect);

        // Category label
        QString label = (i < m_result.categoryLabels.size())
                        ? m_result.categoryLabels[i]
                        : QString::number(i + 1);
        p.setPen(Qt::black);
        p.save();
        p.translate(x + barW * 0.4, area.bottom() + 5);
        p.rotate(45);
        p.drawText(0, 0, 80, 15, Qt::AlignLeft, label);
        p.restore();
    }

    // Y label
    QFont labelFont = font();
    labelFont.setPointSize(9);
    labelFont.setBold(true);
    p.setFont(labelFont);
    p.setPen(Qt::black);
    p.save();
    p.translate(12, area.center().y());
    p.rotate(-90);
    p.drawText(QRect(-area.height()/2, 0, area.height(), 20),
               Qt::AlignCenter, m_result.yLabel);
    p.restore();
}

// ─── Pie chart ───

void ChartWidget::drawPie(QPainter& p, const QRect& area)
{
    if (m_result.series.empty()) return;
    const auto& s = m_result.series[0];
    size_t n = s.y.size();
    if (n == 0) return;

    double total = std::accumulate(s.y.begin(), s.y.end(), 0.0);
    if (total <= 0) return;

    int radius = std::min(area.width(), area.height()) / 2 - 20;
    QPoint center = area.center();
    QRect pieRect(center.x() - radius, center.y() - radius, radius * 2, radius * 2);

    double startAngle = 90 * 16;  // Start from top
    for (size_t i = 0; i < n; i++) {
        double span = (s.y[i] / total) * 360.0 * 16;
        QColor c = kPalette[i % kPaletteSize];
        p.setPen(QPen(Qt::white, 1));
        p.setBrush(c);
        p.drawPie(pieRect, (int)startAngle, (int)span);

        // Label
        double midAngle = (startAngle + span / 2) / 16.0;
        double rad = midAngle * M_PI / 180.0;
        double lx = center.x() + (radius + 20) * std::cos(rad);
        double ly = center.y() - (radius + 20) * std::sin(rad);
        p.setPen(Qt::black);
        QString label = (i < m_result.categoryLabels.size())
                        ? m_result.categoryLabels[i]
                        : QString::number(i + 1);
        double pct = (s.y[i] / total) * 100;
        p.drawText(QRectF(lx - 40, ly - 8, 80, 16), Qt::AlignCenter,
                   QString("%1 (%2%)").arg(label).arg(pct, 0, 'f', 1));

        startAngle += span;
    }
}

// ─── ChartDialog ───

ChartDialog::ChartDialog(const ChartResult& result, QWidget* parent)
    : QDialog(parent)
{
    setWindowTitle(result.title);
    resize(750, 550);

    auto* layout = new QVBoxLayout(this);

    m_chart = new ChartWidget(this);
    m_chart->setChartResult(result);
    layout->addWidget(m_chart, 1);

    // Buttons
    auto* btnLayout = new QHBoxLayout();
    btnLayout->addStretch();

    auto* exportImgBtn = new QPushButton("Export Image...");
    connect(exportImgBtn, &QPushButton::clicked, this, &ChartDialog::exportImage);
    btnLayout->addWidget(exportImgBtn);

    auto* exportCsvBtn = new QPushButton("Export CSV...");
    connect(exportCsvBtn, &QPushButton::clicked, this, &ChartDialog::exportCsv);
    btnLayout->addWidget(exportCsvBtn);

    auto* closeBtn = new QPushButton("Close");
    connect(closeBtn, &QPushButton::clicked, this, &QDialog::accept);
    btnLayout->addWidget(closeBtn);

    layout->addLayout(btnLayout);
}

void ChartDialog::exportImage()
{
    QString path = QFileDialog::getSaveFileName(this, "Export Chart Image",
        QString(), "PNG Images (*.png);;BMP Images (*.bmp);;JPEG Images (*.jpg)");
    if (path.isEmpty()) return;

    QPixmap pixmap(m_chart->size());
    m_chart->render(&pixmap);
    if (!pixmap.save(path)) {
        QMessageBox::warning(this, "Export Failed", "Could not save image to: " + path);
    }
}

void ChartDialog::exportCsv()
{
    QString path = QFileDialog::getSaveFileName(this, "Export Chart Data",
        QString(), "CSV Files (*.csv);;All Files (*)");
    if (path.isEmpty()) return;

    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QMessageBox::warning(this, "Export Failed", "Could not open file: " + path);
        return;
    }

    QTextStream out(&file);
    const auto& result = m_chart->chartResult();

    // Header
    for (size_t i = 0; i < result.series.size(); i++) {
        if (i > 0) out << ",";
        QString label = result.series[i].label.isEmpty()
                        ? QString("Series_%1").arg(i + 1) : result.series[i].label;
        out << label << "_X," << label << "_Y";
    }
    out << "\n";

    // Data rows
    size_t maxLen = 0;
    for (const auto& s : result.series)
        maxLen = std::max(maxLen, std::max(s.x.size(), s.y.size()));

    for (size_t row = 0; row < maxLen; row++) {
        for (size_t i = 0; i < result.series.size(); i++) {
            if (i > 0) out << ",";
            const auto& s = result.series[i];
            if (row < s.x.size()) out << s.x[row];
            out << ",";
            if (row < s.y.size()) out << s.y[row];
        }
        out << "\n";
    }
}

} // namespace aplaceholder
