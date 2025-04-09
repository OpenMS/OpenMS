// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/PainterBase.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cassert>

#include <QPainter>
#include <QPen>
#include <QTransform>

using namespace std;

namespace OpenMS
{

  ShapeIcon PainterBase::toShapeIcon(const String& icon)
  {
    if (icon == "diamond") return ShapeIcon::DIAMOND;
    if (icon == "square")  return ShapeIcon::SQUARE;
    if (icon == "circle")  return ShapeIcon::CIRCLE;
    if (icon == "triangle")  return ShapeIcon::TRIANGLE;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Shape must be one of 'diamond', 'square', 'circle', 'triangle'!", icon);
  }

  void PainterBase::drawDashedLine(const QPoint& from, const QPoint& to, QPainter* painter, const QColor& color)
  {
    QPen pen;
    QVector<qreal> dashes;
    dashes << 5 << 5 << 1 << 5;
    pen.setDashPattern(dashes);
    pen.setColor(color);
    painter->save();
    painter->setPen(pen);
    painter->drawLine(from, to);
    painter->restore();
  }

  void PainterBase::drawCross(const QPoint& pos, QPainter* painter, const int size)
  {
    const int half_size = size / 2;
    painter->drawLine(pos.x(), pos.y() - half_size, pos.x(), pos.y() + half_size);
    painter->drawLine(pos.x() - half_size, pos.y(), pos.x() + half_size, pos.y());
  }

  void PainterBase::drawCaret(const QPoint& caret, QPainter* painter, const int size)
  {
    const int half_size = size / 2;
    painter->drawLine(caret.x(), caret.y(), caret.x() + half_size, caret.y() + half_size);
    painter->drawLine(caret.x(), caret.y(), caret.x() - half_size, caret.y() + half_size);
  }

  void PainterBase::drawDiamond(const QPoint& center, QPainter* painter, const int size)
  {
    const int half_size = size / 2;
    const auto x = center.x();
    const auto y = center.y();
    painter->drawLine(x, y + half_size, x + half_size, y);
    painter->drawLine(x + half_size, y, x, y - half_size);
    painter->drawLine(x, y - half_size, x - half_size, y);
    painter->drawLine(x - half_size, y, x, y + half_size);
  }
  
  void PainterBase::drawIcon(const QPoint& pos, const QRgb& color, const ShapeIcon icon, Size s, QPainter& p)
  {
    p.save();
    p.setPen(color);
    p.setBrush(QBrush(QColor(color), Qt::SolidPattern));

    int s_half = (int)s / 2;

    switch (icon)
    {
      break; case ShapeIcon::DIAMOND:
        {
        QPolygon pol;
        pol.putPoints(0, 4, pos.x() + s_half, pos.y(), pos.x(), pos.y() + s_half, pos.x() - (int)s_half, pos.y(), pos.x(), pos.y() - (int)s_half);
        p.drawConvexPolygon(pol);
        }
      break; case ShapeIcon::SQUARE:
        {
          QPolygon pol;
          pol.putPoints(0, 4, pos.x() + s_half, pos.y() + s_half, pos.x() - s_half, pos.y() + s_half, pos.x() - s_half, pos.y() - s_half, pos.x() + s_half, pos.y() - s_half);
          p.drawConvexPolygon(pol);
        }
      break; case ShapeIcon::CIRCLE:
        {
          p.drawEllipse(QRectF(pos.x() - s_half, pos.y() - s_half, s, s));
        }
      break; case ShapeIcon::TRIANGLE:
        {
          QPolygon pol;
          pol.putPoints(0, 3, pos.x(), pos.y() + s_half, pos.x() + s_half, pos.y() - (int)s_half, pos.x() - (int)s_half, pos.y() - (int)s_half);
          p.drawConvexPolygon(pol);
        }
      break; default:
        assert(false); // should never be reached 
    }
    p.restore();
  }

  QPainterPath PainterBase::getOpenArrow(int arrow_width)
  { // arrow definition
    QPainterPath arrow;
    arrow.moveTo(QPointF(0, 0));
    arrow.lineTo(QPointF(-arrow_width, 4));
    arrow.moveTo(QPointF(0, 0));
    arrow.lineTo(QPointF(-arrow_width, -4));
    return arrow;
  }

  QPainterPath PainterBase::getClosedArrow(int arrow_width)
  { // arrow definition
    QPainterPath arrow;
    arrow.moveTo(QPointF(0, 0));
    arrow.lineTo(QPointF(-arrow_width, 4));
    arrow.lineTo(QPointF(-arrow_width, -4));
    arrow.closeSubpath();
    return arrow;
  }

  QRectF PainterBase::drawLineWithArrows(QPainter* painter, const QPen& pen, const QPoint& start, const QPoint& end, const QPainterPath& arrow_start, const QPainterPath& arrow_end)
  {
    painter->setPen(pen);

    auto line = QLineF(start, end);
    // angle of line
    qreal angle = -line.angle() + 180; // negate since angle() reports counter-clockwise; +180 since painter.rotate() is more intuitive then
    QRectF bounding_rect = QRectF(line.p1(), line.p2()).normalized();
    // draw the actual line
    painter->drawLine(line);

    //painter->save();
    // draw arrow heads
    if (!arrow_start.isEmpty())
    {
      //painter->translate(start);
      //painter->rotate(angle);
      QTransform rotationMatrix;
      rotationMatrix.translate(start.x(), start.y());
      rotationMatrix.rotate(angle);
      QPainterPath path = rotationMatrix.map(arrow_start);
      painter->drawPath(path);
      bounding_rect = bounding_rect.united(path.boundingRect());
      //painter->restore();
    }
    if (!arrow_end.isEmpty())
    {
      QTransform rotationMatrix;
      rotationMatrix.translate(end.x(), end.y());
      rotationMatrix.rotate(angle + 180);
      QPainterPath path = rotationMatrix.map(arrow_end);
      painter->drawPath(path);
      bounding_rect = bounding_rect.united(path.boundingRect());
    }
    return bounding_rect;
  }

} // namespace OpenMS