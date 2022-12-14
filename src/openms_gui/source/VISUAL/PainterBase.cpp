// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      QMatrix rotationMatrix;
      rotationMatrix.translate(start.x(), start.y());
      rotationMatrix.rotate(angle);
      QPainterPath path = rotationMatrix.map(arrow_start);
      painter->drawPath(path);
      bounding_rect = bounding_rect.united(path.boundingRect());
      //painter->restore();
    }
    if (!arrow_end.isEmpty())
    {
      QMatrix rotationMatrix;
      rotationMatrix.translate(end.x(), end.y());
      rotationMatrix.rotate(angle + 180);
      QPainterPath path = rotationMatrix.map(arrow_end);
      painter->drawPath(path);
      bounding_rect = bounding_rect.united(path.boundingRect());
    }
    return bounding_rect;
  }

} // namespace OpenMS