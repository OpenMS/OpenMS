// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit, Chris Bielow $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>

#include <OpenMS/VISUAL/Plot1DCanvas.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  namespace
  {
    Annotation1DDistanceItem p("test", {0, 0}, {0, 0});  
  }

  Annotation1DDistanceItem::Annotation1DDistanceItem(const QString& text, const PointXYType& start_point, const PointXYType& end_point, const bool swap_ends_if_negative) :
      Annotation1DItem(text), start_point_(start_point), end_point_(end_point)
  {
    if (swap_ends_if_negative && start_point_ > end_point_)
    {
      { // make sure the distance is positive when creating the distance item
        start_point_.swap(end_point_);
      }
    }
  }

  void Annotation1DDistanceItem::ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index)
  {
    canvas->pushIntoDataRange(start_point_, layer_index);
    canvas->pushIntoDataRange(end_point_, layer_index);
  }

  void Annotation1DDistanceItem::draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped)
  {
    // translate mz/intensity to pixel coordinates
    QPoint start_p, end_p;
    canvas->dataToWidget(start_point_, start_p, flipped);
    canvas->dataToWidget(end_point_, end_p, flipped);

    // draw arrow heads and the ends if they won't overlap
    const auto arrow = ((start_p - end_p).manhattanLength() > 10) ? Painter1DBase::getClosedArrow(4) : QPainterPath();
    auto line_bb = Painter1DBase::drawLineWithArrows(&painter, painter.pen(), start_p, end_p, arrow, arrow).toRect();

    // find out how much additional space is needed for the text:
    QRect text_bb = painter.boundingRect(QRect(), Qt::AlignCenter, text_);
    // place text in the center of the line's bb, but gravitate it up
    auto new_text_center = canvas->getGravitator().gravitateWith(line_bb.center(), {line_bb.width() / 2 + text_bb.width()/2, line_bb.height() / 2 + text_bb.height()/2});
    text_bb.translate(new_text_center - text_bb.center());
    painter.drawText(text_bb, Qt::AlignHCenter, text_);

    // draw ticks
    if (!ticks_.empty())
    {
      for (auto tick_xy : ticks_)
      {
        tick_xy = canvas->getGravitator().gravitateTo(tick_xy, start_point_); // move to same level as line
        QPoint tick_px;
        canvas->dataToWidget(tick_xy, tick_px, flipped);
        QPoint tick_px_start = canvas->getGravitator().gravitateWith(tick_px, {4, 4});
        QPoint tick_px_end = canvas->getGravitator().gravitateWith(tick_px, {-8, -8});
        painter.drawLine(tick_px_start, tick_px_end);
      }
    }

    // overall bounding box
    bounding_box_ = text_bb.united(line_bb);

    if (selected_)
    {
      drawBoundingBox_(painter);
    }
  }

  void Annotation1DDistanceItem::move(const PointXYType delta, const Gravitator& gr, const DimMapper<2>& /*dim_mapper*/)
  {
    start_point_ = gr.gravitateWith(start_point_, delta); // only change the gravity axis
    end_point_ = gr.gravitateWith(end_point_, delta);     // only change the gravity axis
  }

  double Annotation1DDistanceItem::getDistance() const
  {
    const auto delta = end_point_ - start_point_;
    const auto dist = std::sqrt(delta.getX() * delta.getX() + delta.getY() * delta.getY());
    return std::copysign(1, delta.getX() + delta.getY()) * dist;
  }

  void Annotation1DDistanceItem::setTicks(const std::vector<PointXYType>& ticks)
  {
    ticks_ = ticks;
  }
} // namespace OpenMS
