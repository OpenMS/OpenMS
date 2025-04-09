// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim, Timo Sachsenberg  $
// $Authors: Jihyung Kim, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DVerticalLineItem.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <QtCore/QPoint>
#include <QtGui/QPainter>

using namespace std;

namespace OpenMS
{
  Annotation1DVerticalLineItem::Annotation1DVerticalLineItem(const PointXYType& center_pos, const QColor& color, const QString& text) :
      Annotation1DItem(text), pos_(center_pos),
      color_(color)
  {
  }

  Annotation1DVerticalLineItem::Annotation1DVerticalLineItem(const PointXYType& center_pos, const float width, const int alpha255, const bool dashed_line, const QColor& color, const QString& text)
    : Annotation1DItem(text), pos_(center_pos),
      width_(width),
      alpha255_(alpha255),
      dashed_(dashed_line),
      color_(color)
  {
  }

  void Annotation1DVerticalLineItem::ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index)
  {
    canvas->pushIntoDataRange(pos_, layer_index);
  }

  void Annotation1DVerticalLineItem::draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped)
  {
    painter.save();
    auto pen = painter.pen();

    QColor col = pen.color();
    if (color_.isValid())
    {
      col = color_;
    }
    col.setAlpha(alpha255_);
    // if you try this for larger widths, the dash pattern will scale up automatically (which might look ugly)
    // (if you now think of using { 5 / width_, ... }, to counter the scaling: this only works if width_ < 5, since internally Qt seems to use integer arithmetic...)
    if (dashed_)
    {
      pen.setDashPattern({ 5, 5, 1, 5 });
    }

    // get left/right corner points of the rectangle (line + width); names are as if the line is vertical, but depending on gravity, it could be horizontal as well
    QPoint start_px_left;
    canvas->dataToWidget(pos_, start_px_left, flipped);
    start_px_left = canvas->getGravitator().gravitateMax(start_px_left, canvas->canvasPixelArea());
    QPoint end_px_right = canvas->getGravitator().gravitateMin(start_px_left, canvas->canvasPixelArea());
    QPoint px_width;
    canvas->dataToWidgetDistance(width_, width_, px_width);
    px_width = canvas->getGravitator().gravitateZero(px_width); // make sure that 'height' is 0
    // get width in NON-gravity (=swapped) dimension
    int width = canvas->getGravitator().swap().gravityValue(px_width);
    // if width_==0, only make a 1px line; in any case, it should be 1 px at least
    pen.setWidth(std::max(1, width));
    pen.setColor(col);
    painter.setPen(pen);
    painter.drawLine(start_px_left, end_px_right);

    // compute bounding box on the specified painter
    bounding_box_ = QRectF(QPointF(start_px_left - px_width/2), QPointF(end_px_right + px_width/2)).normalized();

    if (!text_.isEmpty())
    {
      auto top_left_px = bounding_box_.topLeft() + QPoint(5,5);
      // shift gravity axis by text_offset_
      auto final = canvas->getGravitator().gravitateTo(top_left_px.toPoint(), QPoint(text_offset_, text_offset_));
      GUIHelpers::drawText(painter, text_.split('\n'), final, Qt::black);
    }

    painter.restore();
  }

  void Annotation1DVerticalLineItem::move(PointXYType delta, const Gravitator& gr, const DimMapper<2>& /*dim_mapper*/)
  {
    pos_ = gr.swap().gravitateWith(pos_, delta); // only change the non-gravity axis
  }

  void Annotation1DVerticalLineItem::setPosition(const PointXYType& pos)
  {
    pos_ = pos;
  }

  const PointXYType& Annotation1DVerticalLineItem::getPosition() const
  {
    return pos_;
  }

  QRectF Annotation1DVerticalLineItem::getTextRect() const
  {
    int dummy;
    return GUIHelpers::getTextDimension(getText().split('\n'), QFont("Courier"), dummy);
  }

  void Annotation1DVerticalLineItem::setTextOffset(int offset)
  {
    text_offset_ = offset;
  }
} //Namespace
