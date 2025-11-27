// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Johannes Junker, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

namespace OpenMS
{
  namespace
  {
    Annotation1DTextItem p({0, 0}, "test");
  }

  void Annotation1DTextItem::ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index)
  {
    canvas->pushIntoDataRange(position_, layer_index);
  }

  void Annotation1DTextItem::draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped)
  {
    // translate units to pixel coordinates
    QPoint pos_text;
    canvas->dataToWidget(position_, pos_text, flipped);

    // compute bounding box of text_item on the specified painter
    bounding_box_ = painter.boundingRect(QRectF(pos_text, pos_text), flags_, text_);

    painter.drawText(bounding_box_, flags_, text_);
    if (selected_)
    {
      drawBoundingBox_(painter);
    }
  }

  void Annotation1DTextItem::move(const PointXYType delta, const Gravitator& /*gr*/, const DimMapper<2>& /*dim_mapper*/)
  {
    position_ += delta;
  }
} // namespace OpenMS
