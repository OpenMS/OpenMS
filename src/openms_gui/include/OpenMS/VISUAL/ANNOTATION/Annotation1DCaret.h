// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/Painter1DBase.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

#include <QPainter>
#include <QtGui/QColor>
#include <QStaticText> 

#include <vector>

namespace OpenMS
{
  /** @brief An annotation item which paints a set of carets on the canvas.

      Most useful to visualize (theoretical) isotope distributions (one caret per isotope position).
      Additionally, a text annotation can be provided.

      @see Annotation1DItem
  */
  template <class DataPoint>
  class Annotation1DCaret :
    public Annotation1DItem
  {
public:
    typedef std::vector<DataPoint> PositionsType;
    using PointType = DataPoint;

    /// Constructor
    Annotation1DCaret(const PositionsType& caret_positions, const QString& text, const QColor& color, const QColor& connection_line_color) :
        Annotation1DItem(text), caret_positions_(caret_positions), position_(caret_positions[0]), color_(color), connection_line_color_(connection_line_color)
    {
      st_.setText(text);
    }

    /// Copy constructor
    Annotation1DCaret(const Annotation1DCaret& rhs) = default;

    /// Destructor
    ~Annotation1DCaret() override = default;

    // Docu in base class
    void ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index) override
    {
      canvas->pushIntoDataRange(position_, layer_index);
    }


    // Docu in base class
    void draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped = false) override
    {
      painter.save();

      painter.setPen(color_);
      // translate mz/intensity to pixel coordinates
      QPoint position_widget, caret_position_widget;

      auto xy_pos = canvas->getMapper().map(position_);
      auto xy_1stcaret = canvas->getMapper().map(position_);
      canvas->dataToWidget(xy_pos, position_widget, flipped);
      canvas->dataToWidget(xy_1stcaret, caret_position_widget, flipped);

      // draw carets '^'
      for (const auto& pos : caret_positions_)
      {
        auto xy_pos_caret = canvas->getMapper().map(pos);
        QPoint caret;
        canvas->dataToWidget(xy_pos_caret, caret, flipped);
        Painter1DBase::drawCaret(caret, &painter);
      }

      // compute bounding box of text_item on the specified painter
      bounding_box_ = QRectF(position_widget, st_.size());

      // shift pos - annotation should be over peak or, if not possible, next to it
      double vertical_shift = bounding_box_.height() / 2 + 5;
      if (!flipped)
      {
        vertical_shift *= -1;
      }

      bounding_box_.translate(0.0, vertical_shift);

      if (flipped && bounding_box_.bottom() > canvas->height())
      {
        bounding_box_.moveBottom(canvas->height());
        bounding_box_.moveLeft(position_widget.x() + 5.0);
      }
      else if (!flipped && bounding_box_.top() < 0.0)
      {
        bounding_box_.moveTop(0.0);
        bounding_box_.moveLeft(position_widget.x() + 5.0);
      }
      // keep inside canvas
      if (bounding_box_.right() > canvas->width())
      {
        bounding_box_.moveRight(canvas->width());
      }

      // draw connection line between anchor point and current position if pixel coordinates differ significantly
      if ((position_widget - caret_position_widget).manhattanLength() > 2)
      {
        QPointF border_point = GUIHelpers::intersectionPoint(bounding_box_, caret_position_widget);
        if (bounding_box_.center() != border_point)
        {
          painter.save();
          painter.setPen(Qt::DashLine);
          painter.drawLine(caret_position_widget, border_point);
          painter.restore();
        }
      }

      painter.drawStaticText(bounding_box_.topLeft(), st_);

      if (selected_)
      {
        drawBoundingBox_(painter);
      }

      painter.restore();
    }

    // Docu in base class
    void move(const PointXYType delta, const Gravitator& /*gr*/, const DimMapper<2>& dim_mapper) override
    {
      auto xy_before = dim_mapper.map(position_);
      xy_before += delta;
      dim_mapper.fromXY(xy_before, position_);
    }

    /// Returns the positions of the lines (in unit coordinates)
    const PositionsType& getCaretPositions() const
    {
      return caret_positions_;
    }

    /// Sets the position of the label (in unit coordinates)
    void setPosition(const DataPoint& position)
    {
      position_ = position;
    }
    /// Returns the position of the annotated peak (in unit coordinates)
    const DataPoint& getPosition() const
    {
      return position_;
    }

    /// Set the color of the carets (color of text must be set using html)
    void setColor(const QColor& color)
    {
      color_ = color;
    }
    /// Returns the color of the carets
    const QColor& getColor() const
    {
      return color_;
    }

    /// The text to display (optional).
    /// Rendered using QStaticText, so HTML formatting is allowed.
    void setRichText(const QString& text)
    {
      st_.setText(text);
      text_ = text; // this is just to keep the base class consistent.. we don't really use text_
    }

    // Docu in base class
    Annotation1DItem* clone() const override
    {
      return new Annotation1DCaret(*this);
    }

  protected:
    /// The positions of points (in unit coordinates)
    /// Ensure positions are sorted by m/z axis (or equivalent) when assigning
    PositionsType caret_positions_;

    /// The position of the label (in unit coordinates)
    DataPoint position_;

    /// The color of the label
    QColor color_;

    /// The color of the (optional) dashed line connecting peak and label
    QColor connection_line_color_;

    /// Holds the (rich) text
    QStaticText st_;
  };
} // namespace OpenMS

