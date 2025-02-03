// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim, Timo Sachsenberg $
// $Authors: Jihyung Kim, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>

#include <QtGui/QColor>

namespace OpenMS
{
  /** @brief An annotation item which represents a vertical line (or more precisely a line along the gravity axis, i.e. it could also be horizontal) and text label on top.
      @see Annotation1DItem
  */
  class Annotation1DVerticalLineItem :
      public Annotation1DItem
  {
  public:
    /**
      Constructor for a single vertical line of 1px width.

      @param center_pos Center of the line in unit coordinates (only the non-gravity component will be used)
      @param color Optional color. If invalid (=default), the current painter color will be used when this is painted
      @param text Optional text displayed next to the line. Can contain '\n' which will force multiple lines.
    **/ 
    Annotation1DVerticalLineItem(const PointXYType& center_pos, const QColor& color = QColor("as_before"), const QString& text = "");
    /**
      Constructor for a single vertical line of 1px width or a broader line (band) with the given width

      @param center_pos Center of the line in unit coordinates (only the non-gravity component will be used)
      @param width Full width of the band in unit coordinates; use =0 to make a thin line and not a band;
      @param alpha255 A transparency value from 0 (not visible), to 255 (fully opaque)
      @param dashed_line Should the line/band be dashed
      @param color Optional color. If invalid (=default), the current painter color will be used when this is painted
      @param text Optional text displayed next to the line/band. Can contain '\n' which will force multiple lines. Text will be plotted at the very top (modify using setTextYOffset())
    **/
    Annotation1DVerticalLineItem(const PointXYType& center_pos, const float width, const int alpha255 = 128, const bool dashed_line = false, const QColor& color = QColor("as_before"),
                                 const QString& text = "");

    /// Copy constructor
    Annotation1DVerticalLineItem(const Annotation1DVerticalLineItem& rhs) = default;

    /// Destructor
    ~Annotation1DVerticalLineItem() override = default;

    // Docu in base class
    void ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index) override;

    // Docu in base class
    void draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped = false) override;

    // Docu in base class
    void move(const PointXYType delta, const Gravitator& gr, const DimMapper<2>& dim_mapper) override;

    /// Sets the center position of the line (the widths will extend from there)
    void setPosition(const PointXYType& pos);

    /// Returns the position on the non-gravity-axis
    const PointXYType& getPosition() const;

    /// size of the painted text (width and height of the rectangle)
    QRectF getTextRect() const;

    /// offset the text by this much downwards in y-direction (to avoid overlaps etc)
    void setTextOffset(int offset);

    // Docu in base class
    Annotation1DItem* clone() const override
    {
      return new Annotation1DVerticalLineItem(*this);
    }
  protected:
    /// The position of the line (gravity axis is ignored)
    PointXYType pos_;
    /// offset (in pixel coordinates of gravity axis) in for the text (to avoid overlaps)
    int text_offset_ {0};
    /// width of the item in unit coordinates (allowing to show a 'band'; use =0 to make a thin line and not a band)
    float width_ = 0;
    /// transparency 0...255 of the band/line
    int alpha255_ = 128;
    /// is the band/line dashed?
    bool dashed_{false};
    /// The color of the line; if invalid, the current painter color will be used
    QColor color_ = Qt::black;
  };
} // namespace OpenMS
