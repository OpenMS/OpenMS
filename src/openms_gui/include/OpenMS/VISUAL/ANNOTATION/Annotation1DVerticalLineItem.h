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
