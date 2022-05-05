// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <vector>

namespace OpenMS
{
  /** @brief An annotation item which represents a measured distance between two peaks.
      @see Annotation1DItem
  */
  template<class DataPoint>
  class Annotation1DDistanceItem :
    public Annotation1DItem
  {

public:
    /// Constructor
    Annotation1DDistanceItem(const QString & text, const DataPoint& start_point, const DataPoint& end_point)
      : Annotation1DItem(text), start_point_(start_point), end_point_(end_point)
    {
    }
    /// Copy constructor
    Annotation1DDistanceItem(const Annotation1DDistanceItem & rhs) = default;
    /// Destructor
    ~Annotation1DDistanceItem() override = default;

    // Docu in base class
    void ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index) override
    {
      canvas->pushIntoDataRange(start_point_, layer_index);
      canvas->pushIntoDataRange(end_point_, layer_index);
    }

    // Docu in base class
    void draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped = false) override
    {
      // translate mz/intensity to pixel coordinates
      QPoint start_p, end_p;
      canvas->dataToWidget(canvas->getMapper().map(start_point_), start_p, flipped, true);
      canvas->dataToWidget(canvas->getMapper().map(end_point_), end_p, flipped, true);

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
        auto pos_xy = canvas->getMapper().map(start_point_);
        for (const auto tick : ticks_)
        {
          auto tick_xy = canvas->getMapper().map(tick);
          tick_xy = canvas->getGravitator().gravitateTo(tick_xy, pos_xy); // move to same level as line
          QPoint tick_px;
          canvas->dataToWidget(tick_xy, tick_px, flipped, true);
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

    void move(PointXYType delta, const Gravitator& gr, const DimMapper<2>& dim_mapper) override
    {
      auto new_start_xy = gr.gravitateWith(dim_mapper.map(start_point_), delta); // only change the gravity axis
      dim_mapper.fromXY(new_start_xy, start_point_);

      auto new_end_xy = gr.gravitateWith(dim_mapper.map(end_point_), delta); // only change the gravity axis
      dim_mapper.fromXY(new_end_xy, end_point_);
    }
                           
    /// Returns the start point
    const DataPoint& getStartPoint() const
    {
      return start_point_;
    }

    /// Returns the end point
    const DataPoint& getEndPoint() const
    {
      return end_point_;
    }

    /// Set tick lines for the distance item (the gravity dimension is ignored)
    void setTicks(const std::vector<DataPoint>& ticks)
    {
      ticks_ = ticks;
    }

    // Docu in base class
    Annotation1DItem* clone() const override
    {
      return new Annotation1DDistanceItem(*this);
    }

  protected:
    /// The start point of the measured distance line (in data coordinates)
    DataPoint start_point_;
    /// The end point of the measured distance line (in data coordinates)
    DataPoint end_point_;
    /// Additional tick lines for the distance item (the gravity dimension is ignored)
    std::vector<DataPoint> ticks_;

  };
} // namespace OpenMS

