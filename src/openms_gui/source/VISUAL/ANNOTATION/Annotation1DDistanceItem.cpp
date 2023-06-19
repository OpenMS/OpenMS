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
