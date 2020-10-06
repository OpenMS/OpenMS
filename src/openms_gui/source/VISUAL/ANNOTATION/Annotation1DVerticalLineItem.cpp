// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DVerticalLineItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

#include <QtCore/QPoint>
#include <QtGui/QPainter>

using namespace std;

namespace OpenMS
{

  Annotation1DVerticalLineItem::Annotation1DVerticalLineItem(const QString & text, const PointType & start_point, const PointType & end_point) :
      Annotation1DItem(text),
      start_point_(start_point),
      end_point_(end_point)
  {
  }

  Annotation1DVerticalLineItem::Annotation1DVerticalLineItem(const Annotation1DVerticalLineItem & rhs) :
      Annotation1DItem(rhs)
  {
    start_point_ = rhs.getStartPoint();
    end_point_ = rhs.getEndPoint();
  }

  Annotation1DVerticalLineItem::~Annotation1DVerticalLineItem()
  {
  }

  void Annotation1DVerticalLineItem::draw(Spectrum1DCanvas * const canvas, QPainter & painter, bool flipped)
  {
    //translate mz/intensity to pixel coordinates
    QPoint start_p, end_p;
    canvas->dataToWidget(start_point_.getX(), start_point_.getY(), start_p, flipped, true);
    canvas->dataToWidget(end_point_.getX(), end_point_.getY(), end_p, flipped, true);

    // compute bounding box on the specified painter
    if (canvas->isMzToXAxis())
    {
      bounding_box_ = QRectF(QPointF(start_p.x(), start_p.y()), QPointF(end_p.x(), end_p.y() + 4));     // +4 for lower half of arrow heads
    }
    else
    {
      bounding_box_ = QRectF(QPointF(start_p.x() - 4, start_p.y()), QPointF(end_p.x(), end_p.y()));
    }

    // find out how much additional space is needed for the text:
    QRectF text_boundings = painter.boundingRect(QRectF(), Qt::AlignCenter, text_);
    if (canvas->isMzToXAxis())
    {
      bounding_box_.setTop(bounding_box_.top() - text_boundings.height());
    }
    else
    {
      bounding_box_.setRight(bounding_box_.right() + text_boundings.width());
    }
    // if text doesn't fit between peaks, enlarge bounding box:
    if (canvas->isMzToXAxis())
    {
      if (text_boundings.width() > bounding_box_.width())
      {
        float additional_space = (text_boundings.width() - bounding_box_.width()) / 2;
        bounding_box_.setLeft(bounding_box_.left() - additional_space);
        bounding_box_.setRight(bounding_box_.right() + additional_space);
      }
    }
    else
    {
      if (text_boundings.height() > bounding_box_.height())
      {
        float additional_space = (text_boundings.height() - bounding_box_.height()) / 2;
        bounding_box_.setTop(bounding_box_.top() - additional_space);
        bounding_box_.setBottom(bounding_box_.bottom() + additional_space);
      }
    }

    // draw line
    painter.drawLine(start_p, end_p);

    // draw ticks
    if (!ticks_.empty())
    {
      for (vector<double>::iterator it = ticks_.begin(); it != ticks_.end(); ++it)
      {
        QPoint tick;
        canvas->dataToWidget(*it, start_point_.getY(), tick, flipped, true);
        painter.drawLine(tick.x(), tick.y() - 4, tick.x(), tick.y() + 4);
      }
    }

    // draw arrow heads and the ends if they won't overlap
    if ((start_p - end_p).manhattanLength() > 10)
    {
      if (canvas->isMzToXAxis())
      {
        painter.drawLine(start_p, QPoint(start_p.x() + 5, start_p.y() - 4));
        painter.drawLine(start_p, QPoint(start_p.x() + 5, start_p.y() + 4));
        painter.drawLine(end_p, QPoint(end_p.x() - 5, end_p.y() - 4));
        painter.drawLine(end_p, QPoint(end_p.x() - 5, end_p.y() + 4));
      }
      else
      {
        painter.drawLine(start_p, QPoint(start_p.x() + 4, start_p.y() - 5));
        painter.drawLine(start_p, QPoint(start_p.x() - 4, start_p.y() - 5));
        painter.drawLine(end_p, QPoint(end_p.x() + 4, end_p.y() + 5));
        painter.drawLine(end_p, QPoint(end_p.x() - 4, end_p.y() + 5));
      }
    }

    if (!canvas->isMzToXAxis())
    {
      bounding_box_.setWidth(bounding_box_.width() + 10.0);
    }

    painter.drawText(bounding_box_, Qt::AlignHCenter, text_);

    if (selected_)
    {
      drawBoundingBox_(painter);
    }
  }

  void Annotation1DVerticalLineItem::move(const PointType & delta)
  {
    // shift vertical position according to y-component of delta
    start_point_.setY(start_point_.getY() + delta.getY());
    end_point_.setY(end_point_.getY() + delta.getY());
  }

  void Annotation1DVerticalLineItem::setStartPoint(const PointType & p)
  {
    start_point_ = p;
  }

  void Annotation1DVerticalLineItem::setEndPoint(const PointType & p)
  {
    end_point_ = p;
  }

  const Annotation1DVerticalLineItem::PointType & Annotation1DVerticalLineItem::getStartPoint() const
  {
    return start_point_;
  }

  const Annotation1DVerticalLineItem::PointType & Annotation1DVerticalLineItem::getEndPoint() const
  {
    return end_point_;
  }

  void Annotation1DVerticalLineItem::setTicks(const std::vector<double> & ticks)
  {
    ticks_ = ticks;
  }

  void Annotation1DVerticalLineItem::ensureWithinDataRange(Spectrum1DCanvas * const canvas)
  {
    // can only be moved vertically, so check only y-position
    DRange<3> data_range = canvas->getDataRange();
    CoordinateType y_pos = start_point_.getY() * canvas->getPercentageFactor();

    if (y_pos < data_range.minPosition()[1])
    {
      start_point_.setY(data_range.minPosition()[1] / canvas->getPercentageFactor());
      end_point_.setY(start_point_.getY());
    }
    if (y_pos > data_range.maxPosition()[1])
    {
      start_point_.setY(data_range.maxPosition()[1] / canvas->getPercentageFactor());
      end_point_.setY(start_point_.getY());
    }
  }

} //Namespace
