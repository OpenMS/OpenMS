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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DCaret.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

#include <QtGui/QPainter>
#include <QtCore/QPoint>


namespace OpenMS
{

  Annotation1DCaret::Annotation1DCaret(const Annotation1DCaret::PositionsType& caret_positions, 
                                       const QString& text,
                                       const QColor& colour,
                                       const QColor& connection_line_color) :
    Annotation1DItem(text),
    caret_positions_(caret_positions),
    position_(caret_positions[0]),
    color_(colour),
    connection_line_color_(connection_line_color)
  {
    st_.setText(text);
  }

  Annotation1DCaret::Annotation1DCaret(const Annotation1DCaret& rhs) :
    Annotation1DItem(rhs)
  {
    caret_positions_ = rhs.caret_positions_;
    position_ = rhs.position_;
    color_ = rhs.color_;
    connection_line_color_ = rhs.connection_line_color_;
    st_ = rhs.st_;
  }

  Annotation1DCaret::~Annotation1DCaret()
  {
  }

  void Annotation1DCaret::setRichText(const QString& text)
  {
    st_.setText(text);
    text_ = text; // this is just to keep the base class consistent.. we don't really use text_
  }

  void Annotation1DCaret::draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped)
  {
    painter.save();

    painter.setPen(color_);
    // translate mz/intensity to pixel coordinates
    QPoint position_widget, caret_position_widget;

    canvas->dataToWidget(position_.getX(), position_.getY(), position_widget, flipped, true);
    canvas->dataToWidget(caret_positions_[0].getX(), caret_positions_[0].getY(), caret_position_widget, flipped, true);

    //std::cerr << "color" << color_.value() << " ";
    // draw ticks (for now)
    if (!caret_positions_.empty())
    {
      QPoint caret;        // draw ^ to indicate theoretical position
      for (PositionsType::iterator it = caret_positions_.begin(); it != caret_positions_.end(); ++it)
      {
        canvas->dataToWidget(it->getX(), it->getY(), caret, flipped, true);
        painter.drawLine(caret.x(), caret.y(), caret.x()+4, caret.y() + 4);
        painter.drawLine(caret.x(), caret.y(), caret.x()-4, caret.y() + 4);
        //std::cout << "caret: " << caret.x() << "," << caret.y() << "\n";
      }
    }

    // compute bounding box of text_item on the specified painter
    bounding_box_ = QRectF(position_widget, st_.size());
    //std::cout << "posP: " << position_.getX() << "," << position_.getY() << "\n";
    //std::cout << "posW: " << position_widget.x() << "," << position_widget.y() << "\n";
    //std::cout <<"init BB topleft: " << bounding_box_.topLeft().x()  << "," << bounding_box_.topLeft().y() <<"\n";

    double vertical_shift = 0;
    double horizontal_shift = 0;

    if (canvas->isMzToXAxis())
    {
      // shift pos - annotation should be over peak or, if not possible, next to it
      vertical_shift = bounding_box_.height() / 2 + 5;
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
    }
    else
    {
      // annotation should be next to the peak (to its right)
      horizontal_shift = bounding_box_.width() / 2 + 5;
      bounding_box_.translate(horizontal_shift, 0.0);
      if (bounding_box_.right() > canvas->width())
      {
        bounding_box_.moveRight(canvas->width());
      }
    }

    // draw connection line between anchor point and current position if pixel coordinates differ significantly
    if ((position_widget - caret_position_widget).manhattanLength() > 2)
    {
      // check if line crosses bounding box, if so move line startpoint to correct bounding box intersection
      QLineF line(caret_position_widget, position_widget + QPoint(horizontal_shift, vertical_shift));
      QLineF top(bounding_box_.x(), bounding_box_.y(), bounding_box_.x() + bounding_box_.width(), bounding_box_.y());
      QLineF left(bounding_box_.x(), bounding_box_.y(), bounding_box_.x(), bounding_box_.y() + bounding_box_.height());
      QLineF right(bounding_box_.x() + bounding_box_.width(), bounding_box_.y(), bounding_box_.x() + bounding_box_.width(), bounding_box_.y() + bounding_box_.height());
      QLineF bottom(bounding_box_.x(), bounding_box_.y() + bounding_box_.height(), bounding_box_.x() + bounding_box_.width(), bounding_box_.y() + bounding_box_.height());

      QLineF::IntersectType itype;
      QPointF * ip = new QPointF();
      QPointF * closest_ip = new QPointF(-10e10, -10e10);
      bool found_intersection = false;

      // intersection with top
      itype = line.intersect(top, ip);
      if (itype == QLineF::BoundedIntersection &&
          QLineF(caret_position_widget, *ip).length() < QLineF(caret_position_widget, *closest_ip).length())
      {
        found_intersection = true;
        *closest_ip = *ip;
      }
      // intersection with left
      itype = line.intersect(left, ip);
      if (itype == QLineF::BoundedIntersection &&
          QLineF(caret_position_widget, *ip).length() < QLineF(caret_position_widget, *closest_ip).length())
      {
        found_intersection = true;
        *closest_ip = *ip;
      }

      // intersection with right
      itype = line.intersect(right, ip);
      if (itype == QLineF::BoundedIntersection &&
          QLineF(caret_position_widget, *ip).length() < QLineF(caret_position_widget, *closest_ip).length())
      {
        found_intersection = true;
        *closest_ip = *ip;
      }

      // intersection with bottom
      itype = line.intersect(bottom, ip);
      if (itype == QLineF::BoundedIntersection &&
          QLineF(caret_position_widget, *ip).length() < QLineF(caret_position_widget, *closest_ip).length())
      {
        found_intersection = true;
        *closest_ip = *ip;
      }

      painter.save();
      QPen qp(Qt::DashLine);
      qp.setColor(connection_line_color_);
      painter.setPen(qp);
      if (!found_intersection) // no intersection with bounding box of text -> normal drawing
      {
        painter.drawLine(caret_position_widget, position_widget);
        painter.drawLine(caret_position_widget, position_widget);
      }
      else
      {
        painter.drawLine(caret_position_widget, *closest_ip);
        painter.drawLine(caret_position_widget, *closest_ip);
      }
      painter.restore();
      delete(ip);
      delete(closest_ip);
    }

    //painter.drawText(bounding_box_, Qt::AlignLeft, text_);
    //std::cout << "Text to draw: " << st_.text() << " @ " << bounding_box_.topLeft().x() << "," << bounding_box_.topLeft().y() << "\n\n";
    painter.drawStaticText(bounding_box_.topLeft(), st_);
    
    
    //QRect rect = QRect(10, 30, 180, 20);
    //painter.translate( rect.topLeft() );
    //doc_.drawContents( &painter, bounding_box_ );
    //painter.
    
    
    if (selected_)
    {
      drawBoundingBox_(painter);
    }

    painter.restore();
  }

  void Annotation1DCaret::move(const PointType& delta)
  {
    position_.setX(position_.getX() + delta.getX());
    position_.setY(position_.getY() + delta.getY());
  }

  void Annotation1DCaret::setPosition(const Annotation1DCaret::PointType& position)
  {
    position_ = position;
  }

  const Annotation1DCaret::PointType& Annotation1DCaret::getPosition() const
  {
    return position_;
  }

  const Annotation1DCaret::PositionsType& Annotation1DCaret::getCaretPositions() const
  {
    return caret_positions_;
  }

  void Annotation1DCaret::ensureWithinDataRange(Spectrum1DCanvas* const canvas)
  {
    DRange<3> data_range = canvas->getDataRange();

    CoordinateType x_pos = position_.getX();
    CoordinateType y_pos = position_.getY() * canvas->getPercentageFactor();

    if (x_pos < data_range.minPosition()[0])
    {
      position_.setX(data_range.minPosition()[0]);
    }
    if (x_pos > data_range.maxPosition()[0])
    {
      position_.setX(data_range.maxPosition()[0]);
    }
    if (y_pos < data_range.minPosition()[1])
    {
      position_.setY(data_range.minPosition()[1] / canvas->getPercentageFactor());
    }
    if (y_pos > data_range.maxPosition()[1])
    {
      position_.setY(data_range.maxPosition()[1] / canvas->getPercentageFactor());
    }
  }

  void Annotation1DCaret::setColor(const QColor& color)
  {
    color_ = color;
  }

  const QColor& Annotation1DCaret::getColor() const
  {
    return color_;
  }

} // Namespace
