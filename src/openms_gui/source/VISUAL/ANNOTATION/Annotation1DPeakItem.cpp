// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

#include <QtGui/QPainter>
#include <QtCore/QPoint>

namespace OpenMS
{

  Annotation1DPeakItem::Annotation1DPeakItem(const PointType & peak_position, const QString & text, const QColor & color) :
    Annotation1DItem(text),
    peak_position_(peak_position),
    position_(peak_position),
    color_(color)
  {
  }

  Annotation1DPeakItem::Annotation1DPeakItem(const Annotation1DPeakItem & rhs) :
    Annotation1DItem(rhs)
  {
    peak_position_ = rhs.peak_position_;
    position_ = rhs.position_;
    color_ = rhs.color_;
  }

  Annotation1DPeakItem::~Annotation1DPeakItem()
  {
  }

  QRectF Annotation1DPeakItem::calculateBoundingBox(
    const PointType & peak_position,
    const PointType & position,
    const QString & text,
    Spectrum1DCanvas * const canvas,
    bool flipped,
    QPoint & position_widget,
    QPoint & peak_position_widget,
    double & horizontal_shift,
    double & vertical_shift)
  {
    // translate mz/intensity to pixel coordinates
    canvas->dataToWidget(position.getX(), position.getY(), position_widget, flipped, true);
    canvas->dataToWidget(peak_position.getX(), peak_position.getY(), peak_position_widget, flipped, true);

    // compute bounding box of text_item on the specified painter
    QRectF bounding_box = QApplication::fontMetrics().boundingRect(
      position_widget.x(),
      position_widget.y(),
      0, 0,
      Qt::AlignCenter,
      text);

    vertical_shift = 0;
    horizontal_shift = 0;

    if (canvas->isMzToXAxis())
    {
      // shift pos - annotation should be over peak or, if not possible, next to it
      vertical_shift = bounding_box.height() / 2 + 5;
      if (!flipped)
      {
        vertical_shift *= -1;
      }

      bounding_box.translate(0.0, vertical_shift);

      if (flipped && bounding_box.bottom() > canvas->height())
      {
        bounding_box.moveBottom(canvas->height());
        bounding_box.moveLeft(position_widget.x() + 5.0);
      }
      else if (!flipped && bounding_box.top() < 0.0)
      {
        bounding_box.moveTop(0.0);
        bounding_box.moveLeft(position_widget.x() + 5.0);
      }
    }
    else
    {
      // annotation should be next to the peak (to its right)
      horizontal_shift = bounding_box.width() / 2 + 5;
      bounding_box.translate(horizontal_shift, 0.0);
      if (bounding_box.right() > canvas->width())
      {
        bounding_box.moveRight(canvas->width());
      }
    }

    return bounding_box;
  }

  void Annotation1DPeakItem::draw(Spectrum1DCanvas * const canvas, QPainter & painter, bool flipped)
  {
    painter.save();

    painter.setPen(color_);

    QPoint position_widget, peak_position_widget;
    double horizontal_shift(0), vertical_shift(0);

    bounding_box_ = calculateBoundingBox(peak_position_,
      position_,
      text_,
      canvas,
      flipped,
      position_widget,
      peak_position_widget,
      horizontal_shift,
      vertical_shift);

    // draw connection line between anchor point and current position if pixel coordinates differ significantly
    if ((position_widget - peak_position_widget).manhattanLength() > 2)
    {
      // check if line crosses bounding box, if so move line startpoint to correct bounding box intersection
      QLineF line(peak_position_widget, position_widget + QPoint(horizontal_shift, vertical_shift));
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
          QLineF(peak_position_widget, *ip).length() < QLineF(peak_position_widget, *closest_ip).length())
      {
        found_intersection = true;
        *closest_ip = *ip;
      }
      // intersection with left
      itype = line.intersect(left, ip);
      if (itype == QLineF::BoundedIntersection &&
          QLineF(peak_position_widget, *ip).length() < QLineF(peak_position_widget, *closest_ip).length())
      {
        found_intersection = true;
        *closest_ip = *ip;
      }

      // intersection with right
      itype = line.intersect(right, ip);
      if (itype == QLineF::BoundedIntersection &&
          QLineF(peak_position_widget, *ip).length() < QLineF(peak_position_widget, *closest_ip).length())
      {
        found_intersection = true;
        *closest_ip = *ip;
      }

      // intersection with bottom
      itype = line.intersect(bottom, ip);
      if (itype == QLineF::BoundedIntersection &&
          QLineF(peak_position_widget, *ip).length() < QLineF(peak_position_widget, *closest_ip).length())
      {
        found_intersection = true;
        *closest_ip = *ip;
      }

      painter.save();
      painter.setPen(Qt::DashLine);
      if (!found_intersection) // no intersection with bounding box of text -> normal drawing
      {
        painter.drawLine(peak_position_widget, peak_position_widget + QPoint(0, vertical_shift));
        painter.drawLine(peak_position_widget + QPoint(0, vertical_shift), position_widget);
      }
      else
      {
        QPoint vertical_end = QPoint(peak_position_widget.x(), position_widget.y() + 20);

        if (vertical_end.y() < peak_position_widget.y()) // label is above peak
        {
          painter.drawLine(peak_position_widget, vertical_end);
          painter.drawLine(vertical_end, position_widget);
        }
        else // label is below peak
        {
          painter.drawLine(peak_position_widget, *closest_ip);
          painter.drawLine(peak_position_widget, *closest_ip);
        }
      }
      painter.restore();
      delete(ip);
      delete(closest_ip);
    }

    // some pretty printing
    QString text = text_;
    if (!text.contains("<\\")) // don't process HTML strings again
    {
      // extract ion index
      {
        QRegExp reg_exp("[abcdwxyz](\\d+)");
        int match_pos = reg_exp.indexIn(text);

        if (match_pos == 0)
        {
          QString index_str = reg_exp.cap(1);

          // put sub html tag around number
          text = text[match_pos]
                + QString("<sub>") + index_str + QString("</sub>")
                + text.right(text.size() - match_pos - index_str.size() - 1);
        } 
        else // protein-protein XL specific ion names
        {
          QRegExp reg_exp_xlms("(ci|xi)[$][abcxyz](\\d+)");
          match_pos = reg_exp_xlms.indexIn(text);
          if ( (match_pos == 6) || (match_pos == 7))
          {
            // set the match_pos to the position of the ion index
            match_pos += 3;
            QString index_str = reg_exp.cap(1);

            // put sub html tag around number
            text = text.left(match_pos)
                  + text[match_pos]
                  + QString("<sub>") + index_str + QString("</sub>")
                  + text.right(text.size() - match_pos - index_str.size() - 1);
          }
        }
      }

      // common losses
      text.replace("H2O1","H<sub>2</sub>O"); // mind the order with H2O substitution
      text.replace("H2O","H<sub>2</sub>O");
      text.replace("NH3","NH<sub>3</sub>");
      text.replace("H3N1","NH<sub>3</sub>");
      text.replace("C1H4O1S1", "H<sub>4</sub>COS");  // methionine sulfoxide loss

      // nucleotide XL realted losses
      text.replace("H3PO4","H<sub>3</sub>PO<sub>4</sub>");
      text.replace("HPO3","HPO<sub>3</sub>");
      text.replace("C3O","C<sub>3</sub>O");

      // charge format: +z
      QRegExp charge_rx("[\\+|\\-](\\d+)$");
      int match_pos = charge_rx.indexIn(text);
      if (match_pos > 0)
      {
        text = text.left(match_pos)
               + QString("<sup>") + text[match_pos] // + or - 
               + charge_rx.cap(1) + QString("</sup>"); // charge
      }

      // charge format: z+
      charge_rx = QRegExp("(\\d+)[\\+|\\-]$");
      match_pos = charge_rx.indexIn(text);
      if (match_pos > 0)
      {
        text = text.left(match_pos)
               + QString("<sup>") + charge_rx.cap(1) // charge 
               + text[match_pos + charge_rx.cap(1).size()] + QString("</sup>"); // + or -
      }

      text.replace(QRegExp("\\+\\+$"), "<sup>2+</sup>");
      text.replace(QRegExp("\\+$"), "");
      text.replace(QRegExp("\\-\\-$"), "<sup>2-</sup>");
      text.replace(QRegExp("\\-$"), "");

    }

    text = "<font color=\"" + color_.name() + "\">" + text + "</font>";
    QTextDocument td;
    td.setHtml(text);

    //draw html text
    painter.save();
    double w = td.size().width();
    double h = td.size().height();
    painter.translate(position_widget.x() - w/2, position_widget.y() - h);
    td.drawContents(&painter);
    painter.restore();

    //painter.drawText(bounding_box_, Qt::AlignCenter, text_);
    if (selected_) { drawBoundingBox_(painter); }

    painter.restore();
  }

  void Annotation1DPeakItem::move(const PointType & delta)
  {
    position_.setX(position_.getX() + delta.getX());
    position_.setY(position_.getY() + delta.getY());
  }

  void Annotation1DPeakItem::setPosition(const Annotation1DPeakItem::PointType & position)
  {
    position_ = position;
  }

  const Annotation1DPeakItem::PointType & Annotation1DPeakItem::getPosition() const
  {
    return position_;
  }

  const Annotation1DPeakItem::PointType & Annotation1DPeakItem::getPeakPosition() const
  {
    return peak_position_;
  }

  void Annotation1DPeakItem::ensureWithinDataRange(Spectrum1DCanvas * const canvas)
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

  void Annotation1DPeakItem::setColor(const QColor & color)
  {
    color_ = color;
  }

  const QColor & Annotation1DPeakItem::getColor() const
  {
    return color_;
  }

} // Namespace
