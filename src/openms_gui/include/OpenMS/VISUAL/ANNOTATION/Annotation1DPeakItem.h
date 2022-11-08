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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

#include <QtGui/QColor>

namespace OpenMS
{

  /** @brief A peak annotation item
            @see Annotation1DItem
    */
  template <class DataPoint> // e.g.  Peak1D
  class Annotation1DPeakItem :
    public Annotation1DItem
  {
public:
    /// Constructor
    Annotation1DPeakItem(const DataPoint& peak_position, const QString& text, const QColor& color) :
      Annotation1DItem(text), peak_position_(peak_position), position_(peak_position), color_(color)
    {
    }

    /// Copy constructor
    Annotation1DPeakItem(const Annotation1DPeakItem& rhs) = default;

    /// Destructor
    ~Annotation1DPeakItem() override = default;

    // Docu in base class
    void draw(Plot1DCanvas* const canvas, QPainter& painter, bool flipped = false) override
    {
      painter.save();

      painter.setPen(color_);

      QPoint position_widget, peak_position_widget;

      // translate units to pixel coordinates
      canvas->dataToWidget(canvas->getMapper().map(position_), position_widget, flipped);
      canvas->dataToWidget(canvas->getMapper().map(peak_position_), peak_position_widget, flipped);

      // compute bounding box of text_item on the specified painter
      bounding_box_ = QApplication::fontMetrics().boundingRect(position_widget.x(), position_widget.y(), 0, 0, Qt::AlignCenter, getText());

      // draw connection line between anchor point and current position if pixel coordinates differ significantly
      if ((position_widget - peak_position_widget).manhattanLength() > 2)
      {
        QPointF border_point = GUIHelpers::intersectionPoint(bounding_box_, peak_position_widget);
        if (bounding_box_.center() != border_point)
        {
          painter.save();
          painter.setPen(Qt::DashLine);
          painter.drawLine(peak_position_widget, border_point);
          painter.restore();
        }
      }

      // some pretty printing
      QString text = text_;
      if (!text.contains(R"(<\)")) // don't process HTML strings again
      {
        // extract ion index
        {
          QRegExp reg_exp(R"([abcdwxyz](\d+))");
          int match_pos = reg_exp.indexIn(text);

          if (match_pos == 0)
          {
            QString index_str = reg_exp.cap(1);

            // put sub html tag around number
            text = text[match_pos] + QString("<sub>") + index_str + QString("</sub>") + text.right(text.size() - match_pos - index_str.size() - 1);
          }
          else // protein-protein XL specific ion names
          {
            QRegExp reg_exp_xlms(R"((ci|xi)[$][abcxyz](\d+))");
            match_pos = reg_exp_xlms.indexIn(text);
            if ((match_pos == 6) || (match_pos == 7))
            {
              // set the match_pos to the position of the ion index
              match_pos += 3;
              QString index_str = reg_exp.cap(1);

              // put sub html tag around number
              text = text.left(match_pos) + text[match_pos] + QString("<sub>") + index_str + QString("</sub>") + text.right(text.size() - match_pos - index_str.size() - 1);
            }
          }
        }

        // common losses
        text.replace("H2O1", "H<sub>2</sub>O"); // mind the order with H2O substitution
        text.replace("H2O", "H<sub>2</sub>O");
        text.replace("NH3", "NH<sub>3</sub>");
        text.replace("H3N1", "NH<sub>3</sub>");
        text.replace("C1H4O1S1", "H<sub>4</sub>COS"); // methionine sulfoxide loss

        // nucleotide XL related losses
        text.replace("H3PO4", "H<sub>3</sub>PO<sub>4</sub>");
        text.replace("HPO3", "HPO<sub>3</sub>");
        text.replace("C3O", "C<sub>3</sub>O");

        // charge format: +z
        QRegExp charge_rx(R"([\+|\-](\d+)$)");
        int match_pos = charge_rx.indexIn(text);
        if (match_pos > 0)
        {
          text = text.left(match_pos) + QString("<sup>") + text[match_pos] // + or -
                 + charge_rx.cap(1) + QString("</sup>");                   // charge
        }

        // charge format: z+
        charge_rx = QRegExp(R"((\d+)[\+|\-]$)");
        match_pos = charge_rx.indexIn(text);
        if (match_pos > 0)
        {
          text = text.left(match_pos) + QString("<sup>") + charge_rx.cap(1)       // charge
                 + text[match_pos + charge_rx.cap(1).size()] + QString("</sup>"); // + or -
        }

        text.replace(QRegExp(R"(\+\+$)"), "<sup>2+</sup>");
        text.replace(QRegExp(R"(\+$)"), "");
        text.replace(QRegExp(R"(\-\-$)"), "<sup>2-</sup>");
        text.replace(QRegExp(R"(\-$)"), "");
      }

      text = "<font color=\"" + color_.name() + "\">" + text + "</font>";
      QTextDocument td;
      td.setHtml(text);

      // draw html text
      painter.save();
      double w = td.size().width();
      double h = td.size().height();
      painter.translate(position_widget.x() - w / 2, position_widget.y() - h / 2);
      td.drawContents(&painter);
      painter.restore();

      if (selected_)
      {
        drawBoundingBox_(painter);
      }

      painter.restore();
    }

    // Docu in base class
    void move(const PointXYType delta, const Gravitator& /*gr*/, const DimMapper<2>& dim_mapper) override
    {
      auto pos_xy = dim_mapper.map(position_);
      pos_xy += delta;
      dim_mapper.fromXY(pos_xy, position_);
    }

    /// Sets the position of the label
    void setPosition(const DataPoint& position)
    {
      position_ = position;
    }

    /// Returns the position of the label (peak)
    const DataPoint& getPosition() const
    {
      return position_;
    }

    /// Returns the position of the annotated peak
    const DataPoint& getPeakPosition() const
    {
      return peak_position_;
    }

    // Docu in base class
    void ensureWithinDataRange(Plot1DCanvas* const canvas, const int layer_index) override
    {
      canvas->pushIntoDataRange(position_, layer_index);
    }

    /// Set the color of the label
    void setColor(const QColor& color)
    {
      color_ = color;
    }

    /// Returns the color of the label
    const QColor& getColor() const
    {
      return color_;
    }

    /// Convert the 'text()' to a Peptide::PeakAnnotation
    PeptideHit::PeakAnnotation toPeakAnnotation() const
    {
      // add new fragment annotation
      QString peak_anno = this->getText().trimmed();

      // regular expression for a charge at the end of the annotation
      QRegExp reg_exp(R"(([\+|\-]\d+)$)");

      // check for newlines in the label and only continue with the first line for charge determination
      QStringList lines = peak_anno.split(QRegExp("[\r\n]"), QString::SkipEmptyParts);
      if (lines.size() > 1)
      {
        peak_anno = lines[0];
      }

      // read charge and text from annotation item string
      // we support two notations for the charge suffix: '+2' or '++'
      // cut and convert the trailing + or - to a proper charge
      int match_pos = reg_exp.indexIn(peak_anno);
      int tmp_charge(0);
      if (match_pos >= 0)
      {
        tmp_charge = reg_exp.cap(1).toInt();
        peak_anno = peak_anno.left(match_pos);
      }
      else
      {
        // count number of + and - in suffix (e.g., to support "++" as charge 2 annotation)
        int plus(0), minus(0);

        for (int p = (int)peak_anno.size() - 1; p >= 0; --p)
        {
          if (peak_anno[p] == '+')
          {
            ++plus;
            continue;
          }
          else if (peak_anno[p] == '-')
          {
            ++minus;
            continue;
          }
          else // not '+' or '-'?
          {
            if (plus > 0 && minus == 0) // found pluses?
            {
              tmp_charge = plus;
              peak_anno = peak_anno.left(peak_anno.size() - plus);
              break;
            }
            else if (minus > 0 && plus == 0) // found minuses?
            {
              tmp_charge = -minus;
              peak_anno = peak_anno.left(peak_anno.size() - minus);
              break;
            }
            break;
          }
        }
      }

      PeptideHit::PeakAnnotation fa;
      fa.charge = tmp_charge;
      fa.mz = this->getPeakPosition().getMZ();
      fa.intensity = this->getPeakPosition().getIntensity();
      if (lines.size() > 1)
      {
        peak_anno.append("\n").append(lines[1]);
      }
      fa.annotation = peak_anno;

      return fa;
    }

    // Docu in base class
    Annotation1DItem* clone() const override
    {
      return new Annotation1DPeakItem(*this);
    }

  protected:
    /// The position of the anchor (e.g. the Peak1D)
    DataPoint peak_position_;

    /// The position of the label (e.g. the Peak1D)
    DataPoint position_;

    /// The color of the label
    QColor color_;
  };
} // namespace OpenMS
