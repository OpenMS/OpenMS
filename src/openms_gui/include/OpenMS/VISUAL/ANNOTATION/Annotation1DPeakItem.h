// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

      // pre-compute bounding box of text_item
      const auto prebox = QApplication::fontMetrics().boundingRect(position_widget.x(), position_widget.y(), 0, 0, Qt::AlignCenter, getText());
      // Shift position of the widget/text, so it sits 'on top' of the peak
      // We can only do that there, since we do not know the state of 'flipped' in general
      // Compute the delta in data-units, NOT pixels, since the shift (up/down, or even left/right) depends on state of 'flipped' and axis 
      const auto deltaXY_in_units = canvas->widgetToDataDistance(prebox.width(), prebox.height()).abs(); // abs() to make sure y axis is not negative
      const auto delta_gravity_in_units = canvas->getGravitator().swap().gravitateZero(deltaXY_in_units); // only keep gravity dim
      // recompute 'position_widget', shifting the text up by 1/2 box
      canvas->dataToWidget(canvas->getMapper().map(position_) + delta_gravity_in_units / 2, position_widget, flipped);
      // re-compute bounding box of text_item on with new position!
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
      if (!text.contains("<\\")) // don't process HTML strings again
      {
        // extract ion index
        {
          QRegularExpression reg_exp(R"(([abcdwxyz])(\d+))");
          QRegularExpressionMatch match = reg_exp.match(text);
          if (text.indexOf(reg_exp) == 0) // only process if at the beginning of the string
          {
            text.replace(reg_exp, "\\1<sub>\\2</sub>");
          }
          else // protein-protein XL specific ion names
          { // e.g. "[alpha|ci$y1]"
            QRegularExpression reg_exp_xlms(R"((ci|xi)[$][abcxyz](\d+))");
            auto match_pos = text.indexOf(reg_exp_xlms);
            if ((match_pos == 6) || (match_pos == 7))
            {
              // set the match_pos to the position of the ion index
              match_pos += 3; // skip "ci$" or "xi$"
              ++match_pos; // skip the ion type (=captured(1))
              QString charge_str = match.captured(2);
              // put sub html tag around number
              text = text.left(match_pos) + QString("<sub>") + charge_str + QString("</sub>") + text.right(text.size() - match_pos - charge_str.size());
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
        QRegularExpression charge_rx(R"([\+|\-](\d+)$)");
        int match_pos = text.indexOf(charge_rx);
        if (match_pos > 0)
        {
          text = text.left(match_pos) + QString("<sup>") + text[match_pos] // + or -
                 + charge_rx.match(text).captured(1) + QString("</sup>");  // charge
        }

        // charge format: z+
        charge_rx = QRegularExpression(R"((\d+)[\+|\-]$)");
        match_pos = text.indexOf(charge_rx);
        if (match_pos > 0)
        {
          auto charge_match = charge_rx.match(text).captured(1);
          text = text.left(match_pos) + QString("<sup>") + charge_match       // charge
                 + text[match_pos + charge_match.size()] + QString("</sup>"); // + or -
        }

        text.replace(QRegularExpression(R"(\+\+$)"), "<sup>2+</sup>");
        text.replace(QRegularExpression(R"(\+$)"), "");
        text.replace(QRegularExpression(R"(\-\-$)"), "<sup>2-</sup>");
        text.replace(QRegularExpression(R"(\-$)"), "");
      }

      text = "<font color=\"" + color_.name() + "\">" + text + "</font>";

      // draw html text
      {
        QTextDocument td;
        td.setHtml(text);
        painter.save();
        double w = td.size().width();
        double h = td.size().height();
        painter.translate(position_widget.x() - w / 2, position_widget.y() - h / 2);
        td.drawContents(&painter);
        painter.restore();
      }
      
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

      // check for newlines in the label and only continue with the first line for charge determination
      peak_anno.remove('\r');
      QStringList lines = peak_anno.split('\n', Qt::SkipEmptyParts);
      if (lines.size() > 1)
      {
        peak_anno = lines[0];
      }

      // regular expression for a charge at the end of the annotation
      QRegularExpression reg_exp(R"(([\+|\-]\d+)$)");

      // read charge and text from annotation item string
      // we support two notations for the charge suffix: '+2' or '++'
      // cut and convert the trailing + or - to a proper charge
      int match_pos = peak_anno.indexOf(reg_exp);
      int tmp_charge(0);
      if (match_pos >= 0)
      {
        tmp_charge = reg_exp.match(peak_anno).captured(1).toInt();
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
