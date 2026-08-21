// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>

#include <OpenMS/CHEMISTRY/IonNaming.h>

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

#include <QtGui/QColor>
#include <QtGui/QFontMetricsF>

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
    /**
      @brief Constructor

      @p charge is the charge of the annotated peak for producers that keep it in
      PeptideHit::PeakAnnotation::charge rather than spelling it out in @p text. It is drawn next to the
      ion name but deliberately kept out of @p text, so that reading the item back yields the name the
      producer wrote: putting it into the text made every redraw rewrite the stored annotation.
      Pass 0 (the default) when the text is all there is, as for a plain user label.

      @param[in] peak_position The peak this annotation belongs to
      @param[in] text The annotation text, i.e. the ion name as the producer wrote it
      @param[in] color The colour to draw it in
      @param[in] charge Charge of the peak when @p text does not spell one out, 0 otherwise
    */
    Annotation1DPeakItem(const DataPoint& peak_position, const QString& text, const QColor& color, int charge = 0) :
      Annotation1DItem(text), peak_position_(peak_position), position_(peak_position), color_(color), charge_(charge)
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
      const QFontMetricsF fm(QApplication::font());
      const auto prebox = fm.boundingRect(QRectF(position_widget.x(), position_widget.y(), 0, 0), Qt::AlignCenter, getText());
      // Shift position of the widget/text, so it sits 'on top' of the peak
      // We can only do that there, since we do not know the state of 'flipped' in general
      // Compute the delta in data-units, NOT pixels, since the shift (up/down, or even left/right) depends on state of 'flipped' and axis 
      const auto deltaXY_in_units = canvas->widgetToDataDistance(prebox.width(), prebox.height()).abs(); // abs() to make sure y axis is not negative
      const auto delta_gravity_in_units = canvas->getGravitator().swap().gravitateZero(deltaXY_in_units); // only keep gravity dim
      // recompute 'position_widget', shifting the text up by 1/2 box
      canvas->dataToWidget(canvas->getMapper().map(position_) + delta_gravity_in_units / 2, position_widget, flipped);
      // re-compute bounding box of text_item on with new position!
      bounding_box_ = fm.boundingRect(QRectF(position_widget.x(), position_widget.y(), 0, 0), Qt::AlignCenter, getText());


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

        // A charge that is not in the name at all (see the constructor) is drawn from the member, so
        // that it never has to travel through the label and get written back into the stored annotation.
        // Charge 1 is implied and left off, exactly as it is for a named ion like "y3+" whose lone sign
        // is stripped below; only |charge| >= 2 is drawn, and next to the ion name (end of the first
        // line) rather than after a trailing free-text comment.
        if (charge_ != 0 && IonNaming::chargeFromName(fromQString(text_)) == 0)
        {
          // magnitude in a wider type: std::abs(INT_MIN) is undefined behaviour
          const long long magnitude = (charge_ < 0) ? -static_cast<long long>(charge_) : static_cast<long long>(charge_);
          if (magnitude > 1)
          {
            const QString sup = QString("<sup>") + QString::number(magnitude) + ((charge_ < 0) ? "-" : "+") + QString("</sup>");
            static const QRegularExpression first_break(R"([\r\n])");
            const int eol = text.indexOf(first_break);
            text.insert((eol < 0) ? text.size() : eol, sup);
          }
        }

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

      // check for newlines in the label and only continue with the first line for charge determination.
      // KeepEmptyParts: a blank line inside a label is the user's text, and dropping it here would
      // silently rewrite the label every time the annotations are read back.
      // Same rule as the identification view uses when it builds the label: "\r\n" first in the
      // alternation so a CRLF is one break, a lone CR is a break of its own, and KeepEmptyParts so a
      // blank line the user typed survives. IonNaming::chargeFromName likewise ends the ion name at
      // either character, so all three agree on where the first line stops.
      static const QRegularExpression line_break(R"(\r\n|[\r\n])");
      QStringList lines = peak_anno.split(line_break, Qt::KeepEmptyParts);
      if (lines.size() > 1)
      {
        peak_anno = lines[0];
      }

      // Read the charge from the label but leave the label alone. Cutting the charge off used to
      // corrupt the stored annotation: the identification view puts the charge back on redraw, so a
      // round trip through this function turned "y3+" into "y3" and then into "y3++".
      // IonNaming::chargeFromName understands all notations that occur ('+2', '++' and mzPAF's '^2').
      // A charge the user typed into the label wins; otherwise fall back to the one the item was built
      // with. The name itself is returned untouched, so viewing a spectrum cannot rewrite it.
      const int named_charge = IonNaming::chargeFromName(fromQString(peak_anno));
      const int tmp_charge = (named_charge != 0) ? named_charge : charge_;

      PeptideHit::PeakAnnotation fa;
      fa.charge = tmp_charge;
      fa.mz = this->getPeakPosition().getMZ();
      fa.intensity = this->getPeakPosition().getIntensity();
      for (int l = 1; l < lines.size(); ++l)
      { // keep every extra label line, not just the first (they are free text below the ion name)
        peak_anno.append("\n").append(lines[l]);
      }
      fa.annotation = fromQString(peak_anno);

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

    /// Charge of the annotated peak when its name does not spell one out; 0 if unknown or already named
    int charge_ = 0;
  };
} // namespace OpenMS
