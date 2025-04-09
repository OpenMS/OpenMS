// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/AxisPainter.h>
#include <OpenMS/MATH/MathFunctions.h>

using namespace std;

namespace OpenMS
{
  using namespace Math;

  void AxisPainter::paint(QPainter * painter, QPaintEvent *, const double & min, const double & max, const GridVector & grid,
                          const Int width, const Int height, const AxisPainter::Alignment alignment, const UInt margin,
                          bool show_legend, const String& legend, bool shorten_number,
                          bool is_log, bool is_inverse_orientation)
  {
    // position of the widget
    bool horizontal_alignment = (alignment == BOTTOM || alignment == TOP);

    // spacing between ticks and text
    const UInt tick_spacing = 4;

    // determine paintable area without margin
    UInt h = height;
    UInt w = width;
    w = (horizontal_alignment) ? w - margin : w;
    h = (horizontal_alignment) ? h : h - margin;

    // paint axis lines
    switch (alignment)
    {
    case BOTTOM:
      painter->drawLine(0, 0, w - 1, 0);
      break;

    case TOP:
      painter->drawLine(0, h - 1, w - 1, h - 1);
      break;

    case LEFT:
      painter->drawLine(w - 1, 0, w - 1, h - 1);
      break;

    case RIGHT:
      painter->drawLine(0, 0, 0, h - 1);
      break;
    }

    // shrink font size if text does not fit
    UInt font_size = painter->font().pointSize();
    UInt max_width = 0;

    if (!grid.empty()) //check big intervals only
    {
      QFontMetrics metrics(QFont(painter->font().family(), font_size));
      for (Size i = 0; i < grid[0].size(); i++)
      {
        QString tmp;
        if (shorten_number)
        {
          getShortenedNumber_(tmp, scale_(grid[0][i], is_log));
        }
        else
        {
          tmp = QString("%1").arg(scale_(grid[0][i], is_log));
        }
        QRect rect = metrics.boundingRect(tmp);
        max_width = std::max(max_width, (UInt)rect.width());
      }
    }

    //shrink font if too much text it displayed
    UInt overall_required_pixels = 0;
    if (horizontal_alignment)
    {
      Size tick_count = 0;
      for (Size i = 0; i < grid.size(); i++)
      {
        tick_count += grid[i].size();
      }
      overall_required_pixels = (UInt)(max_width * tick_count);
    }
    else // Shrink font if the largest text is too big
    {
      overall_required_pixels = UInt(max_width  + 0.25 * font_size + tick_spacing);
    }

    if (w < overall_required_pixels)
    {
      font_size = UInt(font_size * w / overall_required_pixels);
    }

    // Painting tick levels
    for (Size i = 0; i != grid.size(); i++)
    {
      // Just draw text on big intervals
      if (is_log && i > 0)
      {
        break;
      }

      QColor text_color;
      UInt tick_size = 0;
      QFontMetrics metrics(painter->font());

      if (i == 0) // big intervals
      {
        painter->setFont(QFont(painter->font().family(), UInt(font_size)));
        metrics = QFontMetrics(painter->font());
        tick_size = UInt(0.33 * font_size);
        text_color = QColor(0, 0, 0);
      }
      else // small intervals
      {
        painter->setFont(QFont(painter->font().family(), UInt(0.8 * font_size)));
        metrics = QFontMetrics(painter->font());
        tick_size = UInt(0.25 * font_size);
        text_color = QColor(20, 20, 20);
      }

      // painting all ticks of the level
      UInt i_beg = (horizontal_alignment) ? 0 : h;
      UInt i_end = (horizontal_alignment) ? w : 0;
      for (Size j = 0; j != grid[i].size(); j++)
      {
        UInt tick_pos;
        if (is_inverse_orientation)
        {
          tick_pos = UInt(intervalTransformation(grid[i][j], min, max, i_end, i_beg)) + ((alignment == LEFT || alignment == RIGHT) ? -1 : 1) * margin;
        }
        else
        {
          tick_pos = UInt(intervalTransformation(grid[i][j], min, max, i_beg, i_end));
        }

        // paint ticks
        painter->setPen(QPen(Qt::black));
        switch (alignment)
        {
        case BOTTOM:
          painter->drawLine(tick_pos, 0, tick_pos, tick_size);
          break;

        case TOP:
          painter->drawLine(tick_pos, h, tick_pos, h - tick_size);
          break;

        case LEFT:
          painter->drawLine(w - tick_size, tick_pos + margin, w, tick_pos + margin);
          break;

        case RIGHT:
          painter->drawLine(0, tick_pos + margin, tick_size, tick_pos + margin);
          break;
        }

        // values at axis lines

        // count number of grid lines
        Size gl_count = 0;
        for (GridVector::const_iterator it = grid.begin(); it != grid.end(); ++it)
        {
          gl_count += it->size();
        }
        //  if there are too many grid lines, hide every odd minor grid line value
        if (gl_count >= 30 && i == 1)
        {
          // needed because skips occur in small grid lines at the position of big grid lines
          double dist_small = std::min<double>(fabs(grid[1][1] - grid[1][0]), fabs(grid[1][2] - grid[1][1]));
          UInt n = Math::round((grid[1][j] - grid[0][0]) / dist_small);
          if (n % 2 == 1)
          {
            continue;
          }
        }

        QString text;
        if (shorten_number)
        {
          getShortenedNumber_(text, scale_(grid[i][j], is_log));
        }
        else
        {
          text = QString("%1").arg(scale_(grid[i][j], is_log));
        }

        painter->setPen(QPen(text_color));

        // get bounding rectangle for text we want to layout
        QRect textbound = metrics.boundingRect(text);

        // Calculate text position
        UInt x_pos = 0;
        switch (alignment)
        {
        case BOTTOM:
        case TOP:
          x_pos = tick_pos - UInt(0.5 * textbound.width());
          break;

        case LEFT:
          x_pos = w - (tick_size + tick_spacing) - textbound.width();
          break;

        case RIGHT:
          x_pos = tick_size + tick_spacing;
          break;
        }

        UInt y_pos = 0;
        switch (alignment)
        {
        case BOTTOM:
          y_pos = tick_size + tick_spacing + UInt(0.5 * textbound.height());
          break;

        case TOP:
          y_pos = h - (tick_size + tick_spacing);
          break;

        case LEFT:
        case RIGHT:
          y_pos = tick_pos + margin + UInt(0.25 * textbound.height());
          break;
        }
        painter->drawText(x_pos, y_pos, text);
      }
    }

    // Painting legend
    painter->setFont(QFont(painter->font().family(), font_size));
    if (show_legend && !legend.empty())
    {
      // style settings
      painter->setPen(QPen(Qt::black));

      switch (alignment)
      {
      case BOTTOM:
        painter->drawText(0, 0, w, h, Qt::AlignBottom | Qt::AlignHCenter, legend.c_str());
        break;

      case TOP:
        painter->drawText(0, 0, w, h, Qt::AlignTop | Qt::AlignHCenter, legend.c_str());
        break;

      case LEFT:
        painter->rotate(270);
        painter->drawText(-(int)h, 0, h, w, Qt::AlignHCenter | Qt::AlignTop, legend.c_str());
        break;

      case RIGHT:
        painter->rotate(270);
        painter->drawText(-(int)h, 0, h, w, Qt::AlignHCenter | Qt::AlignBottom, legend.c_str());
        break;
      }
    }
  }

  void AxisPainter::getShortenedNumber_(QString& short_num, double number)
  {
    if (number < 1e3)
    {
      short_num = QString("%1").arg(number);
    }
    else if (number < 1e6)
    {
      short_num = QString("%1k").arg(Math::roundDecimal(number /1e3, -2));
    }
    else if (number < 1e9)
    {
      short_num = QString("%1M").arg(number / 1e6);
    }
    else
    {
      short_num = QString("%1G").arg(number / 1e9);
    }
  }

  double AxisPainter::scale_(double x, bool is_log)
  { // 
    return (is_log) ? Math::roundDecimal(pow(10, x), -8) : Math::roundDecimal(x, -8);
  }

}
