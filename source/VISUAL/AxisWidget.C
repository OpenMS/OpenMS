// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// Qt
#include <QtGui/QPaintEvent>
#include <QtGui/QPainter>
#include <QtGui/QBrush>

// STL
#include <iostream>
#include <algorithm>

// OpenMS
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/VISUAL/AxisPainter.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace std;

namespace OpenMS
{
  using namespace Math;

  AxisWidget::AxisWidget(AxisPainter::Alignment alignment, const char * legend, QWidget * parent) :
    QWidget(parent),
    is_log_(false),
    show_legend_(true),
    alignment_(alignment),
    is_inverse_orientation_(false),
    margin_(0),
    legend_(legend),
    tick_level_(2),
    allow_short_numbers_(false)
  {
    setAxisBounds(0.0, 100.0);

    if (alignment == AxisPainter::RIGHT || alignment == AxisPainter::LEFT)
    {
      setMinimumSize(30, 100);
      setSizePolicy(QSizePolicy::Fixed, QSizePolicy::MinimumExpanding);
    }
    else
    {
      setMinimumSize(100, 30);
      setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
    }

    resize(minimumSize());
  }

  AxisWidget::~AxisWidget()
  {
  }

  void AxisWidget::paintEvent(QPaintEvent * e)
  {
    QPainter painter(this);
    paint(&painter, e);
    painter.end();
  }

  void AxisWidget::paint(QPainter * painter, QPaintEvent * e)
  {
    AxisPainter::paint(painter, e, min_, max_, grid_line_,
                       width(), height(), alignment_, margin_,
                       show_legend_, legend_, allow_short_numbers_, is_log_, is_inverse_orientation_);
  }

  void AxisWidget::setAxisBounds(DoubleReal min, DoubleReal max)
  {
    if (min >= max)
    {
      return;
    }

    if (is_log_)
    {
      //abort if no change
      if (min_ == linear2log(min) && max_ == linear2log(max))
        return;

      min_ = linear2log(min);
      max_ = linear2log(max);

      AxisTickCalculator::calcLogGridLines(min_, max_, grid_line_);
    }
    else
    {
      //abort if no change
      if (min_ == min && max_ == max)
        return;

      min_ = min;
      max_ = max;
      AxisTickCalculator::calcGridLines(min_, max_, grid_line_);
    }

    update();
  }

  void AxisWidget::setLogScale(bool is_log)
  {
    if (is_log_ != is_log)
    {
      is_log_ = is_log;

      if (is_log_)
      {
        setAxisBounds(min_, max_);
      }
      else
      {
        setAxisBounds(log2linear(min_), log2linear(max_));
      }
      update();
    }
  }

  bool AxisWidget::isLogScale()
  {
    return is_log_;
  }

  UInt AxisWidget::margin()
  {
    return margin_;
  }

  void AxisWidget::setMargin(UInt margin)
  {
    if (margin_ != margin)
    {
      margin_ = margin;
      update();
    }
  }

  void AxisWidget::showLegend(bool show_legend)
  {
    if (show_legend_ != show_legend)
    {
      show_legend_ = show_legend;
      if (show_legend_)
      {
        setToolTip("");
      }
      else
      {
        setToolTip(legend_.c_str());
      }
      update();
    }
  }

  bool AxisWidget::isLegendShown() const
  {
    return show_legend_;
  }

  const String & AxisWidget::getLegend()
  {
    return legend_;
  }

  void AxisWidget::setInverseOrientation(bool inverse_orientation)
  {
    if (is_inverse_orientation_ != inverse_orientation)
    {
      is_inverse_orientation_ = inverse_orientation;
      update();
    }
  }

  bool AxisWidget::hasInverseOrientation()
  {
    return is_inverse_orientation_;
  }

  void AxisWidget::setLegend(const String & legend)
  {
    legend_ = legend;
    if (!show_legend_)
    {
      setToolTip(legend_.c_str());
    }
  }

  void AxisWidget::setTickLevel(UInt level)
  {
    if (level == 1 || level == 2)
    {
      tick_level_ = level;
    }
  }

  const AxisWidget::GridVector & AxisWidget::gridLines()
  {
    return grid_line_;
  }

  void AxisWidget::setAllowShortNumbers(bool short_nums)
  {
    allow_short_numbers_ = short_nums;
  }

  DoubleReal AxisWidget::getAxisMinimum() const
  {
    return min_;
  }

  DoubleReal AxisWidget::getAxisMaximum() const
  {
    return max_;
  }

} //Namespace
