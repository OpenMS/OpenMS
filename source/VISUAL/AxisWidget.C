// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
