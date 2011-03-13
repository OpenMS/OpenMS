// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	AxisWidget::AxisWidget(Alignment alignment, const char* legend, QWidget* parent)
		: QWidget( parent),
			is_log_(false),
			show_legend_(true),
			alignment_(alignment),
			inverse_orientation_(false),
			margin_(0),
			legend_(legend),
      tick_level_(2),
			allow_short_numbers_(false)
	{
		setAxisBounds(0.0,100.0);
		
	 	if (alignment==RIGHT || alignment==LEFT )	
	 	{
			setMinimumSize(30,100);
			setSizePolicy(QSizePolicy::Fixed,QSizePolicy::MinimumExpanding);
		}
		else
		{
			setMinimumSize(100, 30);
			setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::Fixed);
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

  void AxisWidget::paint(QPainter* painter, QPaintEvent*)
  {
    //position of the widget
    bool horizontal_alignment = (alignment_ == BOTTOM || alignment_ == TOP);

    //spacing between ticks and text
    static UInt tick_spacing = 4;

    // determine paintable area without margin
    UInt h = height();
    UInt w = width();
    w = (horizontal_alignment) ? w-margin_ : w;
    h = (horizontal_alignment) ? h : h-margin_;

    // paint axis lines
    switch (alignment_)
    {
    case BOTTOM:
      painter->drawLine(0, 0, w-1, 0);
      break;
    case TOP:
      painter->drawLine(0, h-1, w-1, h-1);
      break;
    case LEFT:
      painter->drawLine(w-1, 0 , w-1, h-1);
      break;
    case RIGHT:
      painter->drawLine(0, 0 , 0, h-1);
      break;
    }

    // shrink font size if text does not fit
    UInt font_size = font().pointSize();
    UInt max_width = 0;
    if (grid_line_.size()>=1) //check big intervals only
    {
      QFontMetrics metrics(QFont(font().family(),font_size));
      for (Size i=0; i<grid_line_[0].size(); i++)
      {
        QString tmp;
        getShortenedNumber_(tmp, scale_(grid_line_[0][i]));
        QRect rect = metrics.boundingRect(tmp);
        max_width = max(max_width,(UInt)rect.width());
      }
    }
    //shrink font if too much text it displayed
    UInt overall_required_pixels = 0;
    if (horizontal_alignment)
    {
      Size tick_count = 0;
      for (Size i=0; i<grid_line_.size(); i++)
      {
        tick_count += grid_line_[i].size();
      }
      overall_required_pixels = (UInt)(max_width*tick_count);
    }
    else //shrink font if the largest text is too big
    {
      overall_required_pixels = UInt(max_width  + 0.25*font_size + tick_spacing);
    }

    if (w<overall_required_pixels)
    {
      font_size = UInt(font_size * w / overall_required_pixels);
    }

    //painting tick levels
    for (Size i = 0; i!=grid_line_.size(); i++)
    {
      // iust draw text on big intervalls
      if (is_log_ && i>0) break;

      QColor text_color;
      UInt tick_size = 0;
      QFontMetrics metrics(font());
      if (i==0) //big intervals
      {
        painter->setFont(QFont(font().family(), UInt(font_size)));
        metrics = QFontMetrics(painter->font());
        tick_size = UInt(0.33 * font_size);
        text_color = QColor(0, 0, 0);
      }
      else //small intervals
      {
        painter->setFont(QFont(font().family(),UInt(0.8*font_size)));
        metrics = QFontMetrics(painter->font());
        tick_size = UInt(0.25 * font_size);
        text_color = QColor(20, 20, 20);
      }

      //painting all ticks of the level
      UInt i_beg = (horizontal_alignment)? 0 : h;
      UInt i_end = (horizontal_alignment)? w : 0;
      for (Size j = 0; j!=grid_line_[i].size(); j++)
      {
        UInt tick_pos;
        if (inverse_orientation_)
        {
          tick_pos = UInt(intervalTransformation(grid_line_[i][j], min_, max_, i_end, i_beg)) + ((alignment_==LEFT || alignment_==RIGHT)?-1:1)*margin_;
        }
        else
        {
          tick_pos = UInt(intervalTransformation(grid_line_[i][j], min_, max_, i_beg, i_end));
        }

        //paint ticks
        painter->setPen(QPen(Qt::black));
        switch (alignment_)
        {
          case BOTTOM:
            painter->drawLine(tick_pos, 0, tick_pos, tick_size);
            break;
          case TOP:
            painter->drawLine(tick_pos, h, tick_pos,  h-tick_size);
            break;
          case LEFT:
            painter->drawLine(w-tick_size, tick_pos+margin_, w, tick_pos+margin_);
            break;
          case RIGHT:
            painter->drawLine(0, tick_pos+margin_, tick_size, tick_pos+margin_);
            break;
        }

        // values at axis lines

        // count number of grid lines
        Size gl_count = 0;
        for (GridVector::iterator it = grid_line_.begin(); it != grid_line_.end(); ++it)
        {
          gl_count += it->size();
        }
        //  if there are too many grid lines, hide every odd minor grid line value
        if (gl_count >= 30 && i == 1)
        {
          // needed because skips occur in small grid lines at the position of big grid lines
          DoubleReal dist_small = std::min<DoubleReal>(fabs(grid_line_[1][1] - grid_line_[1][0]), fabs(grid_line_[1][2] - grid_line_[1][1]));
          UInt n = Math::round((grid_line_[1][j] - grid_line_[0][0]) / dist_small);
          if (n % 2 == 1)
          {
            continue;
          }
        }

        QString text;
        getShortenedNumber_(text, scale_(grid_line_[i][j]));
        painter->setPen(QPen(text_color));

        // get bounding rectangle for text we want to layout
        QRect textbound = metrics.boundingRect(text);

        // Calculate text position
        UInt x_pos = 0;
        switch (alignment_)
        {
          case BOTTOM:
          case TOP:
            x_pos = tick_pos - UInt(0.5*textbound.width());
            break;
          case LEFT:
            x_pos = w - (tick_size + tick_spacing) - textbound.width();
            break;
          case RIGHT:
            x_pos = tick_size + tick_spacing;
            break;
        }

        UInt y_pos = 0;
        switch (alignment_)
        {
          case BOTTOM:
            y_pos = tick_size + tick_spacing + UInt(0.5*textbound.height());
            break;
          case TOP:
            y_pos = h - (tick_size + tick_spacing);
            break;
          case LEFT:
          case RIGHT:
            y_pos = tick_pos + margin_ + UInt(0.25*textbound.height());
            break;
        }
        painter->drawText(x_pos, y_pos, text);
      }
    }

    //painting legend
    if (show_legend_ && legend_!="")
    {
      // style settings
      painter->setFont(font());
      painter->setPen(QPen(Qt::black));

      switch (alignment_)
      {
        case BOTTOM:
          painter->drawText(0, 0 ,  w, h, Qt::AlignBottom|Qt::AlignHCenter, legend_.c_str());
          break;
        case TOP:
          painter->drawText(0, 0 ,  w, h, Qt::AlignTop|Qt::AlignHCenter, legend_.c_str());
          break;
        case LEFT:
          painter->rotate(270);
          painter->drawText(-(int)h, 0 ,h ,w, Qt::AlignHCenter|Qt::AlignTop, legend_.c_str());
          break;
        case RIGHT:
          painter->rotate(270);
          painter->drawText(-(int)h, 0 ,h ,w, Qt::AlignHCenter|Qt::AlignBottom, legend_.c_str());
          break;
      }
    }

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
			if (min_== linear2log(min) && max_ == linear2log(max)) return;
			
			min_ = linear2log(min);
			max_ = linear2log(max);

			AxisTickCalculator::calcLogGridLines(min_,max_,grid_line_);
		}
		else
		{
			//abort if no change
			if (min_==min && max_==max) return;
			
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
				setAxisBounds(min_,max_);
			}
			else 
			{
				setAxisBounds(log2linear(min_),log2linear(max_));
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
	  if (margin_ != margin){
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
	
	const String& AxisWidget::getLegend()
	{
		return legend_;
	}
	
	void AxisWidget::setInverseOrientation(bool inverse_orientation)
	{
		if (inverse_orientation_ != inverse_orientation)
		{
	  	inverse_orientation_ = inverse_orientation;
	    update();
	  }
	}
	
	bool AxisWidget::hasInverseOrientation()
	{
		return inverse_orientation_;
	}
	
	void AxisWidget::setLegend(const String& legend)
	{
		legend_ = legend;
		if (!show_legend_)
  	{
		  setToolTip(legend_.c_str());
 	 	}
	}
	
	void AxisWidget::setTickLevel(UInt level)
	{
		if (level==1 || level==2)
		{
			tick_level_ = level;
		}
	}

	const AxisWidget::GridVector& AxisWidget::gridLines()
	{
		return grid_line_;
	}
	
	void AxisWidget::getShortenedNumber_(QString& short_num, DoubleReal number)
	{
		if (!allow_short_numbers_)
		{
			short_num = QString("%1").arg(number);
		}
		else if (number < 1000.0)
		{
			short_num = QString("%1").arg(number);
		}
		else if (number < 1000000.0)
		{
			short_num = QString("%1k").arg(Math::roundDecimal(number/1000.0, -1));
		}
		else if (number < 1000000000.0)
		{
			short_num = QString("%1M").arg(number/1000000.0);
		}
		else
		{
			short_num = QString("%1G").arg(number/1000000000.0);
		}
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

