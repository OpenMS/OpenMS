// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_AXISWIDGET_H
#define OPENMS_VISUAL_AXISWIDGET_H

// QT
#include <qwidget.h>
#include <qpainter.h>
class QPixmap;

// STL
#include <vector>
#include <math.h>

// OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/MappingInfo.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{
	/**
		@brief Widget that represents an axis of a graph.
		
		Additional to ticks and tick values a label e.g. the unit can be displayed.
		It support both linear and logarithmic scale.
	
		@ingroup Visual
	*/
	class AxisWidget : public QWidget
	{
		Q_OBJECT
		
		public:
	
			///Vector of vector of doubles that defines the grid
			typedef std::vector<std::vector<double> > GridVector;
			/// Where the axis is placed
			static enum {TOP, BOTTOM, LEFT, RIGHT} ALIGNMENT_ENUM;
	
	
			/// constructor
			AxisWidget(UnsignedInt alignment, const char* legend="", QWidget* parent = 0, const char* name = "AxisWidget", WFlags f = 0);
			/// destructor
			virtual ~AxisWidget();
	
			/// sets the margin on the top/right side (default is 0)
			void setMargin(UnsignedInt size);
	
			/// returns the margin
			UnsignedInt margin();
	
			/// enable the display of the legend (default true)
			void showLegend(bool show_legend);
	
			// returns true if legend is shown
			bool isLegendShown() const;
	
			/// return constant reference to grid_line_
			/// gridlines are calculated by setAxisBounds or set by the user using setGridLines
			const GridVector& gridLines();
	
			/// set user-defined gridlines (instead of using setAxisBounds to calculate GrdiLines)
			void setGridLines(std::vector<double>&);
	
			/// sets the axis to logarithmic scale
			void setLogScale(bool is_log);
	
			/// returns true if the axis has logarithmic scal
			bool isLogScale();
	
			/// sets the legend text
			void setLegend(const std::string& legend);
	
			/// returns the actual legend text
			std::string getLegend();

			/// set true to display the axis label in inverse order (left to right or bottom to top)
			void setInverseOrientation(bool inverse_orientation);
			bool hasInverseOrientation();
	
			/// see QWidget
			virtual QSize sizeHint () const;
	
	    inline double getAxisMinimum() const
	    {
	    	return min_;
	    };
	    inline double getAxisMaximum() const 
	    {
	    	return max_;
	    };
	
			inline void setPenWidth(int p)
			{
				pen_width_ = p;
				invalidate_();
			}
	
		signals:
			void sendAxisMouseMovement();
	
		public slots:
			///sets min/max of the axis
			void setAxisBounds(double min, double max);
			/// see QWidget
			void resize(int,int);
			/// see QWidget
			void resize(QSize);
			/// set maximum number of tick levels (default=3), reduce number in case of small axis
			bool setTickLevel(UnsignedInt level);
			void mouseMoveEvent( QMouseEvent *e);
	
		protected:
			/// Vector that defines the position of the ticks/gridlines and the shown values on axis
			GridVector grid_line_;
	
	    /// The MappingInfo object holds information about switched/mirrored axis
	    MappingInfo mapping_info_;
	
			/// format of axis scale (linear or logarithmic)
			bool is_log_;
			/// display of legend enabled or not
			bool show_legend_;
			/// Position of the axis (right, left, top, down as defined in ALIGNMENT_ENUM)
			UnsignedInt alignment_;
			/// true if axis label are displayed in inverse order (left to right or bottom to top)
			bool inverse_orientation_;
			/// margin of axis
			UnsignedInt margin_;
			// interval of values on axis
			double min_, max_;
			/// text/unit on axis
			std::string legend_;
			/// maximum number of tick levels (default=3)
			UnsignedInt tick_level_;
	
			///painting buffer
			QPixmap* buffer_;
			///the painter
			QPainter painter_;
			/// drawing thicker lines (e.g. in printing) leads to better results
			UnsignedInt pen_width_;
	
			/// see QWidget
			void resizeEvent( QResizeEvent * );
			/// see QWidget
			void paintEvent( QPaintEvent * );
	
			double grid_line_dist_;
	
		private:
			///repaint the content of the widget to the internal pixmap if content or size changed.
			void invalidate_();
	
			/// Scale axis values to correct value (i.e. reverse log, unit conversion)
			inline double scale_(double x)
			{
				return (is_log_)? Math::round_decimal(pow(x,10),-8) : Math::round_decimal(x,-8);
			}
	
			inline int probeFont_(QString probe, double width, double height, int index=0)
			{
				int probe_font = 10;
				painter_.setFont(QFont("courier",probe_font));
				QRect prb_bound = QFontMetrics(painter_.font()).boundingRect(probe);
				double ratio1 = prb_bound.width()/width;
				double ratio2 = prb_bound.height()/height;
				double res = (ratio1 > ratio2)? ratio1:ratio2;
				if (index==1)
				{
					return static_cast<int>(probe_font/res*3.0/2.4);
				}
				else if (index==2)
				{
					return static_cast<int>(probe_font/res*3.0/2);
				}
				else
				{ 
					return static_cast<int>(probe_font/res);
				}
			}
 	};
} // namespace OpenMS


#endif

