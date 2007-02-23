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


#ifndef OPENMS_VISUAL_HISTOGRAMWIDGET_H
#define OPENMS_VISUAL_HISTOGRAMWIDGET_H

// QT
#include <QtGui/QWidget>
#include <QtGui/QPixmap>
class QPaintEvent;
class QResizeEvent;
class QMouseEvent;

//OpenMS
#include <OpenMS/MATH/STATISTICS/Histogram.h>

namespace OpenMS
{
	class AxisWidget;
	
	/**
		@brief Widget which can visualize a histogram.
		
		It can also be used to define a left and right boundary inside the values.
		
		@ingroup Visual
	*/
	class HistogramWidget
		: public QWidget
	{
		Q_OBJECT
		
		public:
			/// Constructor
			HistogramWidget(const Math::Histogram<UnsignedInt,float>& distribution, QWidget* parent = 0);
			/// Destructor
			virtual ~HistogramWidget();
			
			/// Returns the value f the lower splitter
			float getLeftSplitter();
			/// Returns the value of the upper splitter
			float getRightSplitter();
			/// set x-axis legend; default is 'x'
			void setLegend(const std::string& legend);

		public slots:
			/// Shows the splitters if @p on is true. Hides them otherwise.
			void showSplitters(bool on);
			/// Sets the value of the right splitter
			void setRightSplitter(float pos);
			/// Sets the value of the left splitter
			void setLeftSplitter(float pos);

		protected:
			Math::Histogram<UnsignedInt,float> dist_;
			bool show_splitters_;
			float left_splitter_;
			float right_splitter_;
			/// the splitter that is currently dragged (0=none, 1=left, 2=right)
			UnsignedInt moving_splitter_;
			AxisWidget *bottom_axis_;
			UnsignedInt margin_;
			/// internal buffer for the double buffering
			QPixmap buffer_;
			/// repaints the contents to the buffer and calls update()
			void invalidate_();
			
			/// Factor that is used to for scaling the intensities after the log transformation. 
			/// The higher this factor is the smoother is the diagram.
			float scaling_factor_;
			
			void paintEvent( QPaintEvent * );
			void mousePressEvent( QMouseEvent *);
			void mouseReleaseEvent( QMouseEvent *);
			void mouseMoveEvent( QMouseEvent *);
			void resizeEvent( QResizeEvent *);
	};
} // namespace OpenMS


#endif
