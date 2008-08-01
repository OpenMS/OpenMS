// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
	class String;
	
	/**
		@brief Widget which can visualize a histogram.
		
		@image html HistogramWidget.png
		
		It can also be used to define a left and right boundary inside the values.
		
		@ingroup Visual
	*/
	class HistogramWidget
		: public QWidget
	{
		Q_OBJECT
		
		public:
			/// Constructor
			HistogramWidget(const Math::Histogram<UInt,Real>& distribution, QWidget* parent = 0);
			/// Destructor
			virtual ~HistogramWidget();
			
			/// Returns the value f the lower splitter
			Real getLeftSplitter();
			/// Returns the value of the upper splitter
			Real getRightSplitter();
			/// set axis legends
			void setLegend(const String& legend);

		public slots:
			/// Shows the splitters if @p on is true. Hides them otherwise.
			void showSplitters(bool on);
			/// Sets the value of the right splitter
			void setRightSplitter(Real pos);
			/// Sets the value of the left splitter
			void setLeftSplitter(Real pos);

		protected:
			/// the histogram to display
			Math::Histogram<UInt,Real> dist_;
			/// Flag that indicates if splitters are shown
			bool show_splitters_;
			/// value of the right splitter
			Real left_splitter_;
			/// value of the right splitter
			Real right_splitter_;
			/// the splitter that is currently dragged (0=none, 1=left, 2=right)
			UInt moving_splitter_;
			/// x axis
			AxisWidget *bottom_axis_;
			/// margin around plot
			UInt margin_;
			/// internal buffer for the double buffering
			QPixmap buffer_;
			/// repaints the contents to the buffer and calls update()
			void invalidate_();
			///@name reimplemented Qt events
			//@{
			void paintEvent( QPaintEvent * );
			void mousePressEvent( QMouseEvent *);
			void mouseReleaseEvent( QMouseEvent *);
			void mouseMoveEvent( QMouseEvent *);
			void resizeEvent( QResizeEvent *);
			//@}
	};
} // namespace OpenMS


#endif
