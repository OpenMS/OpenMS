// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_MULTIGRADIENTSELECTOR_H
#define OPENMS_VISUAL_MULTIGRADIENTSELECTOR_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/MultiGradient.h>

//QT
#include <QtGui/QWidget>
class QPaintEvent;
class QMouseEvent;
class QKeyEvent;
class QContextMenuEvent;

//std lib
#include <vector>

namespace OpenMS 
{

	/**
		@brief A widget witch allows constructing gradients of multiple colors.
		
		@image html MultiGradientSelector.png
		
		The above example image shows a MultiGradientSelector.
		
		@ingroup Visual
	*/
	class OPENMS_DLLAPI MultiGradientSelector
		: public QWidget
	{
		Q_OBJECT
		public:
			///Constructor
			MultiGradientSelector( QWidget* parent = 0);
			///Desctructor
			~MultiGradientSelector();
			
			///returns a const reference to the gradient
			const MultiGradient& gradient() const;
			///returns a mutable reference to the gradient
			MultiGradient& gradient();
		
			/// sets the interploation mode
			void setInterpolationMode(MultiGradient::InterpolationMode mode);
			/// returns the interpolaion mode
			MultiGradient::InterpolationMode getInterpolationMode() const;
		
			public slots:
			/// sets what interpolation mode is used
			void stairsInterpolation(bool state);
	
		protected:
			///@name Reimpelented Qt events
			//@{
			void paintEvent(QPaintEvent* e);
			void mousePressEvent(QMouseEvent* e);
			void mouseMoveEvent(QMouseEvent* e);		
			void mouseReleaseEvent(QMouseEvent* e);
			void mouseDoubleClickEvent (QMouseEvent* e);
			void keyPressEvent(QKeyEvent* e);
			void contextMenuEvent(QContextMenuEvent* e);
			//@}
			
			// the actual gradient
			MultiGradient gradient_;
			
			// border margin
			Int margin_;
			// height of the gradient area
			Int gradient_area_width_;
			// heigth of the lever area
			Int lever_area_height_;
			
			//position (0-100) in the vector of the selected lever
			Int selected_;
			//color of the selected lever
			QColor selected_color_;  
			
			//stores if the left button is pressed
			bool left_button_pressed_;
	};

}
#endif // OPENMS_VISUAL_MULTIGRADIENTSELECTOR_H

