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
// $Maintainer: Marc Sturm  $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_RUBBERBAND_H
#define OPENMS_VISUAL_RUBBERBAND_H

#include <qpainter.h>


namespace OpenMS 
{
	/**
		@brief Rubber band to select ranges on the screen.
		
		Implementation as described in the "C++ GUI Programming with Qt 3" book.
		
		@todo replace by the same as used in 2D and 3D? (Marc)
		
		@ingroup Visual
	*/
	class RubberBand
	{
		public:
			RubberBand()
				:is_shown_(false),
				rect_(QRect()) 
			{
				
			}
			
			bool isShown() 
			{ 
				return is_shown_; 
			}
			
			void show() 
			{ 
				is_shown_ = true; 
			}
			
			void hide() 
			{ 
				is_shown_ = false;
			}

      // schedule rectangle area in parent widget for redraw
			void updateRegion(QWidget* parent)
			{
				QRect r = rect_.normalize();
				parent->update(r);
			}

      // draw rubber band using XorRop (xored cyan will appear red on white background)
			void draw(QPainter& painter)
			{
				painter.save();
				painter.setBrush(Qt::cyan);
				painter.setRasterOp(Qt::XorROP);
				painter.drawRect(rect_.normalize());
				painter.restore();
			}

      QRect getRect() 
      { 
      	return rect_.normalize(); 
      }

      void setRect(const QRect& r) 
      { 
      	rect_ = r;	
      }
      
			void setBottomRight(const QPoint& p) 
			{ 
				rect_.setBottomRight(p);
			}
			
			void setTopLeft(const QPoint& p) 
			{ 
				rect_.setTopLeft(p);
			}

    protected:
			bool is_shown_;
			QRect rect_;
	};
}
#endif
