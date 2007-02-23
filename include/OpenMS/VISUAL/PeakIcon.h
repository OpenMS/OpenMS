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

#ifndef OPENMS_VISUAL_PEAKICON_H
#define OPENMS_VISUAL_PEAKICON_H

class QPainter;
class QRect;

namespace OpenMS 
{
	/**
		@brief Class for drawing icons with a QPainter.
		
		
		
		@ingroup Visual
	*/
	class PeakIcon
	{
		public:
			/// Icon names
			enum Icon
			{
				IT_NOICON,
				IT_ELLIPSE,
				IT_TRIANGLE,
				IT_ASTERIX, 
				IT_SQUARE
			};
			
			/// Draws an icon into a boundingbox
      static void drawIcon(Icon icon, QPainter& painter, const QRect& r);
			/// Draws an ellipse into a boundingbox
      static void drawEllipse(QPainter& painter, const QRect& r);
      /// Draws a triangle into a boundingbox
      static void drawTriangle(QPainter& painter, const QRect& r);
      /// Draws an asterix into a boundingbox
      static void drawAsterix(QPainter& painter, const QRect& r);
      /// Draws a rectangle into a boundingbox
      static void drawRectangle(QPainter& painter, const QRect& r);   
        
		protected:
			/// Constructor is protexcted as all methods are static
			PeakIcon();				
	};
}
#endif
