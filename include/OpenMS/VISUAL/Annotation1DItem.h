// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_ANNOTATION1DITEM_H
#define OPENMS_VISUAL_ANNOTATION1DITEM_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <QtCore/QRectF>
#include <QtCore/QString>

class QPainter;

namespace OpenMS
{
	class Spectrum1DCanvas;

	/** @brief An abstract class acting as an interface for the different 1d annotation items.
	
			This is an abstract polymorphic class which acts as an interface between its
			subclasses and all containers and methods that contain or handle Annotation1DItem
			objects.
			
			If you want to add new kinds of annotation items, inherit this class,
			implement the pure virtual methods, and add everything else the item should
			have or be capable of.
			
	*/
	class Annotation1DItem
	{
	
		public:
			
			/// Type of the Points
			typedef DPosition<2> PointType;
			/// Intensity type
			typedef Real IntensityType;
			/// Coordinate type
			typedef DoubleReal CoordinateType;			
			
			/// Destructor
			virtual ~Annotation1DItem();
			
			/// Returns the current bounding box of this item on the canvas where it has last been drawn
			const QRectF& boundingBox() const;
			/// Returns true if this item is currently selected on the canvas, else false
			bool isSelected() const;
			/// Sets whether this item is currently selected on the canvas or not
			void setSelected(bool selected);
			/// Sets the text of the item
			void setText(const QString& text);
			/// Returns the text of the item
			const QString& getText() const;
			/// Ensures that the item has coordinates within the visible area of the canvas
			virtual void ensureWithinDataRange(Spectrum1DCanvas* const canvas) = 0;
			/// Draws the item on @p painter
			virtual void draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped = false) = 0;
			/// Moves the item; behaviour depends on item type and is implemented in the subclasses
			virtual void move(const PointType& delta) = 0;
			
		protected:
		
			/// Constructor
			Annotation1DItem(const QString& text);
			/// Copy constructor
			Annotation1DItem(const Annotation1DItem& rhs);
			/// Draws the bounding_box_
			void drawBoundingBox_(QPainter& painter);
			
			/// The current bounding box of this item on the canvas where it has last been drawn
			QRectF bounding_box_;
			/// Determines whether this item is currently selected on the canvas
			bool selected_;
			/// The displayed text
			QString text_;
			
	};
} // namespace OpenMS

#endif
