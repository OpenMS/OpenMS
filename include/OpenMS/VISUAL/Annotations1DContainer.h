// -*- Mode: C++; tab-width: 2; -*-
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

#ifndef OPENMS_VISUAL_ANNOTATIONS1DCONTAINER_H
#define OPENMS_VISUAL_ANNOTATIONS1DCONTAINER_H

#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/KERNEL/PeakIndex.h>

#include <list>

#include <QtCore/QPoint>
#include <QtCore/QObject>
#include <QtCore/QRectF>
#include <QtGui/QPainter>

namespace OpenMS
{

	class Annotations1DContainer
		: public std::list<Annotation1DItem*>
	{
		public:
		
			/// Iterator for the 1D annotations
			typedef std::list<Annotation1DItem*>::iterator Iterator;
			/// Const iterator for the 1D annotations
			typedef std::list<Annotation1DItem*>::const_iterator ConstIterator;
			/// Type of the Points
			typedef DPosition<2> PointType;
			/// Coordinate type
			typedef DoubleReal CoordinateType;
			
			/// Default constructor
			Annotations1DContainer();
			/// Copy constructor
			Annotations1DContainer(const Annotations1DContainer& rhs);
			/// Destructor
			virtual ~Annotations1DContainer();	
			/// Assignment operator
			Annotations1DContainer& operator= (const Annotations1DContainer& rhs);
			
			/** @brief Returns a pointer to the item at @p pos, or 0, if not existent
					
					If more than one item's bounding box encloses @p pos , the one in the
					foreground is returned.
			*/
			Annotation1DItem* getItemAt(const QPoint& pos) const;
			/// Selects the item at @p pos on the canvas, if it exists.
			void selectItemAt(const QPoint& pos);
			/// Deselects the item at @p pos on the canvas, if it exists.
			void deselectItemAt(const QPoint& pos);
			/// Selects all items
			void selectAll();
			/// Deselects all items
			void deselectAll();
			/// Removes the selected items
			void removeSelectedItems();
	};
	
} // namespace

#endif
