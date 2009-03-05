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
// $Maintainer: Johannes Junker $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASEDGE_H
#define OPENMS_VISUAL_TOPPASEDGE_H

#include <OpenMS/VISUAL/TOPPASVertex.h>

#include <QtGui/QGraphicsItem>

namespace OpenMS
{
	class OPENMS_DLLAPI TOPPASEdge
		: public QGraphicsItem
	{		
		public:
			
			/// Constructor
			TOPPASEdge(TOPPASVertex* from, TOPPASVertex* to);
			
			/// Destructor
			virtual ~TOPPASEdge();
			
			/// Returns the bounding rectangle of this item
			QRectF boundingRect() const;
			
			/// Paints the item
			void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
		
		protected:
		
			/// Pointer to the source of this edge
			TOPPASVertex* from_;
			/// Pointer to the target of this edge
			TOPPASVertex* to_;
	};
}

#endif
