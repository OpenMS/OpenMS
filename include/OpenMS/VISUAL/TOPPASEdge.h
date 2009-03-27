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
		: public QObject,
			public QGraphicsItem
	{
		Q_OBJECT
		
		public:
			
			/// Constructor
			TOPPASEdge(TOPPASVertex* from, const QPointF& hover_pos);
			
			/// Destructor
			virtual ~TOPPASEdge();
			
			/// Returns the bounding rectangle of this item
			QRectF boundingRect() const;
			/// Returns a more precise shape
			QPainterPath shape () const;
			/// Paints the item
			void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			/// Returns the start position of this edge
			QPointF startPos() const;
			/// Returns the end position of this edge
			QPointF endPos() const;
			/// Sets the position of the hovering end while edge is being created
			void setHoverPos(const QPointF& pos);
			/// Sets the target vertex of this edge
			void setTargetVertex(TOPPASVertex* tv);
			
		protected:
		
			/// Pointer to the source of this edge
			TOPPASVertex* from_;
			/// Pointer to the target of this edge
			TOPPASVertex* to_;
			/// Position of hovering end while edge is being created
			QPointF hover_pos_;
	};
}

#endif
