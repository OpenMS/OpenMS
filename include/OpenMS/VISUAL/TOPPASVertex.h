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

#ifndef OPENMS_VISUAL_TOPPASVERTEX_H
#define OPENMS_VISUAL_TOPPASVERTEX_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QGraphicsItem>

namespace OpenMS
{
	class TOPPASEdge;

	class OPENMS_DLLAPI TOPPASVertex
		: public QObject,
			public QGraphicsItem
	{
		Q_OBJECT
		
		public:
			
			enum VertexType
			{
				VT_SOURCE,
				VT_TARGET,
				VT_TOOL
			};
			
			typedef std::vector<TOPPASEdge*> EdgeContainer;
			typedef EdgeContainer::iterator EdgeIterator;
			typedef EdgeContainer::const_iterator ConstEdgeIterator;
			
			
			/// Constructor
			TOPPASVertex(const String& name, const String& type = "", VertexType vt = VT_TOOL);
			
			/// Destructor
			virtual ~TOPPASVertex();
			
			/// Returns the name of the tool
			const String& getName();
			/// Returns the bounding rectangle of this item
			QRectF boundingRect() const;
			/// Returns a more precise shape
			QPainterPath shape () const;
			/// Paints the item
			void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			/// Returns begin() iterator of outgoing edges
			EdgeIterator outEdgesBegin();
			/// Returns end() iterator of outgoing edges
			EdgeIterator outEdgesEnd();
			/// Returns begin() iterator of incoming edges
			EdgeIterator inEdgesBegin();
			/// Returns end() iterator of incoming edges
			EdgeIterator inEdgesEnd();

		signals:
			
			/// Emitted when this item is clicked
			void clicked();
			/// Emitted when this item is double-clicked
			void doubleClicked();
			/// Emitted when the position of the hovering edge changes
			void hoveringEdgePosChanged(const QPointF& new_pos);
			/// Emitted when a new out edge is supposed to be created
			void newHoveringEdge(const QPointF& pos);
			
		protected:
			
			/// The name of the tool
			String name_;
			/// The type of the tool, or "" if it does not have a type
			String type_;
			/// The type of this vertex
			VertexType vertex_type_;
			/// The list of outgoing edges
			EdgeContainer out_edges_;
			/// The list of incoming edges
			EdgeContainer in_edges_;
			/// Indicates whether a new out edge is currently being created
			bool edge_being_created_;
			/// The color of the pen
			QColor pen_color_;
			/// The color of the brush
			QColor brush_color_;
			
			///@name reimplemented Qt events
      //@{
      void mouseReleaseEvent(QGraphicsSceneMouseEvent* e);
      void mousePressEvent(QGraphicsSceneMouseEvent* e);
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
      void mouseMoveEvent(QGraphicsSceneMouseEvent* e);
			//@}
			
	};
}

#endif
