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
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASVERTEX_H
#define OPENMS_VISUAL_TOPPASVERTEX_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>
#include <QtGui/QGraphicsSceneMouseEvent>
#include <QtGui/QGraphicsItem>

namespace OpenMS
{
	class TOPPASEdge;

	/**
		@brief The base class of all different kinds of vertices
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASVertex
		: public QObject,
			public QGraphicsItem
	{
		Q_OBJECT
		
		public:
		
			typedef QList<TOPPASEdge*> EdgeContainer;
			typedef EdgeContainer::iterator EdgeIterator;
			typedef EdgeContainer::const_iterator ConstEdgeIterator;
			
			enum DFS_COLOR
			{
				DFS_WHITE,
				DFS_GRAY,
				DFS_BLACK
			};
			
			/// Default Constructor
			TOPPASVertex();
			/// Copy constructor
			TOPPASVertex(const TOPPASVertex& rhs);
			/// Destructor
			virtual ~TOPPASVertex();
			/// Assignment operator
			TOPPASVertex& operator= (const TOPPASVertex& rhs);
			
			/// Returns the bounding rectangle of this item
			virtual QRectF boundingRect() const = 0;
			/// Returns a more precise shape
			virtual QPainterPath shape () const = 0;
			/// Paints the item
			virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) = 0;
			/// Returns begin() iterator of outgoing edges
			EdgeIterator outEdgesBegin();
			/// Returns end() iterator of outgoing edges
			EdgeIterator outEdgesEnd();
			/// Returns begin() iterator of incoming edges
			EdgeIterator inEdgesBegin();
			/// Returns end() iterator of incoming edges
			EdgeIterator inEdgesEnd();
			/// Returns the number of incoming edges
			Size incomingEdgesCount();
			/// Returns the number of outgoing edges
			Size outgoingEdgesCount();
			/// Adds an incoming edge
			void addInEdge(TOPPASEdge* edge);
			/// Adds an outgoing edge
			void addOutEdge(TOPPASEdge* edge);
			/// Removes an incoming edge
			void removeInEdge(TOPPASEdge* edge);
			/// Removes an outedge
			void removeOutEdge(TOPPASEdge* edge);
			/// Returns the DFS color of this node
			DFS_COLOR getDFSColor();
			/// Sets the DFS color of this node
			void setDFSColor(DFS_COLOR color);
			/// Returns the DFS parent of this node
			TOPPASVertex* getDFSParent();
			/// Sets the DFS parent of this node
			void setDFSParent(TOPPASVertex* parent);
			
		signals:
			
			/// Emitted when this item is clicked
			void clicked();
			/// Emitted when this item is double-clicked
			void doubleClicked();
			/// Emitted when the position of the hovering edge changes
			void hoveringEdgePosChanged(const QPointF& new_pos);
			/// Emitted when a new out edge is supposed to be created
			void newHoveringEdge(const QPointF& pos);
			/// Emitted when the mouse is released after having dragged a new edge somewhere
			void finishHoveringEdge();
			
		protected:
			
			/// The list of incoming edges
			EdgeContainer in_edges_;
			/// The list of outgoing edges
			EdgeContainer out_edges_;
			/// Indicates whether a new out edge is currently being created
			bool edge_being_created_;
			/// The color of the pen
			QColor pen_color_;
			/// The color of the brush
			QColor brush_color_;
			/// The DFS color of this node
			DFS_COLOR dfs_color_;
			/// The DFS parent of this node
			TOPPASVertex* dfs_parent_;
			
			///@name reimplemented Qt events
      //@{
      void mouseReleaseEvent(QGraphicsSceneMouseEvent* e);
      void mousePressEvent(QGraphicsSceneMouseEvent* e);
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
      void mouseMoveEvent(QGraphicsSceneMouseEvent* e);
			//@}
			
			/// Moves the target pos of the edge which is just being created to @p pos
			virtual void moveNewEdgeTo_(const QPointF& pos);
			
	};
}

#endif
