// -*- mode: C++; tab-width: 2; -*-
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
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASEDGE_H
#define OPENMS_VISUAL_TOPPASEDGE_H

#include <OpenMS/config.h>

#include <QtGui/QGraphicsItem>

namespace OpenMS
{
	class TOPPASVertex;
	class TOPPASToolVertex;
	class TOPPASInputFileListVertex;
  

	/**
		@brief An edge representing a data flow in TOPPAS
		
		Like all TOPPASVertex classes, TOPPASEdge is a subclass of QGraphicsItem and thus implements
		methods to draw itself and to react on incoming events such as mouse clicks. It holds
		the data needed to represent an edge between two vertices of a TOPPAS workflow.
		
		@ingroup TOPPAS_elements
	*/
	class OPENMS_GUI_DLLAPI TOPPASEdge
		: public QObject,
			public QGraphicsItem
	{
		Q_OBJECT
		
		public:
			
			/// The status of this edge
			enum EdgeStatus
			{
				ES_VALID,
				ES_NO_TARGET_PARAM,
				ES_NO_SOURCE_PARAM,
				ES_FILE_EXT_MISMATCH,
				ES_MERGER_EXT_MISMATCH,
				ES_MERGER_WITHOUT_TOOL,
				ES_NOT_READY_YET,
        ES_TOOL_API_CHANGED,
				ES_UNKNOWN
			};
			
			/// Standard constructor
			TOPPASEdge();
			/// Constructor
			TOPPASEdge(TOPPASVertex* from, const QPointF& hover_pos);
			/// Copy constructor
			TOPPASEdge(const TOPPASEdge& rhs);
			/// Destructor
			virtual ~TOPPASEdge();
			/// Assignment operator
			TOPPASEdge& operator= (const TOPPASEdge& rhs);
			
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
			/// Sets the source vertex of this edge
			void setSourceVertex(TOPPASVertex* tv);
			/// Sets the target vertex of this edge
			void setTargetVertex(TOPPASVertex* tv);
			/// Returns the source vertex
			TOPPASVertex* getSourceVertex();
			/// Returns the target vertex
			TOPPASVertex* getTargetVertex();
			/// Call this before changing the item geometry
			void prepareResize();
			/// Sets the color
			void setColor(const QColor& color);
			/// Returns the status of this edge
			EdgeStatus getEdgeStatus();
			/// Sets the source output parameter index
			void setSourceOutParam(int out);
			/// Returns the source output parameter index
			int getSourceOutParam();
			/// Sets the target input parameter index
			void setTargetInParam(int in);
			/// Returns the target input parameter index
			int getTargetInParam();
			/// Updates the edge color
			void updateColor();
			/// Emits the somethingHasChanged() signal
			void emitChanged();
			/// Shows the I/O mapping dialog
			void showIOMappingDialog();
		
		public slots:
		
			/// Called by the source vertex when it has changed
			void sourceHasChanged();
		
		signals:
			
			/// Emitted when something has changed
			void somethingHasChanged();
		
		protected:
			
			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
      void contextMenuEvent(QGraphicsSceneContextMenuEvent* event);
			//@}
		
			///@name helper methods of getEdgeStatus()
			//@{
			EdgeStatus getToolToolStatus_(TOPPASToolVertex* source, int source_param_index, TOPPASToolVertex* target, int target_param_index);
		 	EdgeStatus getListToolStatus_(TOPPASInputFileListVertex* source, TOPPASToolVertex* target, int target_param_index);	
			//@}

			/// Returns the point in the @p list that is nearest to @p origin
			QPointF nearestPoint_(const QPointF& origin, const QList<QPointF>& list) const;
			/// Pointer to the source of this edge
			TOPPASVertex* from_;
			/// Pointer to the target of this edge
			TOPPASVertex* to_;
			/// Position of hovering end while edge is being created
			QPointF hover_pos_;
			/// The color
			QColor color_;
			/// The source output parameter index
			int source_out_param_;
			/// The target input parameter index
			int target_in_param_;
	};
}

#endif
