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

#ifndef OPENMS_VISUAL_TOPPASSCENE_H
#define OPENMS_VISUAL_TOPPASSCENE_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>

#include <QtGui/QGraphicsScene>

namespace OpenMS
{
	class TOPPASVertex;
	class TOPPASEdge;
	
	/**
		@brief A container for all visual items of a TOPPAS workflow
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASScene
		:	public QGraphicsScene
	{
		Q_OBJECT
		
		public:
			
			enum ActionMode
      {
      	AM_NEW_EDGE,
      	AM_MOVE
      };
      
      typedef QList<TOPPASEdge*> EdgeContainer;
			typedef EdgeContainer::iterator EdgeIterator;
			typedef EdgeContainer::const_iterator ConstEdgeIterator;
			typedef QList<TOPPASVertex*> VertexContainer;
			typedef VertexContainer::iterator VertexIterator;
			typedef VertexContainer::const_iterator ConstVertexIterator;
			
			/// Constructor
			TOPPASScene(QObject* parent);
			
			/// Destructor
			virtual ~TOPPASScene();
			
			/// Adds a vertex
			void addVertex(TOPPASVertex* tv);
			/// Adds an edge
			void addEdge(TOPPASEdge* te);
			/// Sets the action mode
			void setActionMode(ActionMode mode);
			/// Returns the action mode
			ActionMode getActionMode();
			/// Returns begin() iterator of all vertices
			VertexIterator verticesBegin();
			/// Returns end() iterator of all vertices
			VertexIterator verticesEnd();
			/// Returns begin() iterator of all edges
			EdgeIterator edgesBegin();
			/// Returns end() iterator of all edges
			EdgeIterator edgesEnd();
			/// Removes all currently selected edges and vertices
			void removeSelected();
			/// Updates all edge colors (color of green and yellow edges can change when edges are added/removed)
			void updateEdgeColors();
			/// Runs the pipeline
			void runPipeline();
			
		public slots:
		
			/// Called when an item is clicked
			void itemClicked();
			/// Called when an item is double-clicked
			void itemDoubleClicked();
			/// Called when the position of the hovering edge changes
			void updateHoveringEdgePos(const QPointF& new_pos);
			/// Called when a new out edge is supposed to be created
			void addHoveringEdge(const QPointF& pos);
			/// Called when the new edge is being "released"
			void finishHoveringEdge();
			
		protected:
			
			/// The current action mode
			ActionMode action_mode_;
			/// The list of all vertices
			VertexContainer vertices_;
			/// The list of all edges
			EdgeContainer edges_;
			/// The hovering edge which is currently being created
			TOPPASEdge* hover_edge_;
			/// The current potential target vertex of the hovering edge
			TOPPASVertex* potential_target_;
			
			/// Returns the vertex in the foreground at position @p pos , if existent, otherwise 0.
			TOPPASVertex* getVertexAt_(const QPointF& pos);
			/// Returns whether an edge between node u and v would be allowed
			bool isEdgeAllowed_(TOPPASVertex* u, TOPPASVertex* v);
			/// DFS helper method. Returns true, if a back edge has been discovered
			bool dfsVisit_(TOPPASVertex* vertex);
	};

}

#endif
