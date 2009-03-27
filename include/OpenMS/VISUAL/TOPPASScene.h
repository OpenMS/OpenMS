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

#ifndef OPENMS_VISUAL_TOPPASSCENE_H
#define OPENMS_VISUAL_TOPPASSCENE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QGraphicsScene>

namespace OpenMS
{
	class TOPPASVertex;
	class TOPPASEdge;
	
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
      
      typedef std::vector<TOPPASEdge*> EdgeContainer;
			typedef EdgeContainer::iterator EdgeIterator;
			typedef EdgeContainer::const_iterator ConstEdgeIterator;
			typedef std::vector<TOPPASVertex*> VertexContainer;
			typedef VertexContainer::iterator VertexIterator;
			typedef VertexContainer::const_iterator ConstVertexIterator;
			
			/// Standard constructor
			TOPPASScene();
		
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
			
		public slots:
		
			/// Called when an item is clicked
			void itemClicked();
			/// Called when an item is double-clicked
			void itemDoubleClicked();
			/// Called when the position of the hovering edge changes
			void updateHoveringEdgePos(const QPointF& new_pos);
			/// Called when a new out edge is supposed to be created
			void addHoveringEdge(const QPointF& pos);
			
		protected:
		
			/// The current action mode
			ActionMode action_mode_;
			/// The list of all vertices
			VertexContainer vertices_;
			/// The list of all edges
			EdgeContainer edges_;
			/// The hovering edge which is currently being created
			TOPPASEdge* hover_edge_;
	};

}

#endif
