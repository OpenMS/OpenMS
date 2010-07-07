// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

// ------------- DEBUGGING ----------------

// ---- Uncomment to enable debug mode ----
//#define TOPPAS_DEBUG
// ----------------------------------------

#ifdef TOPPAS_DEBUG
#define __DEBUG_BEGIN_METHOD__	{																																													\
																		for (int dbg_indnt_cntr = 0; dbg_indnt_cntr < global_debug_indent_; ++dbg_indnt_cntr)	\
																		{																																											\
																			std::cout << "  ";																																	\
																		}																																											\
																		std::cout << "BEGIN [" << topo_nr_ << "] " << __PRETTY_FUNCTION__ << std::endl;				\
																		++global_debug_indent_;																																\
																}

#define __DEBUG_END_METHOD__	{																																														\
																		--global_debug_indent_;																																\
																		if (global_debug_indent_ < 0) global_debug_indent_ = 0;																\
																		for (int dbg_indnt_cntr = 0; dbg_indnt_cntr < global_debug_indent_; ++dbg_indnt_cntr)	\
																		{																																											\
																			std::cout << "  ";																																	\
																		}																																											\
																		std::cout << "END [" << topo_nr_ << "] " << __PRETTY_FUNCTION__ << std::endl;					\
																}
#else
#define __DEBUG_BEGIN_METHOD__ {}
#define __DEBUG_END_METHOD__ {}
#endif

// ----------------------------------------

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>
#include <QtGui/QGraphicsSceneMouseEvent>
#include <QtGui/QGraphicsSceneContextMenuEvent>
#include <QtGui/QGraphicsItem>
#include <QtCore/QProcess>
#include <QtGui/QMenu>

namespace OpenMS
{
	class TOPPASEdge;

	/**
		@brief The base class of the different vertex classes.
		
		This class contains the common functionality (such as
		event handling for mouse clicks and drags) and holds the common
		members of all different kinds of vertices (e.g., containers
		for all in and out edges, the vertex ID, the number of a 
		topological sort of the whole graph, etc.)
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASVertex
		: public QObject,
			public QGraphicsItem
	{
		Q_OBJECT
		
		public:
			
			/// The container for in/out edges
			typedef QList<TOPPASEdge*> EdgeContainer;
			/// A mutable iterator for in/out edges
			typedef EdgeContainer::iterator EdgeIterator;
			/// A const iterator for in/out edges
			typedef EdgeContainer::const_iterator ConstEdgeIterator;
			
			/// The color of a vertex during depth-first search
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
			/// Sets whether all tools in the subtree below this node are finished
			void setSubtreeFinished(bool b);
			/// Returns whether the vertex has been marked already (during topological sort)
			bool isTopoSortMarked();
			/// (Un)marks the vertex (during topological sort)
			void setTopoSortMarked(bool b);
			/// Returns the topological sort number
			UInt getTopoNr();
			/// Sets the topological sort number (overridden in tool and output vertices)
			virtual void setTopoNr(UInt nr);
			/// Checks if the tools in the subtree below this node have finished and if yes, propagates this upwards
			virtual void checkIfSubtreeFinished();
			/// Returns if all mergers further upstream in the pipeline have finished merging already
			virtual bool areAllUpstreamMergersFinished();
			/// Resets the status
			virtual void reset(bool reset_all_files = false);
			/// Returns whether all tools in the subtree below this node are finished
			virtual bool isSubtreeFinished();
			/// Resets the subtree below this node
			virtual void resetSubtree(bool including_this_node = true);
			/// Recursive sanity check for mergers and tools
			virtual void checkListLengths(QStringList& unequal_per_round, QStringList& unequal_over_entire_run);
			/// Used in sanity check, returns the overall number of files at this node over all merging rounds
			int getScFilesTotal();
			/// Used in sanity check, returns the number of files at this node per merging round
			int getScFilesPerRound();
			/// Used in sanity check, returns if this node has already been checked
			bool isScListLengthChecked();
			/// Returns whether this node has already been started by an input file node
			bool isAlreadyStarted();
			/// Sets whether this node has already been started by an input file node
			void setAlreadyStarted(bool b);
			/// Marks this node (and everything further downstream) as unreachable. Overridden behavior in mergers.
			virtual void markUnreachable();
			/// Returns whether this node is reachable
			bool isReachable();
		
		public slots:
		
			/// Called by an incoming edge when it has changed
			virtual void inEdgeHasChanged();
		
		signals:
			
			/// Emitted when this item is clicked
			void clicked();
			/// Emitted when this item is released
			void released();
			/// Emitted when the position of the hovering edge changes
			void hoveringEdgePosChanged(const QPointF& new_pos);
			/// Emitted when a new out edge is supposed to be created
			void newHoveringEdge(const QPointF& pos);
			/// Emitted when the mouse is released after having dragged a new edge somewhere
			void finishHoveringEdge();
			/// Emitted when something has changed
			void somethingHasChanged();
			/// Emitted when the item is dragged
			void itemDragged(qreal dx, qreal dy);
			
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
			/// "marked" flag for topological sort
			bool topo_sort_marked_;
			/// The number in a topological sort of the entire graph
			UInt topo_nr_;
			/// Indicates whether all tools in the subtree below this node are finished
			bool subtree_finished_;
			/// Used during sanity check, stores the number of files at this node per merging round
			int sc_files_per_round_;
			/// Used during sanity check, stores the overall number of files at this node over all merging rounds
			int sc_files_total_;
			/// Used during sanity check, stores if this node has been checked already
			bool sc_list_length_checked_;
			/// Stores if the file names for the current call are already known
			bool files_known_;
			/// Stores if this node has already been started by an input file node
			bool already_started_;
			/// Indicates whether this node is reachable (i.e. there is an input node somewhere further upstream)
			bool reachable_;
			
			#ifdef TOPPAS_DEBUG
			// Indentation level for nicer debug output
			static int global_debug_indent_;
			#endif
			
			///@name reimplemented Qt events
      //@{
      void mouseReleaseEvent(QGraphicsSceneMouseEvent* e);
      void mousePressEvent(QGraphicsSceneMouseEvent* e);
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
      void mouseMoveEvent(QGraphicsSceneMouseEvent* e);
      void contextMenuEvent(QGraphicsSceneContextMenuEvent* event);
			//@}
			
			/// Moves the target pos of the edge which is just being created to @p pos
			virtual void moveNewEdgeTo_(const QPointF& pos);
			/// Returns a three character string (i.e. 001 instead of 1) for the given @p number
			String get3CharsNumber_(UInt number);
			/// Removes the specified directory (absolute path). Returns true if successful.
			bool removeDirRecursively_(const QString& dir_name);
			
			/// Displays the debug output @p message, if TOPPAS_DEBUG is defined
			void debugOut_(const String&
			#ifdef TOPPAS_DEBUG
			message
			#endif
			)
			{
				#ifdef TOPPAS_DEBUG
				for (int i = 0; i < global_debug_indent_; ++i)
				{
					std::cout << "  ";
				}
				std::cout << "[" << topo_nr_ << "] " << message << std::endl;
				#endif
			}
	};
}

#endif
