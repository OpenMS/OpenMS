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

#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASScene.h>

#include <QtCore/QDir>

namespace OpenMS
{
	#ifdef TOPPAS_DEBUG
	int TOPPASVertex::global_debug_indent_ = 0;
	#endif
	
	TOPPASVertex::TOPPASVertex()
		:	QObject(),
			QGraphicsItem(),
			in_edges_(),
			out_edges_(),
			edge_being_created_(false),
			pen_color_(),
			brush_color_(),
			dfs_color_(DFS_WHITE),
			dfs_parent_(0),
			id_(0),
			topo_sort_marked_(false),
			topo_nr_(0),
			subtree_finished_(false)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
		setZValue(42);
	}
	
	TOPPASVertex::TOPPASVertex(const TOPPASVertex& rhs)
		:	QObject(),
			QGraphicsItem(),
			in_edges_(rhs.in_edges_),
			out_edges_(rhs.out_edges_),
			edge_being_created_(rhs.edge_being_created_),
			pen_color_(rhs.pen_color_),
			brush_color_(rhs.brush_color_),
			dfs_color_(rhs.dfs_color_),
			dfs_parent_(rhs.dfs_parent_),
			id_(rhs.id_),
			topo_sort_marked_(rhs.topo_sort_marked_),
			topo_nr_(rhs.topo_nr_),
			subtree_finished_(rhs.subtree_finished_)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
		setZValue(42);	
	}
	
	TOPPASVertex::~TOPPASVertex()
	{
		
	}
	
	TOPPASVertex& TOPPASVertex::operator= (const TOPPASVertex& rhs)
	{
		in_edges_ = rhs.in_edges_;
		out_edges_ = rhs.out_edges_;
		edge_being_created_ = rhs.edge_being_created_;
		pen_color_ = rhs.pen_color_;
		brush_color_ = rhs.brush_color_;
		dfs_color_ = rhs.dfs_color_;
		dfs_parent_ = rhs.dfs_parent_;
		id_ = rhs.id_;
		topo_sort_marked_ = rhs.topo_sort_marked_;
		topo_nr_ = rhs.topo_nr_;
		subtree_finished_ = rhs.subtree_finished_;
		
		return *this;
	}
	
	void TOPPASVertex::mousePressEvent(QGraphicsSceneMouseEvent* e)
	{
		if (!(e->modifiers() & Qt::ControlModifier))
		{
			emit clicked();
		}
	}
	
	void TOPPASVertex::mouseReleaseEvent(QGraphicsSceneMouseEvent* e)
	{
		if (edge_being_created_)
		{
			emit finishHoveringEdge();
			edge_being_created_ = false;
		}
		else if (e->modifiers() & Qt::ControlModifier)
		{
			QGraphicsItem::mouseReleaseEvent(e);
		}
		else
		{
			emit released();
			// resize scene rect in case item has been moved outside
			const QRectF& scene_rect = scene()->sceneRect();
			const QRectF& items_bounding = scene()->itemsBoundingRect();
			scene()->setSceneRect(scene_rect.united(items_bounding));
		}
	}
	
	void TOPPASVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e)
	{
		e->ignore();
	}
	
	void TOPPASVertex::contextMenuEvent(QGraphicsSceneContextMenuEvent* e)
	{
		e->ignore();
	}
	
	void TOPPASVertex::mouseMoveEvent(QGraphicsSceneMouseEvent* e)
	{
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		
		if (isSelected())
		{
			QPointF delta = e->pos() - e->lastPos();
			emit itemDragged(delta.x(), delta.y());
		}
		else
		{
			ts->setActionMode(TOPPASScene::AM_NEW_EDGE);
			moveNewEdgeTo_(e->pos());
		}
	}
	
	void TOPPASVertex::moveNewEdgeTo_(const QPointF& pos)
	{	
		if (!edge_being_created_)
		{
			emit newHoveringEdge(mapToScene(pos));
			edge_being_created_ = true;
		}
		
		emit hoveringEdgePosChanged(mapToScene(pos));
	}
	
	TOPPASVertex::EdgeIterator TOPPASVertex::outEdgesBegin()
	{
		return out_edges_.begin();
	}
	
	TOPPASVertex::EdgeIterator TOPPASVertex::outEdgesEnd()
	{
		return out_edges_.end();
	}
	
	TOPPASVertex::EdgeIterator TOPPASVertex::inEdgesBegin()
	{
		return in_edges_.begin();
	}
	
	TOPPASVertex::EdgeIterator TOPPASVertex::inEdgesEnd()
	{
		return in_edges_.end();
	}

	void TOPPASVertex::addInEdge(TOPPASEdge* edge)
	{
		in_edges_.push_back(edge);
	}
	
	void TOPPASVertex::addOutEdge(TOPPASEdge* edge)
	{
		out_edges_.push_back(edge);
	}
	
	void TOPPASVertex::removeInEdge(TOPPASEdge* edge)
	{
		in_edges_.removeAll(edge);
	}
	
	void TOPPASVertex::removeOutEdge(TOPPASEdge* edge)
	{
		out_edges_.removeAll(edge);
	}
	
	TOPPASVertex::DFS_COLOR TOPPASVertex::getDFSColor()
	{
		return dfs_color_;
	}
	
	void TOPPASVertex::setDFSColor(DFS_COLOR color)
	{
		dfs_color_ = color;
	}
	
	TOPPASVertex* TOPPASVertex::getDFSParent()
	{
		return dfs_parent_;
	}
	
	void TOPPASVertex::setDFSParent(TOPPASVertex* parent)
	{
		dfs_parent_ = parent;
	}
	
	Size TOPPASVertex::incomingEdgesCount()
	{
		return in_edges_.size();
	}
	
	Size TOPPASVertex::outgoingEdgesCount()
	{
		return out_edges_.size();
	}
	
	UInt TOPPASVertex::getID()
	{
		return id_;
	}
	
	void TOPPASVertex::setID(UInt id)
	{
		id_ = id;
	}
	
	void TOPPASVertex::inEdgeHasChanged()
	{
		// (overridden behavior in output and tool vertices)
		
		qobject_cast<TOPPASScene*>(scene())->setChanged(true);
		emit somethingHasChanged();
	}
	
	bool TOPPASVertex::isTopoSortMarked()
	{
		return topo_sort_marked_;
	}
	
	void TOPPASVertex::setTopoSortMarked(bool b)
	{
		topo_sort_marked_ = b;
	}
	
	UInt TOPPASVertex::getTopoNr()
	{
		return topo_nr_;
	}
	
	void TOPPASVertex::setTopoNr(UInt nr)
	{
		// (overridden in tool and output vertices)
		topo_nr_ = nr;
	}
	
	String TOPPASVertex::get3CharsNumber_(UInt number)
	{
		String num_str(number);
		int diff = 3 - (int)(num_str.size());
		if (diff <= 0)
		{
			return num_str;
		}
		else
		{
			String res;
			for (int i = 0; i < diff; ++i)
			{
				res += "0";
			}
			res += num_str;
			return res;
		}
	}
	
	bool TOPPASVertex::removeDirRecursively_(const QString& dir_name)
	{
		bool fail = false;
		
		QDir dir(dir_name);
		QStringList files = dir.entryList(QDir::Files | QDir::NoDotAndDotDot);
		foreach (const QString& file_name, files)
		{
			if (!dir.remove(file_name))
			{
				std::cerr << "Could not remove file " << String(file_name) << "!" << std::endl;
				fail = true;
			}
		}
		QStringList contained_dirs = dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
		foreach (const QString& contained_dir, contained_dirs)
		{
			if (!removeDirRecursively_(dir_name+QDir::separator()+contained_dir))
			{
				fail = true;
			}
		}
		
		QDir parent_dir(dir_name);
		if (parent_dir.cdUp())
		{
			if (!parent_dir.rmdir(dir_name))
			{
				std::cerr << "Could not remove directory " << String(dir.dirName()) << "!" << std::endl;
				fail = true;
			}
		}
		
		return !fail;
	}
	
	void TOPPASVertex::checkIfSubtreeFinished()
	{
		__DEBUG_BEGIN_METHOD__
		
		// check if entire subtree below this node is finished now, else return
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* target = (*it)->getTargetVertex();
			if (!target->isSubtreeFinished())
			{
				debugOut_(String("Child ")+target->getTopoNr()+": subtree NOT finished. Returning.");
				__DEBUG_END_METHOD__
				return;
			}
		}
		
		debugOut_("Subtrees of all children finished! Notifying parents...");
		subtree_finished_ = true;
    // tell parents
    for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASVertex* source = (*it)->getSourceVertex();
			debugOut_(String("Notifying parent ")+source->getTopoNr());
			source->checkIfSubtreeFinished();
		}
		
		__DEBUG_END_METHOD__
	}
	
	bool TOPPASVertex::isSubtreeFinished()
	{
		return subtree_finished_;
	}
	
	void TOPPASVertex::setSubtreeFinished(bool b)
	{
		subtree_finished_ = b;
	}
	
	void TOPPASVertex::reset(bool /*reset_all_files*/)
	{
		__DEBUG_BEGIN_METHOD__
		
		subtree_finished_ = false;
		sc_files_per_round_ = 0;
		sc_files_total_ = 0;
		sc_list_length_checked_ = false;
		update(boundingRect());
		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASVertex::resetSubtree(bool including_this_node)
	{
		__DEBUG_BEGIN_METHOD__
		
		if (including_this_node)
		{
			reset(false);
		}
		
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			tv->resetSubtree(true);
		}
		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASVertex::checkListLengths(QStringList& unequal_per_round, QStringList& unequal_over_entire_run, bool merger, bool round_based)
	{
		__DEBUG_BEGIN_METHOD__
		
		//all parents checked?
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			if (!((*it)->getSourceVertex()->sc_list_length_checked_))
			{
				__DEBUG_END_METHOD__
				return;
			}
		}
		
		// do all input lists have equal length?
		int per_round = (*inEdgesBegin())->getSourceVertex()->sc_files_per_round_;
		int total = (*inEdgesBegin())->getSourceVertex()->sc_files_total_;		
		
		if (!merger)
		{
			for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
			{
				if ((*it)->getSourceVertex()->sc_files_per_round_ != per_round)
				{
					unequal_per_round.push_back(QString::number(topo_nr_));
					break;
				}
			}
			
			sc_files_per_round_ = (*inEdgesBegin())->getSourceVertex()->sc_files_per_round_;
			sc_files_total_ = (*inEdgesBegin())->getSourceVertex()->sc_files_total_;
		}
		else // merger
		{
			if (round_based)
			{
				for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
				{
					if ((*it)->getSourceVertex()->sc_files_total_ != total)
					{
						unequal_over_entire_run.push_back(QString::number(topo_nr_));
						break;
					}
				}
				
				// number of in edges determines list length per merging round
				sc_files_per_round_ = inEdgesEnd() - inEdgesBegin();
				// assume all sc_files_total_'s of inputs are equal already
				sc_files_total_ = sc_files_per_round_ * total;
			}
			else
			{
				// no check needed, waiting merger does not care about unequal input list lengths
				sc_files_per_round_ = 0;
				for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
				{
					sc_files_per_round_ += (*it)->getSourceVertex()->sc_files_total_;
				}
				sc_files_total_ = sc_files_per_round_;
			}
		}
		
		sc_list_length_checked_ = true;
		
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			tv->checkListLengths(unequal_per_round, unequal_over_entire_run);
		}
		
		__DEBUG_END_METHOD__
	}
	
}
