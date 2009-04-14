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

#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>

namespace OpenMS
{
	
	TOPPASScene::TOPPASScene()
		:	QGraphicsScene(),
			action_mode_(AM_MOVE),
			vertices_(),
			edges_(),
			hover_edge_(0)
	{
	}
	
	TOPPASScene::~TOPPASScene()
	{
		// Qt should clean up for us..
	}
	
	void TOPPASScene::setActionMode(ActionMode mode)
	{
		action_mode_ = mode;
	}
	
	TOPPASScene::ActionMode TOPPASScene::getActionMode()
	{
		return action_mode_;
	}
	
	TOPPASScene::VertexIterator TOPPASScene::verticesBegin()
	{
		return vertices_.begin();
	}
	
	TOPPASScene::VertexIterator TOPPASScene::verticesEnd()
	{
		return vertices_.end();
	}
	
	TOPPASScene::EdgeIterator TOPPASScene::edgesBegin()
	{
		return edges_.begin();
	}
	
	TOPPASScene::EdgeIterator TOPPASScene::edgesEnd()
	{
		return edges_.end();
	}
	
	void TOPPASScene::addVertex(TOPPASVertex* tv)
	{
		vertices_.push_back(tv);
		addItem(tv);
	}
	
	void TOPPASScene::addEdge(TOPPASEdge* te)
	{
		edges_.push_back(te);
		addItem(te);
	}
	
	void TOPPASScene::itemClicked()
	{
		TOPPASVertex* sender = qobject_cast<TOPPASVertex*>(QObject::sender());
		if (!sender)
		{
			return;
		}
		
		if (getActionMode() == AM_MOVE)
		{
			std::cout << "AM_MOVE" << std::endl;
		}
		else if (getActionMode() == AM_NEW_EDGE)
		{
			std::cout << "AM_NEW_EDGE" << std::endl;
		}
	}
	
	void TOPPASScene::itemDoubleClicked()
	{
		std::cout << "double click!" << std::endl;
	}
	
	void TOPPASScene::updateHoveringEdgePos(const QPointF& new_pos)
	{
		if (!hover_edge_)
		{
			return;
		}
		
		hover_edge_->setHoverPos(new_pos);
	}
	
	void TOPPASScene::addHoveringEdge(const QPointF& pos)
	{
		TOPPASVertex* sender = dynamic_cast<TOPPASVertex*>(QObject::sender());
		if (!sender)
		{
			return;
		}
		TOPPASEdge* new_edge = new TOPPASEdge(sender, pos);
		hover_edge_ = new_edge;
		addEdge(new_edge);
	}
	
	void TOPPASScene::finishHoveringEdge()
	{
		QList<QGraphicsItem*> target_list = items(hover_edge_->endPos());
		bool destroy = true;
		
		// if one of the items at this position is a vertex: use it as target of the edge
		for (QList<QGraphicsItem*>::iterator it = target_list.begin(); it != target_list.end(); ++it)
		{
			TOPPASVertex* target = dynamic_cast<TOPPASVertex*>(*it);
			if (target)
			{
				hover_edge_->setTargetVertex(target);
				TOPPASVertex* source = hover_edge_->getSourceVertex();
				source->addOutEdge(hover_edge_);
				target->addInEdge(hover_edge_);
				
				hover_edge_ = 0;
				destroy = false;
				break;
			}
		}
		
		if (destroy && hover_edge_ != 0)
		{
			removeItem(hover_edge_);
			delete hover_edge_;
			hover_edge_ = 0;
		}
		update();
	}
	
	void TOPPASScene::removeSelected()
	{
		QList<TOPPASVertex*> vertices_to_be_removed;
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			if ((*it)->isSelected())
			{
				// also select all in and out edges (will be deleted below)
				for (TOPPASVertex::EdgeIterator e_it = (*it)->inEdgesBegin(); e_it != (*it)->inEdgesEnd(); ++e_it)
				{
					(*e_it)->setSelected(true);
				}
				for (TOPPASVertex::EdgeIterator e_it = (*it)->outEdgesBegin(); e_it != (*it)->outEdgesEnd(); ++e_it)
				{
					(*e_it)->setSelected(true);
				}
				vertices_to_be_removed.push_back(*it);
				// remove from scene
				removeItem(*it);
			}
		}
		
		QList<TOPPASEdge*> edges_to_be_removed;
		for (EdgeIterator it = edgesBegin(); it != edgesEnd(); ++it)
		{
			if ((*it)->isSelected())
			{
				edges_to_be_removed.push_back(*it);
				removeItem(*it);
			}
		}
		
		TOPPASEdge* edge;
		foreach (edge, edges_to_be_removed)
		{
			edges_.removeAll(edge);
			delete edge;
		}
		TOPPASVertex* vertex;
		foreach (vertex, vertices_to_be_removed)
		{
			vertices_.removeAll(vertex);
			delete vertex;
		}
	}
	
} //namespace OpenMS

