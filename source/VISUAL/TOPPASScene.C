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
		TOPPASVertex* sender = dynamic_cast<TOPPASVertex*>(QObject::sender());
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
		
		hover_edge_->setPos(new_pos);
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
	
} //namespace OpenMS

