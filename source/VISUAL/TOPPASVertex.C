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

namespace OpenMS
{
	TOPPASVertex::TOPPASVertex()
		:	QObject(),
			QGraphicsItem(),
			name_(),
			type_(),
			in_edges_(),
			out_edges_(),
			edge_being_created_(false),
			pen_color_(),
			brush_color_(),
			dfs_color_(DFS_WHITE),
			dfs_parent_(0)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
		setZValue(42);
	}
	
	TOPPASVertex::TOPPASVertex(const String& name, const String& type)
		: QObject(),
			QGraphicsItem(),
			name_(name),
			type_(type),
			in_edges_(),
			out_edges_(),
			edge_being_created_(false),
			pen_color_(),
			brush_color_(),
			dfs_color_(DFS_WHITE),
			dfs_parent_(0)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
		setZValue(42);
	}
	
	TOPPASVertex::TOPPASVertex(const TOPPASVertex& rhs)
		:	QObject(),
			QGraphicsItem(),
			name_(rhs.name_),
			type_(rhs.type_),
			in_edges_(rhs.in_edges_),
			out_edges_(rhs.out_edges_),
			edge_being_created_(rhs.edge_being_created_),
			pen_color_(rhs.pen_color_),
			brush_color_(rhs.brush_color_),
			dfs_color_(rhs.dfs_color_),
			dfs_parent_(rhs.dfs_parent_)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
		setZValue(42);	
	}
	
	TOPPASVertex::~TOPPASVertex()
	{
		foreach (TOPPASEdge* edge, in_edges_)
		{
			edge->setTargetVertex(0);
		}
		foreach (TOPPASEdge* edge, out_edges_)
		{
			edge->setSourceVertex(0);
		}
	}
	
	TOPPASVertex& TOPPASVertex::operator= (const TOPPASVertex& rhs)
	{
		name_ = rhs.name_;
		type_ = rhs.type_;
		in_edges_ = rhs.in_edges_;
		out_edges_ = rhs.out_edges_;
		edge_being_created_ = rhs.edge_being_created_;
		pen_color_ = rhs.pen_color_;
		brush_color_ = rhs.brush_color_;
		dfs_color_ = rhs.dfs_color_;
		dfs_parent_ = rhs.dfs_parent_;
		
		return *this;
	}
	
	const String& TOPPASVertex::getName()
	{
		return name_;
	}
	
	const String& TOPPASVertex::getType()
	{
		return type_;
	}
	
	QRectF TOPPASVertex::boundingRect() const
	{
		return QRectF(-71,-41,142,82);
	}
	
	QPainterPath TOPPASVertex::shape () const
	{
		QPainterPath shape;
		shape.addRoundRect(-70.0, -40.0, 140.0, 80.0, 20, 20);
		
		return shape;
	}
	
	void TOPPASVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		painter->setPen(QPen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
		if (isSelected())
		{
			painter->setBrush(brush_color_.darker(170));
		}
		else
		{
			painter->setBrush(brush_color_);
		}
	
		QPainterPath path;
		path.addRoundRect(-70.0, -40.0, 140.0, 80.0, 20, 20);		
 		painter->drawPath(path);
 		
		if (type_ == "")
		{
			QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, name_.toQString());
			painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), name_.toQString());
		}
		else
		{
			QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, name_.toQString());
			painter->drawText(-(int)(text_boundings.width()/2.0), -(int)(text_boundings.height()/4.0), name_.toQString());
			text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, type_.toQString());
			painter->drawText(-(int)(text_boundings.width()/2.0), +(int)(text_boundings.height()/1.5), type_.toQString());
		}
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
			setSelected(true);
		}
	}
	
	void TOPPASVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
	
	}
	
	void TOPPASVertex::mouseMoveEvent(QGraphicsSceneMouseEvent* e)
	{
		TOPPASScene* ts = static_cast<TOPPASScene*>(scene());
		
		if (isSelected())
		{
			ts->setActionMode(TOPPASScene::AM_MOVE);
			for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
			{
				(*it)->prepareResize();
			}
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				(*it)->prepareResize();
			}
			
			QPointF delta = e->pos() - e->lastPos();
			moveBy(delta.x(), delta.y());
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
	
}
