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

// OpenMS
#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASScene.h>

// Qt
#include <QtGui/QPainter>
#include <QtGui/QPainterPath>
#include <QtGui/QGraphicsSceneMouseEvent>

namespace OpenMS
{
	TOPPASVertex::TOPPASVertex(const String& name, const String& type, VertexType vt)
		: QGraphicsItem(),
			name_(name),
			type_(type),
			vertex_type_(vt),
			out_edges_(),
			edge_being_created_(false)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
	}
	
	TOPPASVertex::~TOPPASVertex()
	{
	
	}
	
	const String& TOPPASVertex::getName()
	{
		return name_;
	}
	
	QRectF TOPPASVertex::boundingRect() const
	{
		return QRectF(-70,-40,140,80);
	}
	
	void TOPPASVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		QPainterPath path;
		path.addRoundRect(-70.0, -40.0, 140.0, 80.0, 20, 20);		
		painter->setPen(QPen(QColor(0, 0, 0), 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
		painter->setBrush(QColor(250,200,0));
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
		
		if (isSelected())
		{
			// draw selection rectangle
		}
	}
	
	void TOPPASVertex::mousePressEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		emit clicked();
	}
	
	void TOPPASVertex::mouseReleaseEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		if (edge_being_created_)
		{
			// TODO construct new edge if pos of new edge on other node
		}
		
		edge_being_created_ = false;
	}
	
	void TOPPASVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		emit doubleClicked();
	}
	
	void TOPPASVertex::mouseMoveEvent(QGraphicsSceneMouseEvent* e)
	{
		TOPPASScene* ts = static_cast<TOPPASScene*>(scene());
		TOPPASScene::ActionMode action_mode = ts->getActionMode();
		
		if (action_mode == TOPPASScene::AM_MOVE)
		{
			QPointF delta = e->pos() - e->lastPos();
			moveBy(delta.x(), delta.y());
			
			if (!ts->collidingItems(this).empty())
			{
				// we collide with other items --> move back to old position
				moveBy(-delta.x(), -delta.y());
			}
		}
		else if (action_mode == TOPPASScene::AM_NEW_EDGE)
		{
			if (!edge_being_created_)
			{
				emit newHoveringEdge(e->pos());
				edge_being_created_ = true;
			}
			
			emit hoveringEdgePosChanged(e->pos());
		}
	}
	
	TOPPASVertex::EdgeIterator TOPPASVertex::outEdgesBegin()
	{
		return out_edges_.begin();
	}
	
	TOPPASVertex::EdgeIterator TOPPASVertex::outEdgesEnd()
	{
		return out_edges_.end();
	}

}
