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
		if (vt == VT_TOOL)
		{
			pen_color_ = Qt::black;
			brush_color_ = QColor(250,200,0);
		}
		else if (vt == VT_SOURCE)
		{
			pen_color_ = Qt::black;
			brush_color_ = Qt::lightGray;
		}
		else if (vt == VT_TARGET)
		{
			pen_color_ = Qt::black;
			brush_color_ = Qt::lightGray;
		}
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
	
	QPainterPath TOPPASVertex::shape () const
	{
		QPainterPath shape;
		shape.addRoundRect(-70.0, -40.0, 140.0, 80.0, 20, 20);
		
		return shape;
	}
	
	void TOPPASVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		QPainterPath path;
		path.addRoundRect(-70.0, -40.0, 140.0, 80.0, 20, 20);		
		painter->setPen(QPen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
		painter->setBrush(brush_color_);
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
			
			for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
			{
				(*it)->update();
			}
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				(*it)->update();
			}
		}
		else if (action_mode == TOPPASScene::AM_NEW_EDGE)
		{
			if (vertex_type_ == VT_TARGET)
			{
				// cannot create outedge for output files
				return;
			}
			
			if (!edge_being_created_)
			{
				emit newHoveringEdge(e->scenePos());
				edge_being_created_ = true;
			}
			
			emit hoveringEdgePosChanged(e->scenePos());
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
	
		TOPPASVertex::EdgeIterator TOPPASVertex::inEdgesBegin()
	{
		return in_edges_.begin();
	}
	
	TOPPASVertex::EdgeIterator TOPPASVertex::inEdgesEnd()
	{
		return in_edges_.end();
	}

}
