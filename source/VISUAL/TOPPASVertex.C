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

// Qt
#include <QtGui/QPainter>
#include <QtGui/QGraphicsSceneMouseEvent>

namespace OpenMS
{
	TOPPASVertex::TOPPASVertex(const String& name, const String& type, VertexType vt)
		: QGraphicsItem(),
			name_(name),
			type_(type),
			vertex_type_(vt)
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
		painter->drawRoundRect(-70,-40,140,80);
		if (type_ == "")
		{
			QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, name_.toQString());
			painter->drawText(-text_boundings.width()/2, text_boundings.height()/4, name_.toQString());
		}
		else
		{
			QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, name_.toQString());
			painter->drawText(-text_boundings.width()/2, -text_boundings.height()/4, name_.toQString());
			text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, type_.toQString());
			painter->drawText(-text_boundings.width()/2, +text_boundings.height()/1.5, type_.toQString());
		}
	}
	
	void TOPPASVertex::mousePressEvent(QGraphicsSceneMouseEvent* e)
	{
		last_mouse_pos_ = e->pos();
	}
	
	void TOPPASVertex::mouseMoveEvent(QGraphicsSceneMouseEvent* e)
	{
		QPointF delta = e->pos() - last_mouse_pos_;
		moveBy(delta.x(), delta.y());
		last_mouse_pos_ = e->pos();
	}

}
