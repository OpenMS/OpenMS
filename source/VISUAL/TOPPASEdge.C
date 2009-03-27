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

#include <OpenMS/VISUAL/TOPPASEdge.h>

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>

namespace OpenMS
{	
	
	TOPPASEdge::TOPPASEdge(TOPPASVertex* from, const QPointF& hover_pos)
		:	QGraphicsItem(),
			from_(from),
			to_(0),
			hover_pos_(hover_pos)
	{
		
	}
	
	TOPPASEdge::~TOPPASEdge()
	{
	
	}
	
	QRectF TOPPASEdge::boundingRect() const
	{
		qreal min_x = startPos().x() < endPos().x() ? startPos().x() : endPos().x();
		qreal min_y = startPos().y() < endPos().y() ? startPos().y() : endPos().y();
		qreal max_x = startPos().x() > endPos().x() ? startPos().x() : endPos().x();
		qreal max_y = startPos().y() > endPos().y() ? startPos().y() : endPos().y();
		
		return QRectF(QPointF(min_x,min_y), QPointF(max_x,max_y));
	}
	
	QPainterPath TOPPASEdge::shape () const
	{
		// this should be slightly more precise..
		QPainterPath shape;
		shape.addRect(boundingRect());
		
		return shape;
	}
	
	void TOPPASEdge::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		painter->drawLine(startPos(),endPos());
	}
	
	QPointF TOPPASEdge::startPos() const
	{
		QPointF position = mapFromScene(from_->scenePos());
		
		return position;
	}
	
	QPointF TOPPASEdge::endPos() const
	{
		QPointF position;
		
		if (!to_)
		{
			// we do not have a target vertex yet
			position = mapFromScene(hover_pos_);
		}
		else
		{
			position = mapFromScene(to_->scenePos());
		}
		
		return position;
	}
	
	void TOPPASEdge::setHoverPos(const QPointF& pos)
	{
		hover_pos_ = pos;
	}
	
	void TOPPASEdge::setTargetVertex(TOPPASVertex* tv)
	{
		to_ = tv;
	}
	
} //namespace


