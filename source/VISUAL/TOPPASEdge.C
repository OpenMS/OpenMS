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

#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>

namespace OpenMS
{	
	
	TOPPASEdge::TOPPASEdge()
		:	QObject(),
			QGraphicsItem(),
			from_(0),
			to_(0),
			hover_pos_(),
			color_()
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
	}
	
	TOPPASEdge::TOPPASEdge(TOPPASVertex* from, const QPointF& hover_pos)
		:	QObject(),
			QGraphicsItem(),
			from_(from),
			to_(0),
			hover_pos_(hover_pos),
			color_()
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
	}
	
	TOPPASEdge::TOPPASEdge(const TOPPASEdge& rhs)
		:	QObject(),
			QGraphicsItem(),
			from_(rhs.from_),
			to_(rhs.to_),
			hover_pos_(rhs.hover_pos_),
			color_(rhs.color_)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
	}
	
	TOPPASEdge& TOPPASEdge::operator= (const TOPPASEdge& rhs)
	{
		from_ = rhs.from_;
		to_ = rhs.to_;
		hover_pos_ = rhs.hover_pos_;
		color_ = rhs.color_;
		
		setFlag(QGraphicsItem::ItemIsSelectable, true);
		
		return *this;
	}
	
	TOPPASEdge::~TOPPASEdge()
	{
		if (from_)
		{
			from_->removeOutEdge(this);
		}
		if (to_)
		{
			to_->removeInEdge(this);
		}
	}
	
	QRectF TOPPASEdge::boundingRect() const
	{
		qreal min_x = startPos().x() < endPos().x() ? startPos().x() : endPos().x();
		qreal min_y = startPos().y() < endPos().y() ? startPos().y() : endPos().y();
		qreal max_x = startPos().x() > endPos().x() ? startPos().x() : endPos().x();
		qreal max_y = startPos().y() > endPos().y() ? startPos().y() : endPos().y();
		
		return QRectF(QPointF(min_x-11.0,min_y-11.0), QPointF(max_x+11.0,max_y+11.0));
	}
	
	QPainterPath TOPPASEdge::shape () const
	{
		// this is not quite correct..
		QPainterPath shape;
		shape.setFillRule(Qt::WindingFill);
		shape.moveTo(startPos());
		shape.lineTo(startPos() - QPointF(10,10));
		shape.lineTo(endPos() - QPointF(10,10));
		shape.lineTo(endPos() + QPointF(10,10));
		shape.lineTo(startPos() + QPointF(10,10));
		shape.closeSubpath();
		shape.addEllipse(endPos().x() - 10, endPos().y() - 10, 20, 20);
		
		return shape;
	}
	
	void TOPPASEdge::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		painter->setBrush(Qt::white);

		QPen pen(color_);
		if (isSelected())
		{
			pen.setWidth(2);
		}
		painter->setPen(pen);
		
		painter->drawLine(startPos(),endPos());
		
		// draw arrow head
		QPointF delta = endPos() - startPos();
		qreal angle;
		if (delta.x() == 0.0)
		{
			angle = endPos().y() > startPos().y() ? 90.0 : 270.0;
		}
		else
		{
			angle = delta.y() / delta.x();
			angle = std::atan(angle);
			angle = (angle / 3.14159265) * 180.0;
			if (delta.x() < 0.0)
			{
				angle += 180;
			}
		}
		
		painter->save();
		painter->translate(endPos());
		painter->rotate(angle);
		QPainterPath path;
		path.moveTo(QPointF(0,0));
		path.lineTo(QPointF(-10,4));
		path.lineTo(QPointF(-10,-4));
		path.closeSubpath();
		painter->drawPath(path);
		painter->restore();
	}
	
	QPointF TOPPASEdge::startPos() const
	{
		if (from_)
		{
			return mapFromScene(from_->scenePos());
		}
		else
		{
			return QPointF();
		}
	}
	
	QPointF TOPPASEdge::endPos() const
	{
		QPointF position;
		
		if (!to_)
		{
			// we do not have a target vertex yet
			position = hover_pos_;
		}
		else
		{
			// we have a target node --> line should end at its border
			
			QList<QPointF> point_list;
			
			QPointF target_pos = mapFromScene(to_->scenePos());
			QRectF target_boundings = mapFromItem(to_, to_->shape()).boundingRect();
			QPointF delta = target_pos - startPos();
			qreal slope;
			if (delta.x() == 0)
			{
				slope = std::numeric_limits<double>::infinity();
			}
			else
			{
				slope = delta.y() / delta.x();
			}
			
			qreal y_1 = startPos().y() + slope * (target_boundings.left() - startPos().x());
			qreal y_2 = startPos().y() + slope * (target_boundings.right() - startPos().x());
			
			slope = 1.0 / slope;
			
			qreal x_3 = startPos().x() + slope * (target_boundings.top() - startPos().y());
			qreal x_4 = startPos().x() + slope * (target_boundings.bottom() - startPos().y());
			
			if (y_1 <= target_boundings.bottom() && y_1 >= target_boundings.top())
			{
				point_list.push_back(QPointF(target_boundings.left(), y_1));
			}
			if (y_2 <= target_boundings.bottom() && y_2 >= target_boundings.top())
			{
				point_list.push_back(QPointF(target_boundings.right(), y_2));
			}
			if (x_3 <= target_boundings.right() && x_3 >= target_boundings.left())
			{
				point_list.push_back(QPointF(x_3, target_boundings.top()));
			}
			if (x_4 <= target_boundings.right() && x_4 >= target_boundings.left())
			{
				point_list.push_back(QPointF(x_4, target_boundings.bottom()));
			}
			
			position = nearestPoint_(startPos(), point_list);
		}
		
		return position;
	}
	
	void TOPPASEdge::setHoverPos(const QPointF& pos)
	{
		prepareResize();
		hover_pos_ = pos;
		update();
	}
	
	void TOPPASEdge::setTargetVertex(TOPPASVertex* tv)
	{
		to_ = tv;
	}
	
	void TOPPASEdge::setSourceVertex(TOPPASVertex* tv)
	{
		from_ = tv;
	}
	
	TOPPASVertex* TOPPASEdge::getSourceVertex()
	{
		return from_;
	}
	
	TOPPASVertex* TOPPASEdge::getTargetVertex()
	{
		return to_;
	}
	
	void TOPPASEdge::prepareResize()
	{
		prepareGeometryChange();
	}
	
	QPointF TOPPASEdge::nearestPoint_(const QPointF& origin, const QList<QPointF>& list) const
	{
		if (list.empty())
		{
			return QPointF();
		}
		QPointF nearest = list.first();
		qreal min_distance = std::numeric_limits<double>::max();
		
		for (QList<QPointF>::const_iterator it = list.begin(); it != list.end(); ++it)
		{
			qreal sqr_distance = (it->x() - origin.x()) * (it->x() - origin.x()) +
												(it->y() - origin.y()) * (it->y() - origin.y());
			if (sqr_distance < min_distance)
			{
				min_distance = sqr_distance;
				nearest = *it;
			}
		}
		
		return nearest;
	}
	
	void TOPPASEdge::setColor(const QColor& color)
	{
		color_ = color;
	}


} //namespace


