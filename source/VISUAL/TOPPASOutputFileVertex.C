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

#include <OpenMS/VISUAL/TOPPASOutputFileVertex.h>

namespace OpenMS
{
	TOPPASOutputFileVertex::TOPPASOutputFileVertex()
		:	TOPPASVertex()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileVertex::TOPPASOutputFileVertex(const TOPPASOutputFileVertex& rhs)
		:	TOPPASVertex(rhs)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileVertex::~TOPPASOutputFileVertex()
	{
	
	}
	
	TOPPASOutputFileVertex& TOPPASOutputFileVertex::operator= (const TOPPASOutputFileVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);		
		
		return *this;
	}
	
	void TOPPASOutputFileVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		// ...
	}
	
	void TOPPASOutputFileVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		QPen pen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin);
		if (isSelected())
		{
			pen.setWidth(2);
			painter->setBrush(brush_color_.darker(130));
		}
		else
		{
			painter->setBrush(brush_color_);
		}
		painter->setPen(pen);
	
		QPainterPath path;
		path.addRoundRect(-70.0, -40.0, 140.0, 80.0, 20, 20);		
 		painter->drawPath(path);
 		
		QString text = "Output file";
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);
	}
	
	QRectF TOPPASOutputFileVertex::boundingRect() const
	{
		return QRectF(-71,-41,142,82);
	}
	
	QPainterPath TOPPASOutputFileVertex::shape () const
	{
		QPainterPath shape;
		shape.addRoundRect(-71.0, -41.0, 142.0, 81.0, 20, 20);
		return shape;
	}
}

