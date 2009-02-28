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
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/TOPPASVertex.h>

// Qt
#include <QtGui/QPainter>

namespace OpenMS
{
	TOPPASVertex::TOPPASVertex(const String& name, VertexType type)
		: QGraphicsItem(),
			name_(name),
			vertex_type_(type)
	{
		// do more
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
		return QRectF(-60,-40,120,80);
	}
	
	void TOPPASVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		painter->drawRoundRect(-60,-40,120,80);
		
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, name_.toQString());
		painter->drawText(-text_boundings.width()/2, text_boundings.height()/4, name_.toQString());
	}

	// here be some methods...
}
