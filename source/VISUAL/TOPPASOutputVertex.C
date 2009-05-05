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

#include <OpenMS/VISUAL/TOPPASOutputVertex.h>

namespace OpenMS
{
	TOPPASOutputVertex::TOPPASOutputVertex()
		:	TOPPASVertex()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputVertex::TOPPASOutputVertex(const String& name, const String& type)
		: TOPPASVertex(name, type)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputVertex::TOPPASOutputVertex(const TOPPASOutputVertex& rhs)
		:	TOPPASVertex(rhs)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputVertex::~TOPPASOutputVertex()
	{
	
	}
	
	TOPPASOutputVertex& TOPPASOutputVertex::operator= (const TOPPASOutputVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);		
		
		return *this;
	}
	
	void TOPPASOutputVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		// ...
	}
}

