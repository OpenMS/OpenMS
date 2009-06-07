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

#include <OpenMS/VISUAL/TOPPASInputFileVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>

namespace OpenMS
{
	TOPPASInputFileVertex::TOPPASInputFileVertex()
		:	TOPPASVertex(),
			file_()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileVertex::TOPPASInputFileVertex(const String& name, const String& type)
		: TOPPASVertex(name, type),
			file_()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileVertex::TOPPASInputFileVertex(const TOPPASInputFileVertex& rhs)
		:	TOPPASVertex(rhs),
			file_(rhs.file_)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileVertex::~TOPPASInputFileVertex()
	{
	
	}
	
	TOPPASInputFileVertex& TOPPASInputFileVertex::operator= (const TOPPASInputFileVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);
		
		file_ = rhs.file_;
		
		return *this;
	}
	
	void TOPPASInputFileVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		TOPPASInputFileDialog tifd(file_);
		if (tifd.exec())
		{
			file_ = tifd.getFilename();
		}
	}
	
	const QString& TOPPASInputFileVertex::getFilename()
	{
		return file_;
	}
}
