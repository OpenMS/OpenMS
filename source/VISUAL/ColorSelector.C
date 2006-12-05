// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/ColorSelector.h>

//qt includes
#include <qpainter.h>
#include <qcolordialog.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

	ColorSelector::ColorSelector( QWidget * parent, const char * name) : QWidget(parent,name),color_(255,255,255)
	{
		setMinimumSize(12,12);
	}
	
	ColorSelector::~ColorSelector()
	{
		
	}
	
	void ColorSelector::paintEvent(QPaintEvent* /*e*/)
	{
		SignedInt size=QMIN(width(),height());
		QPainter painter(this);
		painter.setPen(QColor(0,0,0));
		painter.drawRect(0,0,size,size);
		painter.setPen(QColor(255,255,255));
		painter.drawRect(1,1,size-2,size-2);	
	
		painter.fillRect(2,2,size-4,size-4,color_);
	}
		
	void ColorSelector::mousePressEvent(QMouseEvent* e)
	{
		if ( e->button() != LeftButton ) 
		{
			e->ignore();
			return;
	  } 	
		color_ = QColorDialog::getColor(color_,this, "Color dialog");
		repaint();
	}	
	
	const QColor& ColorSelector::getColor()
	{
		return color_;
	}
	
	void ColorSelector::setColor(const QColor& col)
	{
		color_=col;
		repaint();
	}

} //namespace OpenMS
