// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/VISUAL/Annotation1DItem.h>

namespace OpenMS
{

	Annotation1DItem::Annotation1DItem(const QString& text, const QPen& pen)
		: bounding_box_(),
		  selected_(true),
		  text_(text)
	{
		setPen(pen);
	}
	
	Annotation1DItem::Annotation1DItem(const Annotation1DItem& rhs)
	{
		bounding_box_ = rhs.boundingBox();
		selected_ = rhs.isSelected();
		text_ = rhs.getText();
		pen_ = rhs.getPen();
		selected_pen_ = rhs.getSelectedPen();
	}
	
	Annotation1DItem::~Annotation1DItem()
	{
	}
	
	void Annotation1DItem::drawBoundingBox_(QPainter& painter)
	{
		// draw additional filled rectangles to highlight bounding box of selected distance_item
		painter.fillRect(bounding_box_.topLeft().x()-3, bounding_box_.topLeft().y()-3, 3, 3, selected_pen_.color());
		painter.fillRect(bounding_box_.topRight().x(), bounding_box_.topRight().y()-3, 3, 3, selected_pen_.color());
		painter.fillRect(bounding_box_.bottomRight().x(), bounding_box_.bottomRight().y(), 3, 3, selected_pen_.color());
		painter.fillRect(bounding_box_.bottomLeft().x()-3, bounding_box_.bottomLeft().y(), 3, 3, selected_pen_.color());
	}
	
	const QRectF& Annotation1DItem::boundingBox() const
	{
		return bounding_box_;
	}
	
	void Annotation1DItem::setSelected(bool selected)
	{
		selected_ = selected;
	}
	
	bool Annotation1DItem::isSelected() const
	{
		return selected_;
	}
	
	void Annotation1DItem::setText(const QString& text)
	{
		text_ = text;
	}
	
	const QString& Annotation1DItem::getText() const
	{
		return text_;
	}
	
	void Annotation1DItem::setPen(const QPen& pen)
	{
		pen_ = pen;
		selected_pen_ = pen;
		
		//make selected items a little brighter
		int sel_red = pen.color().red() + 50;
		int sel_green = pen.color().green() + 50;
		int sel_blue = pen.color().blue() + 50;
		//check if rgb out of bounds
		sel_red = sel_red > 255 ? 255 : sel_red;
		sel_green = sel_green > 255 ? 255 : sel_green;
		sel_blue = sel_blue > 255 ? 255 : sel_blue;
		
		selected_pen_.setColor(QColor(sel_red, sel_green, sel_blue));
	}
	
	const QPen& Annotation1DItem::getPen() const
	{
		return pen_;
	}
	
	void Annotation1DItem::setSelectedPen(const QPen& pen)
	{
		selected_pen_ = pen;
	}
	
	const QPen& Annotation1DItem::getSelectedPen() const
	{
		return selected_pen_;
	}
  
}//Namespace
