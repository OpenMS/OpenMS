// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>

#include <QtGui/QPainter>

namespace OpenMS
{

	Annotation1DItem::Annotation1DItem(const QString& text)
    : bounding_box_(),
		  selected_(true),
		  text_(text)
	{
	}
	
	Annotation1DItem::Annotation1DItem(const Annotation1DItem& rhs)
  {
    bounding_box_ = rhs.boundingBox();
    selected_ = rhs.isSelected();
    text_ = rhs.getText();
	}
	
	Annotation1DItem::~Annotation1DItem()
	{
	}
	
	void Annotation1DItem::drawBoundingBox_(QPainter& painter)
	{
		// draw additional filled rectangles to highlight bounding box of selected distance_item
		painter.fillRect((int)(bounding_box_.topLeft().x())-3, (int)(bounding_box_.topLeft().y())-3, 3, 3, painter.pen().color());
		painter.fillRect((int)(bounding_box_.topRight().x()), (int)(bounding_box_.topRight().y())-3, 3, 3, painter.pen().color());
		painter.fillRect((int)(bounding_box_.bottomRight().x()), (int)(bounding_box_.bottomRight().y()), 3, 3, painter.pen().color());
		painter.fillRect((int)(bounding_box_.bottomLeft().x())-3, (int)(bounding_box_.bottomLeft().y()), 3, 3, painter.pen().color());
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
  
}//Namespace
