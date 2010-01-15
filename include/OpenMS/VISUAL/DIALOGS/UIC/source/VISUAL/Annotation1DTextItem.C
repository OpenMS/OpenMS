// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/VISUAL/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

#include <QtGui/QPainter>
#include <QtCore/QPoint>


namespace OpenMS
{	

	Annotation1DTextItem::Annotation1DTextItem(const PointType& position, const QString& text)
		: Annotation1DItem(text),
			position_(position)
	{
	}
	
	Annotation1DTextItem::Annotation1DTextItem(const Annotation1DTextItem& rhs)
		: Annotation1DItem(rhs)
	{
		position_ = rhs.getPosition();
	}
	
	Annotation1DTextItem::~Annotation1DTextItem()
	{
	}
	
	void Annotation1DTextItem::draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped)
	{
		//translate mz/intensity to pixel coordinates
		QPoint pos;
		canvas->dataToWidget(position_.getX(), position_.getY(), pos, flipped, true);
		
		// compute bounding box of text_item on the specified painter
		bounding_box_ = painter.boundingRect(QRectF(pos, pos), Qt::AlignCenter, text_);
		
		painter.drawText(bounding_box_, Qt::AlignCenter, text_);
		if (selected_)
		{
			drawBoundingBox_(painter);
		}
	}
	
	void Annotation1DTextItem::move(const PointType& delta)
	{
		position_.setX(position_.getX()+delta.getX());
		position_.setY(position_.getY()+delta.getY());
	}
	
	void Annotation1DTextItem::setPosition(const Annotation1DTextItem::PointType& position)
	{
		position_ = position;
	}
	
 	const Annotation1DTextItem::PointType& Annotation1DTextItem::getPosition() const
 	{
 		return position_;
 	}
 	
 	void Annotation1DTextItem::ensureWithinDataRange(Spectrum1DCanvas* const canvas)
	{
		DRange<3> data_range = canvas->getDataRange();
		
		CoordinateType x_pos = position_.getX();
		CoordinateType y_pos = position_.getY() * canvas->getPercentageFactor();
		
		if (x_pos < data_range.min()[0])
		{
			position_.setX(data_range.min()[0]);
		}
		if (x_pos > data_range.max()[0])
		{
			position_.setX(data_range.max()[0]);
		}
		if (y_pos < data_range.min()[1])
		{
			position_.setY(data_range.min()[1] / canvas->getPercentageFactor());
		}
		if (y_pos > data_range.max()[1])
		{
			position_.setY(data_range.max()[1] / canvas->getPercentageFactor());
		}
	}
	
}//Namespace




