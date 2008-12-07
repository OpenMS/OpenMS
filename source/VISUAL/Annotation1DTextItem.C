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

#include <OpenMS/VISUAL/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

namespace OpenMS
{	

	Annotation1DTextItem::Annotation1DTextItem(const PointType& position, const QString& text, const QPen& pen)
		: Annotation1DItem(text, pen),
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
		canvas->dataToWidget(position_.getX(), position_.getY(), pos, flipped);
		
		// compute bounding box of text_item on the specified painter
		bounding_box_ = painter.boundingRect(QRectF(pos, pos), Qt::AlignCenter, text_);
		
		if (selected_)
		{
			painter.setPen(selected_pen_);
			drawBoundingBox_(painter);
		}
		else
		{
			painter.setPen(pen_);
		}
		
		painter.drawText(bounding_box_, Qt::AlignCenter, text_);
	}
	
	void Annotation1DTextItem::setPosition(const Annotation1DTextItem::PointType& position)
	{
		position_ = position;
	}
	
 	const Annotation1DTextItem::PointType& Annotation1DTextItem::getPosition() const
 	{
 		return position_;
 	}
	
}//Namespace




