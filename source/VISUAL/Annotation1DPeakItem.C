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

#include <OpenMS/VISUAL/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

namespace OpenMS
{	

	Annotation1DPeakItem::Annotation1DPeakItem(const PointType& position, const PeakIndex& peak, const QString& text, const QPen& pen)
		: Annotation1DItem(text, pen),
			position_(position),
			peak_(peak)
	{
	}
	
	Annotation1DPeakItem::Annotation1DPeakItem(const Annotation1DPeakItem& rhs)
		: Annotation1DItem(rhs)
	{
		position_ = rhs.getPosition();
		peak_ = rhs.getPeak();
	}
	
	Annotation1DPeakItem::~Annotation1DPeakItem()
	{
	}
	
	void Annotation1DPeakItem::draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped)
	{
		//translate mz/intensity to pixel coordinates
		QPoint pos;
		canvas->dataToWidget(position_.getX(), position_.getY(), pos, flipped);
		
		// compute bounding box of text_item on the specified painter
		bounding_box_ = painter.boundingRect(QRectF(pos, pos), Qt::AlignCenter, text_);
		// shift pos, annotation should be next to the peak and not overlap it
		bounding_box_.translate(bounding_box_.width()/2.0 + 10.0, -15.0);
		
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
		painter.drawLine(bounding_box_.bottomLeft(), pos);
	}
	
	void Annotation1DPeakItem::setPosition(const Annotation1DPeakItem::PointType& position)
	{
		position_ = position;
	}
	
 	const Annotation1DPeakItem::PointType& Annotation1DPeakItem::getPosition() const
 	{
 		return position_;
 	}
	
	void Annotation1DPeakItem::setPeak(const PeakIndex& peak)
	{
		peak_ = peak;
	}
	
	const PeakIndex& Annotation1DPeakItem::getPeak() const
	{
		return peak_;
	}
	
}//Namespace




