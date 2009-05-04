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

#include <OpenMS/VISUAL/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

#include <QtGui/QPainter>
#include <QtCore/QPoint>

namespace OpenMS
{	

	Annotation1DPeakItem::Annotation1DPeakItem(const PointType& position, const QString& text)
		: Annotation1DItem(text),
			position_(position)
	{
	}
	
	Annotation1DPeakItem::Annotation1DPeakItem(const Annotation1DPeakItem& rhs)
		: Annotation1DItem(rhs)
	{
		position_ = rhs.getPosition();
	}
	
	Annotation1DPeakItem::~Annotation1DPeakItem()
	{
	}
	
	void Annotation1DPeakItem::draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped)
	{
		//translate mz/intensity to pixel coordinates
		QPoint pos;
		canvas->dataToWidget(position_.getX(), position_.getY(), pos, flipped, true);
		
		// compute bounding box of text_item on the specified painter
		bounding_box_ = painter.boundingRect(QRectF(pos, pos), Qt::AlignCenter, text_);
		
		if (canvas->isMzToXAxis())
		{
			// shift pos - annotation should be over peak or, if not possible, next to it
			DoubleReal vertical_shift = bounding_box_.height()/2 + 5;
			if (!flipped)
			{
				vertical_shift *= -1;
			}
			bounding_box_.translate(0.0, vertical_shift);
			if (flipped && bounding_box_.bottom() > canvas->height())
			{
				bounding_box_.moveBottom(canvas->height());
				bounding_box_.moveLeft(pos.x() + 5.0);
			}
			else if (!flipped && bounding_box_.top() < 0.0)
			{
				bounding_box_.moveTop(0.0);
				bounding_box_.moveLeft(pos.x() + 5.0);
			}
		}
		else
		{
			// annotation should be next to the peak (to its right)
			DoubleReal horizontal_shift = bounding_box_.width()/2 + 5;
			bounding_box_.translate(horizontal_shift, 0.0);
			if (bounding_box_.right() > canvas->width())
			{
				bounding_box_.moveRight(canvas->width());
			}
		}
		
		painter.drawText(bounding_box_, Qt::AlignCenter, text_);
		if (selected_)
		{
			drawBoundingBox_(painter);
		}
	}
	
	void Annotation1DPeakItem::move(const PointType& /*delta*/)
	{
		// do nothing, peak annotations cannot be moved
	}
	
	void Annotation1DPeakItem::setPosition(const Annotation1DPeakItem::PointType& position)
	{
		position_ = position;
	}
	
 	const Annotation1DPeakItem::PointType& Annotation1DPeakItem::getPosition() const
 	{
 		return position_;
 	}
 	
 	void Annotation1DPeakItem::ensureWithinDataRange(Spectrum1DCanvas* const /*canvas*/)
	{
		// peak items are never outside data range (cannot be moved)
	}
	
}//Namespace




