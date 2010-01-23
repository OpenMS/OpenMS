// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_VISUAL_ANNOTATION1DPEAKITEM_H
#define OPENMS_VISUAL_ANNOTATION1DPEAKITEM_H

#include <OpenMS/VISUAL/Annotation1DItem.h>

namespace OpenMS
{
	/** @brief A peak annotation item
			@see Annotation1DItem
	*/
	class Annotation1DPeakItem
		: public Annotation1DItem
	{
		
		public:
			
			/// Constructor
			Annotation1DPeakItem(const PointType& position, const QString& text);
			/// Copy constructor
			Annotation1DPeakItem(const Annotation1DPeakItem& rhs);
			/// Destructor
			virtual ~Annotation1DPeakItem();
			// Docu in base class
			virtual void ensureWithinDataRange(Spectrum1DCanvas* const /*canvas*/);
			// Docu in base class
			virtual void draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped = false);
			// Docu in base class
			virtual void move(const PointType& /*delta*/);
			/// Sets the position of the peak (in MZ / intensity coordinates)
			void setPosition(const PointType& position);
			/// Returns the position of the peak (in MZ / intensity coordinates)
 			const PointType& getPosition() const;
			
		protected:
		
			/// The position of the peak (in MZ / intensity coordinates)
			PointType position_;
			
	};
} // namespace OpenMS

#endif
