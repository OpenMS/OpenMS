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
// $Id: PeakShape.h,v 1.11 2006/04/11 15:29:39 elange Exp $
// $Author: elange $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKSHAPE_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKSHAPE_H

#include <math.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/DPickedPeak.h>

namespace OpenMS 
{

		/** @brief This class is a internal representation (used by the DPeakPickerCWT) of a peak shape.

				It defines an asymmetric lorentzian and asymmetric hyperbolic squared secan function. 

				@todo write test (Eva)
		*/
		//@{
	
	class PeakShape
		{
		public:

			PeakShape() 
				: height(0),
					mz_position(0),
					rt_position(0),
					left_width(0),
					right_width(0),
					area(0),
					r_value(0),
					signal_to_noise(0.),
					type(PeakShapeType::UNDEFINED)
			{}
			///
			PeakShape(double height_, 
								double mz_position_, 
								double rt_position_,
								double left_width_, 
								double right_width_,
								double area_,
							  PeakShapeType::Enum type_);
			///	
			PeakShape(const PeakShape& peakshape);
			///
			virtual ~PeakShape(){}
			///
			PeakShape& operator = (const PeakShape& peakshape);
			
			/// compute the real value of the fitted Peak at position x
			double operator() (const double x) const;
			///
			double getSymmetricMeasure() const;
			///
			double getFWHM() const;
			///
			double height;
			/// 
			double mz_position;
			///
			double rt_position;
			///
			double left_width;
			///
			double right_width;
			///
			double area;
			///
			double r_value;
			///
			double signal_to_noise;
			///
			PeakShapeType::Enum type;

			class PositionLess
			{
			public:
				
				PositionLess(Index i) : dimension_(i) {}
				PositionLess() : dimension_(-1) {}
				~PositionLess() {}
				
				bool operator () (const PeakShape& a, const PeakShape& b)
				{
					if (dimension_==2)
						{
							return ((a.rt_position < b.rt_position) || (!(a.rt_position > b.rt_position) && (a.mz_position < b.mz_position)));
						}
					else
						{
							return (a.mz_position < b.mz_position);
						}
				}
				
			protected:
				Index dimension_;
			};

	};
	
	

} // namespace OpenMS

#endif 
