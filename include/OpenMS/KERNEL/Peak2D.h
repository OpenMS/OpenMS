// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_KERNEL_PEAK2D_H
#define OPENMS_KERNEL_PEAK2D_H

#include <OpenMS/KERNEL/RawDataPoint2D.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{

	/**	
		@brief A 2-dimensional peak.
		
		This datastructure is indended for picked peaks, which have no information
		from peak picking annotated. If you want to handle peaks that have such
		information, use PickedPeak2D.

		The intensity of a peak is defined as the maximum of the model fitted to the raw data during peak picking
		 i.e. aproximately the height of the highest raw data point.
	
		@ingroup Kernel
	*/
	class Peak2D 
		: public RawDataPoint2D, 
			public MetaInfoInterface
	{
		public:
		
		/// Default constructor
		Peak2D() 
			: RawDataPoint2D(),
				MetaInfoInterface()
		{
			
		}
		
		/// Copy constructor
		inline Peak2D(const Peak2D& p) 
			: RawDataPoint2D(p),
				MetaInfoInterface(p)
		{
		}

		/// Destructor
		~Peak2D() 
		{
		}
    
		/// Assignment operator
		inline Peak2D& operator = (const Peak2D& rhs)
		{
			if (this==&rhs) return *this;
			
			RawDataPoint2D::operator = (rhs);
			MetaInfoInterface::operator = (rhs);
			
			return *this;
		}

		/// Equality operator
		inline bool operator == (const Peak2D& rhs) const
		{
			return 
				RawDataPoint2D::operator == (rhs) &&
				MetaInfoInterface::operator == (rhs)
			;
		}

		/// Equality operator
		inline bool operator != (const Peak2D& rhs) const
		{
			return !(operator == (rhs));
		}
	
	};

} // namespace OpenMS

#endif // OPENMS_KERNEL_DPEAK_H
