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

#ifndef OPENMS_KERNEL_DPEAK_H
#define OPENMS_KERNEL_DPEAK_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/DRawDataPoint.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <functional>
#include <sstream>

namespace OpenMS
{

	/**	
		@brief A D-dimensional peak.
		
		This datastructure is indended for picked peaks, which have no information
		from peak picking annotated. If you want to handle peaks that have such
		information, use DPickedPeak.
		<BR>
		The intensity of a peak is defined as the maximum of the model fitted to the raw data during peak picking
		 i.e. aproximately the height of the highest raw data point.
	
		@ingroup Kernel
	*/
	template <UnsignedInt D>
	class DPeak 
		: public DRawDataPoint <D>, 
			public MetaInfoInterface
	{
		public:
		
		/** @name Type definitions
		*/
		//@{	
		enum { DIMENSION = D };
		typedef DPosition<D> PositionType;
		typedef DoubleReal CoordinateType;
		typedef DoubleReal IntensityType;
		//@}

		/** @name Constructors and Destructor
		*/
		//@{
		/// Default constructor
		DPeak() 
			: DRawDataPoint<D>(),
				MetaInfoInterface()
		{
			
		}
		
		/// Copy constructor
		inline DPeak(DPeak const& p) 
			: DRawDataPoint<D>(p),
				MetaInfoInterface(p)
		{
		}

		/// Destructor
		~DPeak() 
		{
		}
		//@}	

		/// Assignment operator
		DPeak& operator = (const DPeak& rhs)
		{
			if (this==&rhs) return *this;
			
			DRawDataPoint<D>::operator = (rhs);
			MetaInfoInterface::operator = (rhs);
			
			return *this;
		}

		/// Equality operator
		bool operator == (const DPeak& rhs) const
		{
			return 
				DRawDataPoint<D>::operator == (rhs) &&
				MetaInfoInterface::operator == (rhs)
			;
		}

		/// Equality operator
		bool operator != (const DPeak& rhs) const
		{
			return !(operator == (rhs));
		}
	};

	///Print the contents to a stream.
	template <UnsignedInt D>
	std::ostream& operator << (std::ostream& os, const DPeak<D>& peak)
	{
		os <<(DRawDataPoint<D>)peak;
		return os;
	}


} // namespace OpenMS

#endif // OPENMS_KERNEL_DPEAK_H
