// -*- mode: C++; tab-width: 2; -*-
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

#ifndef OPENMS_KERNEL_DRICHPEAK_H
#define OPENMS_KERNEL_DRICHPEAK_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <functional>
#include <sstream>

namespace OpenMS
{
	
	/**
		@brief A D-dimensional raw data point or peak mith meta information.
	 
		This datastructure is intended for continuous data or peak data.
		If wou do not need to annotated single peaks with meta data, use DPeak instead.
		
		@ingroup Kernel
	*/
	template <UInt D>
	class DRichPeak
		: public DPeak <D>,
	public MetaInfoInterface
	{
		public:
			
			///@name Type definitions
			//@{
			///Dimensionality
			enum
			{
				DIMENSION = D
			};
			/// Intensity type
			typedef Real IntensityType;
			/// Coordinate type (of the position)
			typedef DoubleReal CoordinateType;
			/// Position type
			typedef DPosition<D> PositionType;
			//@}
			
			///@name Constructors and Destructor
			//@{
			/// Default constructor
			DRichPeak()
				: DPeak<D>(),
					MetaInfoInterface()
			{
				
			}
			
			/// Copy constructor
			inline DRichPeak(DRichPeak const& p)
				: DPeak<D>(p),
					MetaInfoInterface(p)
			{
			}
			
			/// Destructor
			~DRichPeak()
			{
			}
			//@}
			
			/// Assignment operator
			DRichPeak& operator = (const DRichPeak& rhs)
			{
				if (this==&rhs) return *this;
				
				DPeak<D>::operator = (rhs);
				MetaInfoInterface::operator = (rhs);
				
				return *this;
			}
			
			/// Equality operator
			bool operator == (const DRichPeak& rhs) const
			{
				return
					DPeak<D>::operator == (rhs) &&
					MetaInfoInterface::operator == (rhs)
					;
			}
			
			/// Equality operator
			bool operator != (const DRichPeak& rhs) const
			{
				return !(operator == (rhs));
			}
	};
	
	///Print the contents to a stream.
	template <UInt D>
	std::ostream& operator << (std::ostream& os, const DRichPeak<D>& peak)
	{
		os <<(DPeak<D>)peak;
		return os;
	}
	
	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DRICHPEAK_H
