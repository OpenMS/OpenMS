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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_KERNEL_KERNELTRAITS_H
#define OPENMS_KERNEL_KERNELTRAITS_H

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

	// forward declarations ...

	struct DoubleKernelTraits;
	struct FloatKernelTraits;

	/**
		 @brief Default traits class for the KERNEL classes.
		 
		 This struct defines the default types used for coordinates, intensity, charge, area, ...
		 
		 @ingroup Kernel
	*/
	typedef DoubleKernelTraits KernelTraits;
	
	/**
		 @brief Double traits class for the KERNEL classes.  See KernelTraits.
		 
		 This struct defines the types used for coordinates, intensity, charge, area, ...
		 This traits struct uses <code>double</code> for real numbers.
		 
		 @ingroup Kernel
	*/
	struct DoubleKernelTraits 
	{
	 public:
		
		/**\name Typedefs for primitive types used in the kernel classes
		 */
		//@{	
		/// Type for floating point numbers
		typedef DoubleReal RealType;
		/// Type for a coordinate
		typedef RealType CoordinateType;
		/// Type for an intensity, e.g. of a peak
		typedef RealType IntensityType;
		/// Type for a quality value
		typedef RealType QualityType;
		/// Type for a quality value (deprecated name, use QualityType instead)
		typedef QualityType RValueType;
		/// Type for a probability
		typedef RealType ProbabilityType;
		/// Type for priority
		typedef RealType PriorityType;
		/// Type for an area
		typedef RealType AreaType;
		/// Type for a full-width-at-half-max
		typedef RealType FullWidthHalfMaxType;
		/// Type for a width type
		typedef RealType WidthType;
		/// Type for signal to noise
		typedef RealType SignalToNoiseType;
		/// Type for a charge (signed)
		typedef SignedInt ChargeType;
		//@}

	 private:
		/// Constructor intentionally declared private -- instantiating this class makes no sense.
		DoubleKernelTraits();
		/// Constructor intentionally declared private -- instantiating this class makes no sense.
		DoubleKernelTraits(DoubleKernelTraits const &);
	};

	/**
		 @brief Float traits class for the KERNEL classes.  See KernelTraits.
		 
		 This struct defines the types used for coordinates, intensity, charge, area, ...
		 This traits struct uses <code>float</code> for real numbers.
		 
		 @ingroup Kernel
	*/
  struct FloatKernelTraits
  {
	 public:

		/**\name Typedefs for primitive types used in the kernel classes
		 */
		//@{
		/// Type for floating point numbers
		typedef Real RealType;
		/// Type for a coordinate
		typedef RealType CoordinateType;
		/// Type for an intensity, e.g. of a peak
		typedef RealType IntensityType;
		/// Type for a quality value
		typedef RealType QualityType;
		/// Type for a quality value (deprecated name, use QualityType instead)
		typedef QualityType RValueType;
		/// Type for a probability
		typedef RealType ProbabilityType;
		/// Type for priority
		typedef RealType PriorityType;
		/// Type for an area
		typedef RealType AreaType;	
		/// Type for a full-width-at-half-max
		typedef RealType FullWidthHalfMaxType;
		/// Type for a width type
		typedef RealType WidthType;
		/// Type for signal to noise
		typedef RealType SignalToNoiseType;
		/// Type for a charge (signed)
		typedef	SignedInt ChargeType;
		//@}    

	 private:
		/// Constructor intentionally declared private -- instantiating this class makes no sense.
		FloatKernelTraits();
		/// Constructor intentionally declared private -- instantiating this class makes no sense.
		FloatKernelTraits(FloatKernelTraits const &);
  };
  
} // namespace OpenMS

#endif // OPENMS_KERNEL_KERNELTRAITS_H
