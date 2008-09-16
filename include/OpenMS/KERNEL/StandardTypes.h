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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_STANDARDTYPES_H
#define OPENMS_KERNEL_STANDARDTYPES_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/KERNEL/RichPeak2D.h>

namespace OpenMS
{
	//@{	
	/**
		@brief Spectrum consisting of raw data points or peaks.
		
		Meta information includes retention time and MS-level.
		
		@ingroup Kernel
	*/
	typedef MSSpectrum<Peak1D> PeakSpectrum;
	/**
		@brief Two-dimensional map of raw data points or peaks.
		
		@ingroup Kernel
	*/
	typedef MSExperiment<Peak1D> PeakMap;

	/**
		@brief Spectrum consisting of raw data points or peaks with meta information.
		
		@ingroup Kernel
	*/
	typedef MSSpectrum<RichPeak1D> RichPeakSpectrum;
	
	/**
		@brief  Two-dimensional map of raw data points or peaks with meta information.
		
		@ingroup Kernel
	*/
	typedef MSExperiment<RichPeak1D> RichPeakMap;

	//@}

	/**
	@brief Metafunction to choose among Peak1D respectively Peak2D through a
	template argument.  The result is accessible via typedef Type.

	- @c PeakXD<1>::Type is @c Peak1D
	- @c PeakXD<2>::Type is @c Peak2D
	.

	Example:
  @code
	template class BaseModel<UInt D>
	{
	;;;
	// BaseModel<D>::PeakType is either Peak1D or Peak2D, depending on D
	typedef typename PeakXD<D>::Type PeakType;
	;;;
	};
	@endcode

	@internal @note PeakXD is default-defined as an empty class, rather than declared as a template only, since doxygen won't document this otherwise.  :-(

	*/
	template <UInt Dimensions>
	struct PeakXD
	{};
	
	// we do not want the template specializations to show up in normal docu:
	/// @cond INTERNAL_INFO

	/// @c PeakXD<1>::Type is @c Peak1D
	template <>
	struct PeakXD <1>
	{
		typedef Peak1D Type;
	};

	/// @c PeakXD<2>::Type is @c Peak2D
	template <>
	struct PeakXD <2>
	{
		typedef Peak2D Type;
	};
	
	/// @endcond

	/**
	@brief Metafunction to choose among RichPeak1D respectively RichPeak2D through a
	template argument.  The result is accessible via typedef Type.

	- @c RichPeakXD<1>::Type is @c RichPeak1D
	- @c RichPeakXD<2>::Type is @c RichPeak2D
	.

	Example: see @c PeakXD

	@internal @note RichPeakXD is default-defined as an empty class, rather than declared as a template only, since doxygen won't document this otherwise.  :-(

	*/
	template <UInt Dimensions>
	struct RichPeakXD
	{};
	
	// we do not want the template specializations to show up in normal docu:
	/// @cond INTERNAL_INFO

	/// @c RichPeakXD<1>::Type is @c RichPeak1D
	template <>
	struct RichPeakXD <1>
	{
		typedef RichPeak1D Type;
	};

	/// @c RichPeakXD<2>::Type is @c RichPeak2D
	template <>
	struct RichPeakXD <2>
	{
		typedef RichPeak2D Type;
	};
	
	/// @endcond

}

#endif // OPENMS_KERNEL_STANDARDTYPES_H
