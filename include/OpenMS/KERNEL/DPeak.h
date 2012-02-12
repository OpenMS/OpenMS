// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_DPEAK_H
#define OPENMS_KERNEL_DPEAK_H

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>

namespace OpenMS
{
	
	/**
		@brief Metafunction to choose among Peak1D respectively Peak2D through a template argument.
		
		The result is accessible via typedef Type:
		- @c DPeak<1>::Type is @c Peak1D
		- @c DPeak<2>::Type is @c Peak2D
	
		Example:
	  @code
		template class BaseModel<UInt D>
		{
			// BaseModel<D>::PeakType is either Peak1D or Peak2D, depending on D
			typedef typename DPeak<D>::Type PeakType;
		};
		@endcode
	*/
	template <UInt dimensions>
	struct DPeak
	{
	};

	// We do not want these classes to show up in the docu
	/// @cond HIDDENSTUFF
	
	template <>
	struct DPeak <1>
	{
		typedef Peak1D Type;
	};

	template <>
	struct DPeak <2>
	{
		typedef Peak2D Type;
	};

	/// @endcond

} // namespace OpenMS

#endif // OPENMS_KERNEL_DPEAK_H
