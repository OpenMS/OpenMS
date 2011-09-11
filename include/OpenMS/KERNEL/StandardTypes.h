// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_KERNEL_STANDARDTYPES_H
#define OPENMS_KERNEL_STANDARDTYPES_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/KERNEL/MSExperiment.h>

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


	/**
		@brief Chromatogram consisting of raw data points or peaks

		@ingroup Kernel
	*/
	typedef MSChromatogram<ChromatogramPeak> Chromatogram;
	//@}

}

#endif // OPENMS_KERNEL_STANDARDTYPES_H
