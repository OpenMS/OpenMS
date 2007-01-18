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

#ifndef OPENMS_KERNEL_STANDARDTYPES_H
#define OPENMS_KERNEL_STANDARDTYPES_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DFeatureMap.h>

namespace OpenMS
{
	/**
		@defgroup default_types Default MS data types
		
		@brief OpenMS default data types
		
    The following classes and type definitions provide some
		default data types to handle one- and two-dimensional
		raw data, peak (stick) data, and feature data (maps).
		Depending on your application, the use of other internal
		data structures might be preferable, although we recommend
		to use these data types if you have no urgent need to use					
		something else.
		
		@ingroup Kernel
	*/

	//@{	
	/**
		 @brief One-dimensional raw data point
	 */
	typedef DRawDataPoint<1,KernelTraits> RawDataPoint;

	/**
		 @brief Two-dimensional raw data point
	*/
	typedef DRawDataPoint<2,KernelTraits> RawDataPoint2D;

	/**
		 @brief Spectrum consisting of raw data points, with meta information.

		 Meta information includes retention time and MS-level.
	*/
	typedef MSSpectrum<RawDataPoint> RawSpectrum;
	/**
		 @brief Two-dimensional map of raw data points, with meta information about experimental settings.
	*/
	typedef MSExperiment<RawDataPoint> RawMap;

	/**
		 @brief One-dimensional peak
	*/
	typedef DPeak<1, KernelTraits> Peak;
	/**
		 @brief Two-dimensional peak
	*/
	typedef DPeak<2, KernelTraits> Peak2D;
	/**
		 @brief Spectrum consisting of peaks with meta information.
	*/
	typedef MSSpectrum<Peak> PeakSpectrum;
	/**
		 @brief  Two-dimensional map of peaks, with meta information about experimental settings.
	*/
	typedef MSExperiment<Peak> PeakMap;

	/**
		 @brief Two-dimensional feature.
	*/
	typedef DFeature<2, KernelTraits> Feature;
	/**
		 @brief  Two-dimensional map of features, with meta information about experimental settings.
	*/
	typedef DFeatureMap<2, Feature> FeatureMap;
	//@}

}

#endif // OPENMS_KERNEL_STANDARDTYPES_H
