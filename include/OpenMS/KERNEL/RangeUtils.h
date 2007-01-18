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

#ifndef OPENMS_KERNEL_RANGEUTILS_H
#define OPENMS_KERNEL_RANGEUTILS_H

#include <functional>
#include <algorithm>
#include <vector>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS 
{
	/**
		@defgroup Kernel Kernel
		
		@brief Kernel datastructures
	*/	
	
	/**
		@defgroup RangeUtils RangeUtils
		
		@brief Predicates for range operations
		
		@ingroup Kernel
		A group of predicates that can be used to perform range operations on MS data.
		They operate on classes that have the save interface as DSpectrum or DRawDataPoint respectively.
		<BR>
		<BR>
		The code for the removal of spectra in a certain retention time range from a vector of spectra might look like this:
		
		@code
		//data
		std::vector< DSpectrum<1> > spectra;
		
		//... spectra are added to the vector ...
		
		//range from 0.0 to 36.0 s
		RTRange< DSpectrum<1> > range(0.0, 36.0);

		//remove the range		
		spectra.erase(remove_if(spectra.begin(), spectra.end(), range), spectra.end());	
		@endcode
	
		The code for the removal of peaks in s certain intensity range from a spectrum might look like this:
		
		@code
		//data
		DSpectrum<1, DPeak<1> > spectrum;
		
		//... peaks are added to the spectrum ...
		
		//range from 0.0 to 5000.0 intensity
		IntensityRange< DPeak<1> >(0.0, 5000.0);
		
		//remove the range
		spectrum.erase(remove_if(spectrum.begin(), spectrum.end(), range), spectrum.end());
		@endcode
		
	*/
	
	/**
		@brief Predicate that determines if a spectrum lies inside/outside a specific retention time range
		
		SpectrumType must be a DSpectrum or have the same interface
		
		@ingroup RangeUtils
	*/
	template <class SpectrumType>
	class RTRange
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param min lower boundary
				@param max upper boundary
				@param reverse if @p reverse is true, operator() return true if the spectrum lies outside the range
			*/
			RTRange(double min, double max, bool reverse = false)
				: 
				min_(min),
				max_(max),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				double tmp = s.getRetentionTime();
				if (reverse_)
				{
					return ( min_ > tmp || max_ < tmp ); 
				}
				return ( min_ <= tmp && max_ >= tmp );
			}
		
		protected:
			double min_, max_;
			bool reverse_;
	};

	/**
		@brief Predicate that determines if a spectrum lies inside/outside a specific MS level set
		
		SpectrumType must be a DSpectrum or have the same interface
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class MSLevelRange
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param levels an array of MS levels
				@param reverse if @p reverse is true, operator() return true if the spectrum lies outside the set
			*/
			MSLevelRange(const std::vector<UnsignedInt>& levels, bool reverse = false)
				: 
				levels_(levels),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				UnsignedInt tmp = s.getMSLevel();
				if (reverse_)
				{
					return ( std::find(levels_.begin(), levels_.end(), tmp) == levels_.end() ); 
				}
				return ( std::find(levels_.begin(), levels_.end(), tmp) != levels_.end() );
			}
		
		protected:
			std::vector<UnsignedInt> levels_;
			bool reverse_;
	};

	/**
		@brief Predicate that determines if a spectrum has a certain scan mode
		
		SpectrumType must be a DSpectrum or have the same interface (SpectrumSettings)
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class ScanModePredicate
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param mode scan mode
				@param reverse if @p reverse is true, operator() return true if the spectrum has a different scan mode
			*/
			ScanModePredicate(SignedInt mode, bool reverse = false)
				: 
				mode_(mode),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				if (reverse_)
				{
					return s.getInstrumentSettings().getScanMode() != mode_ ; 
				}
				return s.getInstrumentSettings().getScanMode() == mode_ ; 
			}
		
		protected:
			SignedInt mode_;
			bool reverse_;
	};

	/**
		@brief Predicate that determines if a spectrum is empty.
		
		SpectrumType must have a size() member function
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class SpectrumEmptyPredicate
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param reverse if @p reverse is true, operator() return true if the spectrum is not empty
			*/
			SpectrumEmptyPredicate(bool reverse = false)
				: reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				if (reverse_)
				{
					return s.size() != 0 ; 
				}
				return s.size() == 0 ; 
			}
		
		protected:
			bool reverse_;
	};

	/**
		@brief Predicate that determines if a peak lies inside/outside a specific m/z range
		
		PeakType must be a DRawDataPoint or have the same interface.
		
		@note It is assumed that the m/z dimension is dimension 0!
		
		@ingroup RangeUtils
	*/		
	template <class PeakType>
	class MzRange
		: std::unary_function<PeakType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param min lower boundary
				@param max upper boundary
				@param reverse if @p reverse is true, operator() return true if the peak lies outside the set
			*/
			MzRange(double min, double max, bool reverse = false)
				: 
				min_(min),
				max_(max),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const PeakType& p) const
			{
				double tmp = p.getPosition()[0];
				if (reverse_)
				{
					return ( min_ > tmp || max_ < tmp ); 
				}
				return ( min_ <= tmp && max_ >= tmp );
			}
		
		protected:
			double min_, max_;
			bool reverse_;
	};
	
	/**
		@brief Predicate that determines if a peak lies inside/outside a specific intensity range
		
		PeakType must be a DRawDataPoint or have the same interface 
		
		@ingroup RangeUtils
	*/	
	template <class PeakType>
	class IntensityRange
		: std::unary_function<PeakType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param min lower boundary
				@param max upper boundary
				@param reverse if @p reverse is true, operator() return true if the peak lies outside the set
			*/
			IntensityRange(double min, double max, bool reverse = false)
				: 
				min_(min),
				max_(max),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const PeakType& p) const
			{
				double tmp = p.getIntensity();
				if (reverse_)
				{
					return ( min_ > tmp || max_ < tmp ); 
				}
				return ( min_ <= tmp && max_ >= tmp );
			}
		
		protected:
			double min_, max_;
			bool reverse_;
	};

} // namespace OpenMS

#endif // OPENMS_KERNEL_RANGEUTILS_H



