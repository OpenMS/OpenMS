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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_RANGEUTILS_H
#define OPENMS_KERNEL_RANGEUTILS_H

#include <functional>
#include <algorithm>
#include <vector>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/IntList.h>
#include <OpenMS/METADATA/Precursor.h>

namespace OpenMS 
{	
	/**
		@defgroup RangeUtils RangeUtils
		
		@brief Predicates for range operations
		
		@ingroup Kernel
		A group of predicates that can be used to perform range operations on MS data.
		They operate on classes that have the save interface as Spectrum or Peak1D or Peak2D, respectively.
		<BR>
		<BR>
		The code for the removal of spectra in a certain retention time range from a vector of spectra might look like this:
		
		@code
		//data
		std::vector< MSSpectrum<> > spectra;
		
		//... spectra are added to the vector ...
		
		//range from 0.0 to 36.0 s
		InRTRange< MSSpectrum<> > range(0.0, 36.0);

		//remove the range		
		spectra.erase(remove_if(spectra.begin(), spectra.end(), range), spectra.end());	
		@endcode
	
		The code for the removal of peaks within certain intensity range from a spectrum might look like this:
		
		@code
		//data
		MSSpectrum<> spectrum;
		
		//... peaks are added to the spectrum ...
		
		//range from 0.0 to 5000.0 intensity
		InIntensityRange range< Peak1D >(0.0, 5000.0);
		
		//remove the range
		spectrum.erase(remove_if(spectrum.begin(), spectrum.end(), range), spectrum.end());
		@endcode
	*/


	/**
		@brief Predicate that determines if a class has a certain metavalue
		
		MetaContainer must be a MetaInfoInterface or have the same interface
		
		@ingroup RangeUtils
	*/
	template <class MetaContainer>
	class HasMetaValue
		: std::unary_function<MetaContainer, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param metavalue MetaValue that needs to be present
				@param reverse if @p reverse is true, operator() returns true if the metavalue does not exist
			*/
			HasMetaValue(String metavalue, bool reverse = false)
				: 
				metavalue_key_(metavalue),
				reverse_(reverse)
			{
			}
		
			inline bool operator()(const MetaContainer& s) const
			{
				bool has_meta_value = s.metaValueExists(metavalue_key_);
				if (reverse_)
				{
					return !has_meta_value; 
				}
				return has_meta_value;
			}
		
		protected:
			String metavalue_key_;
			bool reverse_;
	};


	/**
		@brief Predicate that determines if a spectrum lies inside/outside a specific retention time range
		
		SpectrumType must be a Spectrum or have the same interface
		
		@ingroup RangeUtils
	*/
	template <class SpectrumType>
	class InRTRange
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param min lower boundary
				@param max upper boundary
				@param reverse if @p reverse is true, operator() returns true if the spectrum lies outside the range
			*/
			InRTRange(DoubleReal min, DoubleReal max, bool reverse = false)
				: 
				min_(min),
				max_(max),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				DoubleReal tmp = s.getRT();
				if (reverse_)
				{
					return ( min_ > tmp || max_ < tmp ); 
				}
				return ( min_ <= tmp && max_ >= tmp );
			}
		
		protected:
			DoubleReal min_, max_;
			bool reverse_;
	};

	/**
		@brief Predicate that determines if a spectrum lies inside/outside a specific MS level set
		
		SpectrumType must be a Spectrum or have the same interface
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class InMSLevelRange
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param levels an array of MS levels
				@param reverse if @p reverse is true, operator() returns true if the spectrum lies outside the set
			*/
			InMSLevelRange(const IntList& levels, bool reverse = false)
				: 
				levels_(levels),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				Int tmp = s.getMSLevel();
				if (reverse_)
				{
					return ( std::find(levels_.begin(), levels_.end(), tmp) == levels_.end() ); 
				}
				return ( std::find(levels_.begin(), levels_.end(), tmp) != levels_.end() );
			}
		
		protected:
			IntList levels_;
			bool reverse_;
	};

	/**
		@brief Predicate that determines if a spectrum has a certain scan mode
		
		SpectrumType must be a Spectrum or have the same interface (SpectrumSettings)
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class HasScanMode
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param mode scan mode
				@param reverse if @p reverse is true, operator() returns true if the spectrum has a different scan mode
			*/
			HasScanMode(Int mode, bool reverse = false)
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
			Int mode_;
			bool reverse_;
	};

	/**
		@brief Predicate that determines if a spectrum is empty.
		
		SpectrumType must have a empty() method
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class IsEmptySpectrum
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param reverse if @p reverse is true, operator() returns true if the spectrum is not empty
			*/
			IsEmptySpectrum(bool reverse = false)
				: reverse_(reverse)
			{
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				if (reverse_)
				{
					return !s.empty(); 
				}
				return s.empty(); 
			}
		
		protected:
			bool reverse_;
	};

	/**
		@brief Predicate that determines if a spectrum is a zoom (enhanced resolution) spectrum.
		
		SpectrumType must have a getInstrumentSettings() method
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class IsZoomSpectrum
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param reverse if @p reverse is true, operator() returns true if the spectrum is not a zoom spectrum
			*/
			IsZoomSpectrum(bool reverse = false)
				: reverse_(reverse)
			{
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				if (reverse_)
				{
					return !s.getInstrumentSettings().getZoomScan(); 
				}
				return s.getInstrumentSettings().getZoomScan(); 
			}
		
		protected:
			bool reverse_;
	};


	/**
		@brief Predicate that determines if a spectrum was generated using any activation method given in the C'tor list
		
		SpectrumType must have a getPrecursors() method
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class HasActivationMethod
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param methods List of methods that is compared against precursor activation methods
				@param reverse if @p reverse is true, operator() returns true if the spectrum is not using one of the specified activation methods
			*/
			HasActivationMethod(const StringList& methods, bool reverse = false)
				: methods_(methods),
					reverse_(reverse)
			{
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
				for (std::vector< Precursor >::const_iterator it = s.getPrecursors().begin(); it!=s.getPrecursors().end(); ++it)
				{
					for (std::set< Precursor::ActivationMethod >::const_iterator it_a = it->getActivationMethods().begin();
							 it_a != it->getActivationMethods().end();
							 ++it_a)
					{
						if (methods_.contains(Precursor::NamesOfActivationMethod[*it_a]))
						{
							// found matching activation method
							if (reverse_)	return false;
							else return true; 
						}
					}
				}
				
				if (reverse_)	return true;
				else return false; 
			}
		
		protected:
			StringList methods_;
			bool reverse_;
	};


	/**
		@brief Predicate that determines if a spectrum has a certain precursor charge as given in the C'tor list
		
		SpectrumType must have a getPrecursors() method
		
		@ingroup RangeUtils
	*/	
	template <class SpectrumType>
	class HasPrecursorCharge
		: std::unary_function<SpectrumType, bool>
	{
		public:
			/**
				@brief Constructor
				
        @param charges List of charges that is compared against precursor charge
				@param reverse if @p reverse is true, operator() returns true if the spectrum has not one of the specified precursor charges
			*/
			HasPrecursorCharge(const IntList& charges, bool reverse = false)
				: charges_(charges),
					reverse_(reverse)
			{
			}
		
			inline bool operator()(const SpectrumType& s) const
			{
        bool match = false;
				for (std::vector< Precursor >::const_iterator it = s.getPrecursors().begin(); it!=s.getPrecursors().end(); ++it)
				{

          Int tmp = it->getCharge();
          match = match || (std::find(charges_.begin(), charges_.end(), tmp) != charges_.end());
				}

				if (reverse_)	return (!match);
				else return match; 
			}
		
		protected:
			IntList charges_;
			bool reverse_;
	};


	/**
		@brief Predicate that determines if a peak lies inside/outside a specific m/z range
		
		PeakType must have a getPosition() method.
		
		@note It is assumed that the m/z dimension is dimension 0!
		
		@ingroup RangeUtils
	*/		
	template <class PeakType>
	class InMzRange
		: std::unary_function<PeakType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param min lower boundary
				@param max upper boundary
				@param reverse if @p reverse is true, operator() returns true if the peak lies outside the range
			*/
			InMzRange(DoubleReal min, DoubleReal max, bool reverse = false)
				: 
				min_(min),
				max_(max),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const PeakType& p) const
			{
				DoubleReal tmp = p.getPosition()[0];
				if (reverse_)
				{
					return ( min_ > tmp || max_ < tmp ); 
				}
				return ( min_ <= tmp && max_ >= tmp );
			}
		
		protected:
			DoubleReal min_, max_;
			bool reverse_;
	};
	
	/**
		@brief Predicate that determines if a peak lies inside/outside a specific intensity range
		
		PeakType must have a getIntensity() method.
		
		@ingroup RangeUtils
	*/	
	template <class PeakType>
	class InIntensityRange
		: std::unary_function<PeakType, bool>
	{
		public:
			/**
				@brief Constructor
				
				@param min lower boundary
				@param max upper boundary
				@param reverse if @p reverse is true, operator() returns true if the peak lies outside the set
			*/
			InIntensityRange(DoubleReal min, DoubleReal max, bool reverse = false)
				: 
				min_(min),
				max_(max),
				reverse_(reverse)
			{
				
			}
		
			inline bool operator()(const PeakType& p) const
			{
				DoubleReal tmp = p.getIntensity();
				if (reverse_)
				{
					return ( min_ > tmp || max_ < tmp ); 
				}
				return ( min_ <= tmp && max_ >= tmp );
			}
		
		protected:
			DoubleReal min_, max_;
			bool reverse_;
	};

} // namespace OpenMS

#endif // OPENMS_KERNEL_RANGEUTILS_H



