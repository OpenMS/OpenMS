// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FILTERING_BASELINE_MORPHOLOGICALFILTER_H
#define OPENMS_FILTERING_BASELINE_MORPHOLOGICALFILTER_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <algorithm>
#include <iterator>

namespace OpenMS
{

	namespace Internal
	{
		/**
			@brief An iterator wrapper to access peak intensities instead of the peak itself.
			
			It is using unary operator *, and the like.  This is not a full implementation of the
			iterator concept, it can only do what is needed for MorphologicalFilter.
		*/
		template <typename IteratorT>
		class /* OPENMS_DLLAPI */ IntensityIteratorWrapper
			: public std::iterator<std::forward_iterator_tag, typename IteratorT::value_type::IntensityType>
		{
			public:
				typedef typename IteratorT::value_type::IntensityType value_type;
				typedef typename IteratorT::value_type::IntensityType& reference;
				typedef typename IteratorT::value_type::IntensityType* pointer;
				typedef typename IteratorT::difference_type difference_type;

				IntensityIteratorWrapper( const IteratorT& rhs )
					: base(rhs)
				{
				}

				value_type operator*()
				{
					return base->getIntensity();
				}
				
				template <typename IndexT> value_type operator[](const IndexT& index)
				{
					return base[index].getIntensity();
				}

				difference_type operator-(IntensityIteratorWrapper& rhs) const
				{
					return base - rhs.base;
				}

				IntensityIteratorWrapper& operator++()
				{
					++base;
					return *this;
				}

				IntensityIteratorWrapper operator++(int)
				{
					IteratorT tmp = *this;
					++(*this);
					return tmp;
				}

				bool operator!=(const IntensityIteratorWrapper&	rhs) const
				{
					return base!=rhs.base;
				}

			protected:
				IteratorT base;
		};

		/// make-function so that we need no write out all those type names to get the wrapped iterator.
		template <typename IteratorT>
		IntensityIteratorWrapper<IteratorT> intensityIteratorWrapper(const IteratorT& rhs)
		{
			return IntensityIteratorWrapper<IteratorT>(rhs);
		}

	}

	/**
		@brief This class implements baseline filtering operations using methods
		from mathematical morphology.
	
		The fundamental operations are erosion and dilation.  These are defined with
		respect to a structuring element.  In our case, this is just a straight line
		and the definitions can be given as follows:
	
		Assume that the input is \f$x_0, x_1, x_2, ...\f$.  Then the <i>erosion</i>
		of \f$x\f$ contains the minima of a sliding window of size struc_size around
		\f$ i \f$, i.e. \f[ \mathrm{erosion}_i = \min\{x_{i-\mathrm{struc\_size}/2},
		\ldots, x_{i+\mathrm{struc\_size}/2}\} \f].  The <i>dilation</i> of \f$x\f$
		contains the maxima of a sliding window of size struc_size around \f$ i \f$,
		i.e. \f[ \mathrm{dilation}_i = \max\{x_{i-\mathrm{struc\_size}/2}, \ldots,
		x_{i+\mathrm{struc\_size}/2}\} \f].
	
		For morphological baseline filtering the <i>tophat</i> method is used.  The
		tophat transform is defined as signal minus opening, where the opening is
		the dilation of the erosion of the signal.
	
		@image html MorphologicalFilter_tophat.png
	
		Several other morphological operations are implemented as well.  See the
		image below and the documentation for further explanation.
	
		@image html MorphologicalFilter_all.png
	
		@note The class #MorphologicalFilter is designed for uniformly spaced profile data.
		
		@note The data must be sorted according to ascending m/z!
		
		@htmlinclude OpenMS_MorphologicalFilter.parameters
	
		@ingroup SignalProcessing
	*/
	class OPENMS_DLLAPI MorphologicalFilter
		:	public ProgressLogger,
			public DefaultParamHandler
	{
		public:
	
			/// Constructor
			MorphologicalFilter()
				: ProgressLogger(),
					DefaultParamHandler("MorphologicalFilter"),
					struct_size_in_datapoints_(0)
			{
				//structuring element
				defaults_.setValue("struc_elem_length",3.0,"Length of the structuring element. This should be wider than the expected peak width.");
				defaults_.setValue("struc_elem_unit","Thomson","The unit of the 'struct_elem_length'.");
				defaults_.setValidStrings("struc_elem_unit",StringList::create("Thomson,DataPoints"));
				//methods
				defaults_.setValue("method","tophat","Method to use, the default is 'tophat'.  Do not change this unless you know what you are doing.  The other methods may be useful for tuning the parameters, see the class documentation of MorpthologicalFilter.");
				defaults_.setValidStrings("method", StringList::create("identity,erosion,dilation,opening,closing,gradient,tophat,bothat,erosion_simple,dilation_simple"));
	
				defaultsToParam_();
			}
	
			/// Destructor
			virtual ~MorphologicalFilter()
			{
			}
	
			/** @brief Applies the morphological filtering operation to an iterator range.
			 
			Input and output range must be valid, i.e. allocated before.
			InputIterator must be a random access iterator type.
	
			@param input_begin the begin of the input range
			@param input_end  the end of the input range
			@param output_begin the begin of the output range
	
			@exception Exception::IllegalArgument The given method is not one of the values defined in the @em method paramter.
			*/
			template < typename InputIterator, typename OutputIterator >
			void filterRange(InputIterator input_begin, InputIterator input_end, OutputIterator output_begin)
			{
				// the buffer is static only to avoid reallocation
				static std::vector< typename InputIterator::value_type > buffer;
				const UInt size = input_end - input_begin;
				
				//determine the struct size in data points if not already set
				if (struct_size_in_datapoints_==0)
				{
					struct_size_in_datapoints_ = (UInt)(DoubleReal)param_.getValue("struc_elem_length");
				}
				
				//apply the filtering
				String method = param_.getValue("method");			
				if (method=="identity")
				{
						std::copy(input_begin,input_end,output_begin);
				}
				else if (method=="erosion")
				{
					applyErosion_(struct_size_in_datapoints_,input_begin,input_end,output_begin);
				}
				else if (method=="dilation")
				{
					applyDilation_(struct_size_in_datapoints_,input_begin,input_end,output_begin);
				}
				else if (method=="opening")
				{
					if ( buffer.size() < size ) buffer.resize(size);
					applyErosion_(struct_size_in_datapoints_,input_begin,input_end,buffer.begin());
					applyDilation_(struct_size_in_datapoints_,buffer.begin(),buffer.begin()+size,output_begin);
				}
				else if (method=="closing")
				{
					if ( buffer.size() < size ) buffer.resize(size);
					applyDilation_(struct_size_in_datapoints_,input_begin,input_end,buffer.begin());
					applyErosion_(struct_size_in_datapoints_,buffer.begin(),buffer.begin()+size,output_begin);
				}
				else if (method=="gradient")
				{
					if ( buffer.size() < size ) buffer.resize(size);
					applyErosion_(struct_size_in_datapoints_,input_begin,input_end,buffer.begin());
					applyDilation_(struct_size_in_datapoints_,input_begin,input_end,output_begin);
					for ( UInt i = 0; i < size; ++i ) output_begin[i] -= buffer[i];
				}
				else if (method=="tophat")
				{
					if ( buffer.size() < size ) buffer.resize(size);
					applyErosion_(struct_size_in_datapoints_,input_begin,input_end,buffer.begin());
					applyDilation_(struct_size_in_datapoints_,buffer.begin(),buffer.begin()+size,output_begin);
					for ( UInt i = 0; i < size; ++i ) output_begin[i] = input_begin[i] - output_begin[i];
				}
				else if (method=="bothat")
				{
					if ( buffer.size() < size ) buffer.resize(size);
					applyDilation_(struct_size_in_datapoints_,input_begin,input_end,buffer.begin());
					applyErosion_(struct_size_in_datapoints_,buffer.begin(),buffer.begin()+size,output_begin);
					for ( UInt i = 0; i < size; ++i ) output_begin[i] = input_begin[i] - output_begin[i];
				}
				else if (method=="erosion_simple")
				{
					applyErosionSimple_(struct_size_in_datapoints_,input_begin,input_end,output_begin);
				}
				else if (method=="dilation_simple")
				{
					applyDilationSimple_(struct_size_in_datapoints_,input_begin,input_end,output_begin);
				}
				
				struct_size_in_datapoints_ = 0;
			}
	
	
			/**
				@brief Applies the morphological filtering operation to an MSSpectrum.
	
				If the size of the structuring element is given in 'Thomson', the number of data points for
				the structuring element is computed as follows:
				<ul>
					<li>The data points are assumed to be uniformly spaced.  We compute the
						average spacing from the position of the first and the last peak and the
						total number of peaks in the input range.
					<li>The number of data points in the structuring element is computed
						from struc_size and the average spacing, and rounded up to an odd
						number.
				</ul>
			*/
			template <typename PeakType>
			void filter(MSSpectrum<PeakType>& spectrum)
			{
				//make sure the right peak type is set
				spectrum.setType(SpectrumSettings::RAWDATA);
				
				//Abort if there is nothing to do
				if ( spectrum.size() <= 1 ) return;
	
				//Determine structuring element size in datapoints (depending on the unit)
				if ( (String)(param_.getValue("struc_elem_unit"))=="Thomson" )
				{
					struct_size_in_datapoints_ =
						UInt (
									ceil (
												(DoubleReal)(param_.getValue("struc_elem_length"))
												*
												DoubleReal( spectrum.size() - 1 )
												/
												( spectrum.back().getMZ() - spectrum.begin()->getMZ() )
											 )
								 );
				}
				else
				{
					struct_size_in_datapoints_ = (UInt)(DoubleReal)param_.getValue("struc_elem_length");
				}
				//make it odd (needed for the algorithm)
				if ( ! Math::isOdd(struct_size_in_datapoints_) ) ++struct_size_in_datapoints_;
				
				//apply the filtering and overwrite the input data
				std::vector<typename PeakType::IntensityType> output(spectrum.size());
				filterRange( Internal::intensityIteratorWrapper(spectrum.begin()),
										 Internal::intensityIteratorWrapper(spectrum.end()),
										 output.begin()
									 );
				
				//overwrite output with data
				for (Size i=0; i<spectrum.size(); ++i)
				{
					spectrum[i].setIntensity(output[i]);
				}
			}
	
	
			/**
				@brief Applies the morphological filtering operation to an MSExperiment.
				
				The size of the structuring element is computed for each spectrum individually, if it is given in 'Thomson'.
				See the filtering method for MSSpectrum for details.
			*/
			template <typename PeakType>
			void filterExperiment(MSExperiment<PeakType>& exp)
			{
				startProgress(0,exp.size(),"filtering baseline");
				for ( UInt i = 0; i < exp.size(); ++i )
				{
					filter(exp[i]);
					setProgress(i);
				}
				endProgress();
			}
		
		protected:
		
			///Member for struct size in data points
			UInt struct_size_in_datapoints_;
	
			/** @brief Applies erosion.  This implementation uses van Herk's method.
			Only 3 min/max comparisons are required per data point, independent of
			struc_size.
			*/
			template < typename InputIterator, typename OutputIterator >
				void applyErosion_( Int struc_size, InputIterator input, InputIterator input_end, OutputIterator output )
			{
				typedef typename InputIterator::value_type ValueType;
				const Int size = input_end - input;
				const Int struc_size_half = struc_size / 2; // yes, integer division
	
				static std::vector<ValueType> buffer;
				if ( Int(buffer.size()) < struc_size) buffer.resize(struc_size);
	
				Int anchor; // anchoring position of the current block
				Int i;      // index relative to anchor, used for 'for' loops
				Int ii = 0; // input index
				Int oi = 0; // output index
				ValueType current; // current value
	
				// we just can't get the case distinctions right in these cases, resorting to simple method.
				if ( size <= struc_size || size <= 5 )
				{
					applyErosionSimple_(struc_size,input,input_end,output);
					return;
				}
				{
					// lower margin area
					current = input[0];
					for ( ++ii; ii < struc_size_half; ++ii ) if ( current > input[ii] ) current = input[ii];
					for ( ; ii < std::min(Int(struc_size),size); ++ii, ++oi )
					{
						if ( current > input[ii] ) current = input[ii];
						output[oi] = current;
					}
				}
				{
					// middle (main) area
					for ( anchor = struc_size;
								anchor <= size - struc_size;
								anchor += struc_size
							)
					{
						ii = anchor;
						current = input[ii];
						buffer[0] = current;
						for ( i = 1; i < struc_size; ++i, ++ii )
						{
							if ( current > input[ii] ) current = input[ii];
							buffer[i] = current;
						}
						ii = anchor - 1;
						oi = ii + struc_size_half;
						current = input[ii];
						for ( i = 1; i < struc_size; ++i, --ii, --oi )
						{
							if ( current > input[ii] ) current = input[ii];
							output[oi] = std::min( buffer[struc_size-i], current );
						}
						if ( current > input[ii] ) current = input[ii];
						output[oi] = current;
					}
				}
				{
					// higher margin area
					ii = size - 1;
					oi = ii;
					current = input[ii];
					for ( --ii; ii >= size - struc_size_half; --ii ) if ( current > input[ii] ) current = input[ii];
					for ( ; ii >= std::max(size - Int(struc_size),0); --ii, --oi )
					{
						if ( current > input[ii] ) current = input[ii];
						output[oi] = current;
					}
					anchor = size - struc_size;
					ii = anchor;
					current = input[ii];
					buffer[0] = current;
					for ( i = 1; i < struc_size; ++i, ++ii )
					{
						if ( current > input[ii] ) current = input[ii];
						buffer[i] = current;
					}
					ii = anchor - 1;
					oi = ii + struc_size_half;
				  current = input[ii];
					for ( i = 1; (ii >= 0) && (i < struc_size) ; ++i, --ii, --oi )
					{
						if ( current > input[ii] ) current = input[ii];
						output[oi] = std::min( buffer[struc_size-i], current );
					}
					if ( ii >= 0 )
					{
						if ( current > input[ii] ) current = input[ii];
						output[oi] = current;
					}
				}
				return;
			}
	
			/** @brief Applies dilation.  This implementation uses van Herk's method.
			Only 3 min/max comparisons are required per data point, independent of
			struc_size.
			*/
			template < typename InputIterator, typename OutputIterator >
			void applyDilation_( Int struc_size, InputIterator input, InputIterator input_end, OutputIterator output )
			{
				typedef typename InputIterator::value_type ValueType;
				const Int size = input_end - input;
				const Int struc_size_half = struc_size / 2; // yes, integer division
	
				static std::vector<ValueType> buffer;
				if ( Int(buffer.size()) < struc_size) buffer.resize(struc_size);
	
				Int anchor; // anchoring position of the current block
				Int i;      // index relative to anchor, used for 'for' loops
				Int ii = 0; // input index
				Int oi = 0; // output index
				ValueType current; // current value
	
				// we just can't get the case distinctions right in these cases, resorting to simple method.
				if ( size <= struc_size || size <= 5 )
				{
					applyDilationSimple_(struc_size,input,input_end,output);
					return;
				}
				{
					// lower margin area
					current = input[0];
					for ( ++ii; ii < struc_size_half; ++ii ) if ( current < input[ii] ) current = input[ii];
					for ( ; ii < std::min(Int(struc_size),size); ++ii, ++oi )
					{
						if ( current < input[ii] ) current = input[ii];
						output[oi] = current;
					}
				}
				{
					// middle (main) area
					for ( anchor = struc_size;
								anchor <= size - struc_size;
								anchor += struc_size
							)
					{
						ii = anchor;
						current = input[ii];
						buffer[0] = current;
						for ( i = 1; i < struc_size; ++i, ++ii )
						{
							if ( current < input[ii] ) current = input[ii];
							buffer[i] = current;
						}
						ii = anchor - 1;
						oi = ii + struc_size_half;
						current = input[ii];
						for ( i = 1; i < struc_size; ++i, --ii, --oi )
						{
							if ( current < input[ii] ) current = input[ii];
							output[oi] = std::max( buffer[struc_size-i], current );
						}
						if ( current < input[ii] ) current = input[ii];
						output[oi] = current;
					}
				}
				{
					// higher margin area
					ii = size - 1;
					oi = ii;
					current = input[ii];
					for ( --ii; ii >= size - struc_size_half; --ii ) if ( current < input[ii] ) current = input[ii];
					for ( ; ii >= std::max(size - Int(struc_size),0); --ii, --oi )
					{
						if ( current < input[ii] ) current = input[ii];
						output[oi] = current;
					}
					anchor = size - struc_size;
					ii = anchor;
					current = input[ii];
					buffer[0] = current;
					for ( i = 1; i < struc_size; ++i, ++ii )
					{
						if ( current < input[ii] ) current = input[ii];
						buffer[i] = current;
					}
					ii = anchor - 1;
					oi = ii + struc_size_half;
				  current = input[ii];
					for ( i = 1; (ii >= 0) && (i < struc_size) ; ++i, --ii, --oi )
					{
						if ( current < input[ii] ) current = input[ii];
						output[oi] = std::max( buffer[struc_size-i], current );
					}
					if ( ii >= 0 )
					{
						if ( current < input[ii] ) current = input[ii];
						output[oi] = current;
					}
				}
				return;
			}
	
			/// Applies erosion.  Simple implementation, possibly faster if struc_size is very small, and used in some special cases.
			template < typename InputIterator, typename OutputIterator >
			void applyErosionSimple_( Int struc_size, InputIterator input_begin, InputIterator input_end, OutputIterator output_begin )
			{
				typedef typename InputIterator::value_type ValueType;
				const int size = input_end - input_begin;
				const Int struc_size_half = struc_size / 2; // yes integer division
				for ( Int index = 0; index < size; ++ index )
				{
					Int start = std::max( 0,    index - struc_size_half );
					Int stop  = std::min( size - 1, index + struc_size_half );
					ValueType value = input_begin[start];
					for ( Int i = start + 1; i <= stop; ++i ) if ( value > input_begin[i] ) value = input_begin[i];
					output_begin[index] = value;
				}
				return;
			}
	
			/// Applies dilation.  Simple implementation, possibly faster if struc_size is very small, and used in some special cases.
			template < typename InputIterator, typename OutputIterator >
			void applyDilationSimple_( Int struc_size, InputIterator input_begin, InputIterator input_end, OutputIterator output_begin )
			{
				typedef typename InputIterator::value_type ValueType;
				const int size = input_end - input_begin;
				const Int struc_size_half = struc_size / 2; // yes integer division
				for ( Int index = 0; index < size; ++ index )
				{
					Int start = std::max( 0,    index - struc_size_half );
					Int stop   = std::min( size - 1, index + struc_size_half );
					ValueType value = input_begin[start];
					for ( Int i = start + 1; i <= stop; ++i ) if ( value < input_begin[i] ) value = input_begin[i];
					output_begin[index] = value;
				}
				return;
			}
	
		private:
			
			/// copy constructor not implemented
			MorphologicalFilter(const MorphologicalFilter& source);

	};

}// namespace OpenMS

#endif
