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
//

#ifndef OPENMS_FILTERING_BASELINE_MORPHOLOGICALFILTER_H
#define OPENMS_FILTERING_BASELINE_MORPHOLOGICALFILTER_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <algorithm>

namespace OpenMS
{

	namespace Internal
	{
		/**@brief An iterator wrapper to access peak intensities using unary
		operator *, and the like.  This is not really an implementation of the
		iterator concept, it can only do what is needed for MorphologicalFilter.
		*/
		template <typename IteratorT>
		struct /* OPENMS_DLLAPI */ IntensityIteratorWrapper : IteratorT
		{
			typedef typename IteratorT::value_type::IntensityType value_type;
		
			/// "Why this?" - We need getIntensity return a reference, not a value.
			struct Peak1D_ : Peak1D
			{
				Peak1D::IntensityType& getIntensity() { return intensity_; }
			};

			IntensityIteratorWrapper( const IteratorT& rhs ) : IteratorT(rhs) {}

			// To avoid further complication this provides mutable access even if IteratorT is const.
			value_type & operator*()
			{
				return static_cast<Peak1D_&>(const_cast<Peak1D&>(IteratorT::operator*())).getIntensity();
			}

			// To avoid further complication this provides mutable access even if IteratorT is const.
			template <typename IndexT> value_type & operator[](const IndexT& rhs)
			{
				return static_cast<Peak1D_&>(const_cast<Peak1D&>(IteratorT::operator[](rhs))).getIntensity();
			}
			
			template <typename IndexT> IntensityIteratorWrapper operator+(const IndexT& rhs) const
			{
				return IntensityIteratorWrapper(IteratorT::operator+(rhs));
			}
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
	image below and the documentation of #Method for further explanation.

	@image html MorphologicalFilter_all.png

	@note The class #MorphologicalFilter is designed for uniformly spaced raw
	data.

	@ingroup SignalProcessing

	*/
	class /* OPENMS_DLLAPI */ MorphologicalFilter
		:	public ProgressLogger
	{
	public:

		/// Constructor
		MorphologicalFilter ()
		{
		}

		/// Destructor
		virtual ~MorphologicalFilter()
		{
		}

		///@name Available methods
		//@{

		/** @brief List of available morphological filtering methods

		@note This list has to be in sync with #method_names.
		*/
		enum Method
		{
			IDENTITY, ///< identity, uses std::copy algorithm
			EROSION, ///< erosion
			DILATION, ///< dilation
			OPENING, ///<  opening, i.e., erosion followed by dilation
			CLOSING, ///< closing, i.e., dilation followed by erosion
			GRADIENT, ///< gradient, i.e., dilation minus erosion
			TOPHAT, ///< "white" top hat, i.e., signal minus opening
			BOTHAT, ///< "black" top hat, i.e., signal minus closing
			EROSION_SIMPLE, ///< erosion, simple implementation
			DILATION_SIMPLE, ///< dilation, simple implementation
			NUMBER_OF_METHODS
			/**< @brief What the name says.  This enum value must be the last one.
			The list of enumerated values must be a consecutive range starting from
			0.
			*/
		};
		
		/**@brief List of names of available morphological filtering methods.
		These are: "identity", "erosion", "dilation", "opening", "closing", "gradient",
		"tophat", "bothat", "erosion_simple", "dilation_simple".

		@note This list has to be in sync with enum #Method.  
		*/
		static const char* method_names[MorphologicalFilter::NUMBER_OF_METHODS] ;

		/**@brief Convert a method name to its corresponding enum value.

		@param rhs a method name written in characters, e.g. "erosion", see #method_names.
		@return the corresponding enum value, e.g. EROSION, see #Method.

		@exception Exception::IllegalArgument The given argument does not name a method.
		*/
		static Method method( const std::string &rhs )
		{
			Method result = (Method)( std::find( method_names, method_names + NUMBER_OF_METHODS, rhs) - method_names );
			if ( result == NUMBER_OF_METHODS )
			{
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,rhs);
			}
			return result;
		}

		//@}

		/** @brief Applies the morphological filtering operation to an iterator
		range. Input and output range must be valid, i.e. allocated before.
		InputIterator must be a random access iterator type.

		@param method specifies the morphological filtering operation to be applied, see #Method, #method_names, and #method()
		@param struc_size_in_datapoints specifies the size of the 'structuring element' in units of data points.  Should  be an odd number.
		@param input_begin the begin of the input range
		@param input_end  the end of the input range
		@param output_begin the begin of the output range

		@exception Exception::IllegalArgument The given method is not one of the values defined in #Method.
		*/
		template < typename InputIterator, typename OutputIterator >
			void filterRange( Method method, Int struc_size_in_datapoints,
												InputIterator input_begin, InputIterator input_end, OutputIterator output_begin
											)
		{
			typedef typename InputIterator::value_type ValueType;
			// the buffer is static only to avoid reallocation
			static std::vector<ValueType> buffer;
			const UInt size = input_end - input_begin;
			switch ( method )
			{
			case IDENTITY:
				std::copy(input_begin,input_end,output_begin);
				break;
			case EROSION:
				applyErosion_(struc_size_in_datapoints,input_begin,input_end,output_begin);
				break;
			case DILATION:
				applyDilation_(struc_size_in_datapoints,input_begin,input_end,output_begin);
				break;
			case OPENING:
				{
					if ( buffer.size() < size ) buffer.resize(size);
					applyErosion_(struc_size_in_datapoints,input_begin,input_end,buffer.begin());
					applyDilation_(struc_size_in_datapoints,buffer.begin(),buffer.begin()+size,output_begin);
				}
				break;
			case CLOSING:
				{
					if ( buffer.size() < size ) buffer.resize(size);
 					applyDilation_(struc_size_in_datapoints,input_begin,input_end,buffer.begin());
					applyErosion_(struc_size_in_datapoints,buffer.begin(),buffer.begin()+size,output_begin);
				}
				break;
			case GRADIENT:
				{
					if ( buffer.size() < size ) buffer.resize(size);
					applyErosion_(struc_size_in_datapoints,input_begin,input_end,buffer.begin());
					applyDilation_(struc_size_in_datapoints,input_begin,input_end,output_begin);
					for ( UInt i = 0; i < size; ++i ) output_begin[i] -= buffer[i];
				}
				break;
			case TOPHAT:
				{
					if ( buffer.size() < size ) buffer.resize(size);
					applyErosion_(struc_size_in_datapoints,input_begin,input_end,buffer.begin());
					applyDilation_(struc_size_in_datapoints,buffer.begin(),buffer.begin()+size,output_begin);
					for ( UInt i = 0; i < size; ++i ) output_begin[i] = input_begin[i] - output_begin[i];
				}
				break;
			case BOTHAT:
				{
					if ( buffer.size() < size ) buffer.resize(size);
 					applyDilation_(struc_size_in_datapoints,input_begin,input_end,buffer.begin());
					applyErosion_(struc_size_in_datapoints,buffer.begin(),buffer.begin()+size,output_begin);
					for ( UInt i = 0; i < size; ++i ) output_begin[i] = input_begin[i] - output_begin[i];
				}
				break;
			case EROSION_SIMPLE:
				applyErosionSimple_(struc_size_in_datapoints,input_begin,input_end,output_begin);
				break;
			case DILATION_SIMPLE:
				applyDilationSimple_(struc_size_in_datapoints,input_begin,input_end,output_begin);
				break;
			default:
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(method));
				break;
			}
			return;
		}

		/** @brief Applies the morphological filtering operation to a range of
		Peak1D.  The filter is applied to the intensities; the positions are
		copied.  Input and output range must be valid, i.e. allocated before.
		InputIterator must be a random access iterator type.

		@param method specifies the morphological filtering operation to be
		applied, see #Method, #method_names, and #method()

		@param struc_size specifies the size of the 'structuring element'.  This
		can be given in Thomson or by the number of data points, see
		#is_struc_size_in_thomson.

		@param is_struc_size_in_thomson If true, then struc_size is interpreted as
		the width of the structuring element in units of Thomson.  If false, then
		struc_size is taken as the width of the 'structuring elemement' in units
		of data points.  (Further details are given below.)

		@param input_begin the begin of the input range

		@param input_end  the end of the input range

		@param output_begin the begin of the output range

		If is_struc_size_in_thomson is true, then the number of data points for
		the structuring element is computed as follows:
		<ul>
		<li>The data points are assumed to be uniformly spaced.  We compute the
		average spacing from the position of the first and the last peak and the
		total number of peaks in the input range.</li>
		<li>The number of data points in the structuring element is computed
		from struc_size and the average spacing, and rounded up to an odd
		number.</li>
		</ul>
		
		For uniformly spaced data you can compute the number of data points in the
		structuring element as follows:
		@code
		UInt struc_size_datapoints = UInt ( ceil ( struc_size / spacing ) );
		if ( !Math::isOdd(struc_size_datapoints) ) ++struc_size_datapoints;
		@endcode

		@exception Exception::IllegalArgument The given method is not one of the
		values defined in #Method.

		*/
		template < typename InputIterator, typename OutputIterator >
			void filterPeak1DRange( Method method, DoubleReal struc_size, bool is_struc_size_in_thomson,
															InputIterator input_begin, InputIterator input_end, OutputIterator output_begin
														)
		{
			if ( input_end - input_begin <= 1 ) return;
			UInt struc_size_in_datapoints;
			if ( is_struc_size_in_thomson )
			{
				struc_size_in_datapoints =
					UInt (
								ceil (
											struc_size
											*
											DoubleReal( (input_end-input_begin) - 1 )
											/
											( (input_end+(-1))->getPos() - input_begin->getPos() )
										 )
							 );
				// the number has to be odd
				if ( ! Math::isOdd(struc_size_in_datapoints) ) ++struc_size_in_datapoints;
			}
			else
			{
				struc_size_in_datapoints = struc_size;
			}
			filterRange( method, struc_size_in_datapoints,
									 Internal::intensityIteratorWrapper(input_begin),
									 Internal::intensityIteratorWrapper(input_end),
									 Internal::intensityIteratorWrapper(output_begin)
								 );
			for ( ; input_begin != input_end; ++input_begin, ++output_begin )
			{
				output_begin->setPos(input_begin->getPos());
			}
			return;
		}

		/** @brief Applies the morphological filtering operation to a spectrum.

		This is designed for spectra of type #MSSpectrum which contain peaks of
		type #Peak1D.  (But it may work for other types as well).  The
		SpectrumSettings are copied.
		*/
		template <typename InputPeakContainer, typename OutputPeakContainer >
			void filterMSSpectrum( Method method, DoubleReal struc_size, bool is_struc_size_in_thomson, 
														 const InputPeakContainer& input_peak_container, OutputPeakContainer& output_peak_container
													 )
		{
			static_cast<SpectrumSettings&>(output_peak_container) = input_peak_container;
			output_peak_container.setType(SpectrumSettings::RAWDATA);
			output_peak_container.resize(input_peak_container.size());
			filterPeak1DRange( method, struc_size, is_struc_size_in_thomson,
												 input_peak_container.begin(),
												 input_peak_container.end(),
												 output_peak_container.begin()
											 );
			return;
		}


		/** @brief Applies the morphological filtering operation to an
		MSExperiment.  The size of the structuring element is computed for each
		spectrum individually, see #filterPeak1DRange() and #filterMSSpectrum().
		*/
		template <typename PeakType, typename AllocType>
			void filterMSExperiment( Method method, DoubleReal struc_size, bool is_struc_size_in_thomson,
															 MSExperiment<PeakType, AllocType>& ms_exp
														 )
		{
			startProgress(0,ms_exp.size(),"filtering baseline");
			for ( UInt i = 0; i < ms_exp.size(); ++i )
			{
				typename MSExperiment<PeakType>::SpectrumType spectrum;
				filterMSSpectrum( method, struc_size, is_struc_size_in_thomson, ms_exp[i], spectrum );
				ms_exp[i].swap(spectrum);
				setProgress(i);
			}
			endProgress();
			return;
		}

	protected:

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

	};
	
}// namespace OpenMS

#endif
