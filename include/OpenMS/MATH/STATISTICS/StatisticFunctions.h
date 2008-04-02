// -*- Mode: C++; tab-width: 2; -*-
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

#include <numeric>
#include <OpenMS/CONCEPT/Types.h>

// #include <iostream> // debugging

#ifndef OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H
#define OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H

namespace OpenMS 
{

	namespace Math
	{

		/**
		@brief Calculates the mean square error for the values in [begin_a, end_a)
		and [begin_b, end_b)

		Calculates the mean square error for the data given by the two iterator
		ranges. If one of the ranges contains a smaller number of values the rest
		of the longer range is omitted.

		@ingroup Math

		*/
		template < typename IteratorType1, typename IteratorType2 >
		inline static
		DoubleReal meanSquareError ( IteratorType1 begin_a, const IteratorType1 end_a,
																 IteratorType2 begin_b, const IteratorType2 end_b
															 )
		{
			Int count = 0;
			DoubleReal error = 0;
			IteratorType1 & it_a = begin_a;
			IteratorType2 & it_b = begin_b;

			while(it_a != end_a && it_b != end_b)
			{
				DoubleReal diff = *it_a - *it_b;
				error += diff * diff;
				++count;
				++it_a;
				++it_b;
			}

			return error / count;

		}

		/**
		@brief Calculates the classification rate for the values in [begin_a,
		end_a) and [begin_b, end_b)

		Calculates the classification rate for the data given by the two iterator ranges. If
		one of the ranges contains a smaller number of values the rest of the longer range is omitted.

		@ingroup Math

		*/
		template < typename IteratorType1, typename IteratorType2 >
		inline static
		Real classificationRate ( IteratorType1 begin_a, const IteratorType1 end_a,
															IteratorType2 begin_b, const IteratorType2 end_b
														)
		{
			Int count = 0;
			DoubleReal error = 0;
			IteratorType1 & it_a = begin_a;
			IteratorType2 & it_b = begin_b;

			if (it_a == end_a || it_b == end_b)
			{
				return 0;
			}

			while(it_a != end_a && it_b != end_b)
			{
				if ((*it_a < 0 && *it_b >= 0)
						|| (*it_a >= 0 && *it_b < 0))
				{
					error += 1;
				}
				++count;
				++it_a;
				++it_b;
			}

			return (count - error) / count;

		}

		/**@brief Calculates the Matthews Correlation Coefficient for the values
		in [begin_a, end_a) and [begin_b, end_b)

		Calculates the Matthews Correlation Coefficient for the data given by the
		two iterator ranges. If one of the ranges contains a smaller number of
		values the rest of the longer range is omitted.  The values in [begin_a,
		end_a) have to be the predicted labels and the values in [begin_b, end_b)
		have to be the real labels.

		Formerly called mcc, renamed for obvious reason.

		@ingroup Math

		*/
		template < typename IteratorType1, typename IteratorType2 >
		inline static  /* mcc */ 
		DoubleReal matthewsCorrelationCoefficient( IteratorType1 begin_a, const IteratorType1 end_a,
																						 IteratorType2 begin_b, const IteratorType2 end_b
																					 )
		{
			IteratorType1 & it_a = begin_a;
			IteratorType2 & it_b = begin_b;
			DoubleReal tp = 0;
			DoubleReal fp = 0;
			DoubleReal tn = 0;
			DoubleReal fn = 0;

			if (it_a == end_a || it_b == end_b)
			{
				return 0;
			}

			while(it_a != end_a && it_b != end_b)
			{
				if (*it_a < 0 && *it_b >= 0)
				{
					++fn;
				}
				else if (*it_a < 0 && *it_b < 0)
				{
					++tn;
				}
				else if (*it_a >= 0 && *it_b >= 0)
				{
					++tp;
				}
				else if (*it_a >= 0 && *it_b < 0)
				{
					++fp;
				}

				++it_a;
				++it_b;
			}

			return ((tp * tn - fp * fn) / 
							sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)));

		}

		/**
		@brief calculates the pearson correlation coefficient for the values in
		[begin_a, end_a) and [begin_b, end_b)

		Calculates the linear correlation coefficient for the data given by the
		two iterator ranges.

		If the iterator ranges are not of the same length or empty an exception is
		thrown.
				
		If one of the ranges contains only the same values 'nan' is returned.

		@ingroup Math

		*/
		template < typename IteratorType1, typename IteratorType2 >
		inline static
		DoubleReal pearsonCorrelationCoefficient ( IteratorType1 begin_a, IteratorType1 end_a,
																						 IteratorType2 begin_b, IteratorType2 end_b
																					 )
			throw (Exception::InvalidRange)
		{
			UInt count = end_a-begin_a;
			//no data or different lengths
			if (count==0 || end_a-begin_a!=end_b-begin_b)
			{
				throw Exception::InvalidRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
				
			//calculate average
			DoubleReal avg_a = std::accumulate(begin_a,end_a,0.0) / count;
			DoubleReal avg_b = std::accumulate(begin_b,end_b,0.0) / count;

			DoubleReal numerator = 0;
			DoubleReal denominator_a = 0;
			DoubleReal denominator_b = 0;
			while(begin_a != end_a)
			{
				DoubleReal temp_a = *begin_a - avg_a;
				DoubleReal temp_b = *begin_b - avg_b;
				numerator += (temp_a * temp_b);
				denominator_a += (temp_a * temp_a);
				denominator_b += (temp_b * temp_b);
				++begin_a;
				++begin_b;
			}
				
			return numerator / sqrt(denominator_a * denominator_b);
		}

	} // namespace Math
} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H
