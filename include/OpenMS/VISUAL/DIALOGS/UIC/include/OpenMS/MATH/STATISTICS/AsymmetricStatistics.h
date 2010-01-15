// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_MATH_STATISTICS_ASYMMETRICSTATISTICS_H
#define OPENMS_MATH_STATISTICS_ASYMMETRICSTATISTICS_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

#include <vector>
#include <ostream>
#include <cmath>

namespace OpenMS
{
  namespace Math
	{

		/**	
			@brief Internal class for asymmetric distributions
				
			Internal class for asymmetric distributions
			used for consistency with BasisStatistic class
			
		*/

		template < typename Real = DoubleReal > class AsymmetricStatistics :
		 	public BasicStatistics<Real>
		{
			/// The real type and basic statistics specified as template argument.
			typedef BasicStatistics<Real> Base;
			typedef typename Base::RealType RealType;

			Base::clear;
			Base::sum_;
			Base::mean_;
			Base::variance_;
			
		 public:

			/// Default constructor.
			AsymmetricStatistics()
				: BasicStatistics<>(),
					variance1_(0),
					variance2_(0)
			{}

			/// "variance to the left hand side"
			RealType variance1() const
		 	{
			 	return variance1_;
			}
			
			/// "variance to the right hand side"
			RealType variance2() const
			{
			 	return variance2_;
			}

			/// You can call this as often as you like, using different input vectors.
			template < typename ProbabilityIterator, typename CoordinateIterator > void update( ProbabilityIterator const probability_begin,
					ProbabilityIterator const probability_end,
					CoordinateIterator  const coordinate_begin)
			{
				// reuse...
				Base::update(probability_begin, probability_end, coordinate_begin);
				
				const RealType stdev = std::sqrt(variance_);	
				
				RealType sum1 = 0;
				RealType sum2 = 0;
				variance1_ = 0;
				variance2_ = 0;
				ProbabilityIterator prob_iter = probability_begin;
				CoordinateIterator  coord_iter = coordinate_begin;
				for ( ; prob_iter != probability_end; ++prob_iter, ++coord_iter )
				{
					RealType diff = *coord_iter - mean_;
					RealType diff_squared = diff * diff;

					if (diff_squared > variance_)
					{
						if (	*coord_iter < mean_ )
						{
							variance1_ += (*prob_iter * diff_squared);
							sum1  += *prob_iter;
						}
						else // ( *coord_iter > mean_ )
						{
							variance2_ += (*prob_iter * diff_squared);
							sum2  += *prob_iter;
						}
					}
					else
					{
						RealType frac = ( diff / stdev + 1. ) / 2.;
						RealType prob_frac = frac * *prob_iter;
						variance2_ += prob_frac	* diff_squared;
						sum2  += prob_frac;
						prob_frac = *prob_iter * (1. - frac);
						variance1_ += prob_frac * diff_squared;
						sum1  += prob_frac;
					}
				}
				variance1_ /= sum1;
				variance2_ /= sum2;
				return;
			}

		 protected:
			/// @name Protected Members
			RealType variance1_, variance2_;
		};

	} // namespace Math 

} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_ASYMMETRICSTATISTICS_H
