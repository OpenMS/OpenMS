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

		/**	@brief Internal class for asymmetric distributions
				
		Internal class for asymmetric distributions
		used for consistency with BasisStatistic class
		
		@todo Various performance improvements are possible, see code. (Clemens)
		@todo Write test and improve documentation. (Clemens)

		*/
		template < typename Real = double >
		class AsymmetricStatistics : public BasicStatistics<Real>
		{
			typedef BasicStatistics<Real> Base;
			typedef typename Base::RealType RealType;

			Base::clear;
			Base::sum_;
			Base::mean_;
			Base::variance_;
			
		 public:
			AsymmetricStatistics()
				: BasicStatistics<>(),
					variance1_(1.0),
					variance2_(1.0)
			{}

			RealType variance1() const throw(){ return variance1_; }
			RealType variance2() const throw(){ return variance2_; }

			template < typename ProbabilityIterator, typename CoordinateIterator >
			void update( ProbabilityIterator const probability_begin,
									 ProbabilityIterator const probability_end,
									 CoordinateIterator  const coordinate_begin)
			{
				clear();
				variance1_ = 0;
				variance2_ = 0;
				ProbabilityIterator prob_iter = probability_begin;
				CoordinateIterator  coord_iter = coordinate_begin;

				for ( ; prob_iter != probability_end; ++prob_iter, ++coord_iter )
				{
					sum_  += *prob_iter;
					mean_ += *prob_iter * *coord_iter;
				}
				mean_ /= sum_;

				RealType sum1(0), sum2(0);
				for ( prob_iter = probability_begin, coord_iter = coordinate_begin;
							prob_iter != probability_end;
							++prob_iter, ++coord_iter
						)
				{
					RealType diff = *coord_iter - mean_;
					diff *= diff;
					variance_ += *prob_iter * diff;

					// TODO Improvements - to be discussed with maintainer - :
					// Why not simply have two cases: < vs. >=  ?
					// Factor 2 can be multiplied after (outside) summation loop.
					// Compute mean_ +/- precision outside this loop.
					const double precision = 0.0001;
					if (*coord_iter < mean_-precision)
					{
						variance1_ += 2* (*prob_iter * diff);
						sum1  += *prob_iter *2;
					}
					else if (*coord_iter > mean_+precision)
					{
						variance2_ += 2* (*prob_iter * diff);
						sum2  += *prob_iter *2;
					}
					else // coord_iter == mean
					{
						variance1_ += *prob_iter * diff;
						sum1  += *prob_iter;
						variance2_ += *prob_iter * diff;
						sum2  += *prob_iter;
					}
				}
				variance1_ /= sum1;
				variance2_ /= sum2;
				variance_ /= sum_;
			}
		 protected:
			RealType variance1_, variance2_;
		};
		
	} // namespace Math 
 
} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_ASYMMETRICSTATISTICS_H
