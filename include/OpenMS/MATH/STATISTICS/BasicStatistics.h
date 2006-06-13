// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: BasicStatistics.h,v 1.6 2006/03/28 16:19:59 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_STATISTICS_BASICSTATISTICS_H
#define OPENMS_MATH_STATISTICS_BASICSTATISTICS_H

#include <vector>
#include <ostream>
#include <cmath>

namespace OpenMS
{
  
  /**
  	@brief Calculates some basic statistical parameters of a distribution: sum, mean, variance, and provides the normal approximation.

  	@todo add to Math namespace (Clemens)
  	
  	@ingroup Math
  */
  template < typename Coordinate = double, typename Probability = Coordinate >
  class
  BasicStatistics
  {

   public:

    /// The type used for probabilities, intensities, etc.
    typedef Probability probability_type;
		typedef probability_type ProbabilityType;

    /// The type used for coordinate, positions, etc.
    typedef Coordinate coordinate_type;
		typedef coordinate_type CoordinateType;

		typedef std::vector < probability_type > probability_container;
		typedef std::vector < coordinate_type > coordinate_container;

    /// Default constructor.
    BasicStatistics ()
      : mean_(0),
        variance_(0),
        sum_(0)
    {}

    /// Copy constructor.
    BasicStatistics ( BasicStatistics const & arg )
      : mean_(arg.mean_),
        variance_(arg.variance_),
        sum_(arg.sum_)
    {}

		/// Assignment.
		BasicStatistics & operator = ( BasicStatistics const & arg )
		{
			mean_      = arg.mean_;
			variance_  = arg.variance_;
			sum_       = arg.sum_;
			return *this;
		}

    /// Use this constructor if you want to update() immediately.
		template < typename ProbabilityIterator >
    BasicStatistics ( ProbabilityIterator const & probability_begin,
											ProbabilityIterator const & probability_end
										)
    {
      update ( probability_begin, probability_end );
      return;
    }

    /// Use this constructor if you want to update() immediately.
		template < typename ProbabilityIterator, typename CoordinateIterator >
    BasicStatistics ( ProbabilityIterator const & probability_begin,
											ProbabilityIterator const & probability_end,
                      CoordinateIterator  const & coordinate_begin
										)
    {
      update ( probability_begin, probability_end, coordinate_begin );
    }

		/// Set sum, mean, and variance to zero.
		void clear ()
		{
			mean_ = 0;
      variance_ = 0;
      sum_ = 0;
		}

    /// This does the actual calculation.
    /** You can call this as often as you like, using different input vectors. */
		template < typename ProbabilityIterator >
    void update ( ProbabilityIterator probability_begin,
									ProbabilityIterator const probability_end
								)
    {
			clear();
			unsigned pos = 0;
			ProbabilityIterator iter = probability_begin;
			
      for ( ; iter != probability_end; ++iter, ++pos )
			{
        sum_  += *iter;
        mean_ += *iter * pos;
      }
      mean_ /= sum_;
      
      for ( iter = probability_begin, pos = 0; iter != probability_end; ++iter, ++pos )
			{
        coordinate_type diff = coordinate_type(pos) - mean_;
        diff *= diff;
        variance_ += *iter * diff;
      }
      variance_ /= sum_;
    }

    /// This does the actual calculation.
    /** You can call this as often as you like, using different input vectors. */
		template < typename ProbabilityIterator, typename CoordinateIterator >
    void update ( ProbabilityIterator const probability_begin,
									ProbabilityIterator const probability_end,
									CoordinateIterator  const coordinate_begin
								)
    {
			clear();
			ProbabilityIterator prob_iter = probability_begin;
			CoordinateIterator  coord_iter = coordinate_begin;

      for ( ; prob_iter != probability_end; ++prob_iter, ++coord_iter )
			{
        sum_  += *prob_iter;
        mean_ += *prob_iter * *coord_iter;
      }
      mean_ /= sum_;

      for ( prob_iter = probability_begin, coord_iter = coordinate_begin;
						prob_iter != probability_end;
						++prob_iter, ++coord_iter
					)
			{
        coordinate_type diff = *coord_iter - mean_;
        diff *= diff;
        variance_ += *prob_iter * diff;
      }
      variance_ /= sum_;
      return;
    }

    /// Returns the mean.
    coordinate_type mean ()     const throw() { return mean_; }
		void setMean ( coordinate_type const & mean ) throw() { mean_ = mean; }

    /// Returns the variance.
    coordinate_type variance() const throw() { return variance_; }
		void setVariance ( coordinate_type const & variance ) throw() { variance_ = variance; }

    /// Returns the sum.
    probability_type sum()     const throw() { return sum_; }
		void setSum ( probability_type const & sum ) throw() { sum_ = sum; }


		/**@brief Returns the density of the normal approximation at point,
			 multiplied by sqrt( 2 * pi ).  This saves a division operation compared
			 to normalDensity()
		*/
		probability_type normalDensity_sqrt2pi ( coordinate_type coordinate ) const throw()
		{
			coordinate -= mean();
			coordinate *= coordinate;
			return exp ( - coordinate / coordinate_type(2) / variance() );
		}

		/// Returns sqrt( 2 * pi ), which is useful to normalize the result of normalDensity_sqrt2pi().
		static probability_type sqrt2pi () throw() { return 2.50662827463100050240; }

		/**@brief See normalDensity_sqrt2pi().  Returns the density of the normal
			 distribution at point.
		*/
		inline probability_type normalDensity ( coordinate_type & coordinate ) const throw()
		{
			return normalDensity_sqrt2pi ( coordinate ) / sqrt2pi() ;
		}

    /**@brief The argument \c probability is filled with values according to
       the normal approximation.  Its \c size() is not changed.  The
       approximation takes place at coordinate positions 0, 1, ..., size()-1.
    */
    void normalApproximation ( probability_container & probability )
    {
      normalApproximation_internal ( probability, probability.size() );
      return;
    }

    /** The argument \c probability is filled with values according to the
        normal approximation.  Its size() is set to \c _size.  The
        approximation takes place at coordinate positions 0, 1, ..., _size-1.
    */
    void normalApproximation ( probability_container & probability,
                               typename probability_container::size_type const _size
                             )
    {
      probability.resize ( _size );
      normalApproximation_internal ( probability, _size );
      return;
    }

    /**@brief The argument probability is filled with values according to the
       normal approximation.  The second argument coordinate contains the
       positions where the approximation takes place.  probability.size() is
       set to coordinate.size().
    */
    void normalApproximation ( probability_container & probability,
															 coordinate_container  const & coordinate
                             )
    { 
      probability.resize ( coordinate.size() );
      normalApproximation_internal ( probability, coordinate );
      return;
    }

    /// A convenient overload for debugging purposes.
    friend std::ostream & operator << ( std::ostream & os, BasicStatistics & arg )
    {
      os << "BasicStatistics:  mean=" << arg.mean() << "  variance=" << arg.variance() << "  sum=" << arg.sum();
      return os;
    }

   protected:

    /// @name Protected Members
    //@{ 

    coordinate_type mean_, variance_;
    probability_type sum_;

		private:
    //@}

    /// @name Private Methods
    //@{ 

    void normalApproximation_internal ( probability_container & probability,
                                        typename probability_container::size_type const _size
																			)
    {
      probability_type gaussSum = 0;
      typename coordinate_container::size_type i;
      
      // precondition _size == probability.size() is guaranteed by wrappers.
      for ( i = 0; i < _size; ++i ) {
				gaussSum += normalDensity_sqrt2pi ( i );
      }

      for ( i = 0; i < _size; ++i ) {
				probability [ i ] = normalDensity_sqrt2pi ( i ) / gaussSum * sum();
      }
      return;
    }

    void normalApproximation_internal ( probability_container & probability,
																				coordinate_container const & coordinate
																			)
    {
      probability_type gaussSum = 0;
      typename coordinate_container::size_type i;
      typename coordinate_container::size_type const size = coordinate.size();

      for ( i = 0; i < size; ++i ) {
				gaussSum += normalDensity_sqrt2pi ( coordinate[i] );
      }

      for ( i = 0; i < size; ++i ) {
				probability [ i ] = normalDensity_sqrt2pi ( coordinate[i] ) / gaussSum * sum();
      }
      return;
    }

    //@}

  };
  
} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_BASICSTATISTICS_H
