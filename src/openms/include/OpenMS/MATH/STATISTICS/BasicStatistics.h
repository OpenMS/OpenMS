// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <vector>
#include <ostream>
#include <cmath>
#include <numeric>

namespace OpenMS
{
  namespace Math
  {

    /**
         @brief Calculates some basic statistical parameters of a distribution:
         sum, mean, variance, and provides the normal approximation.

         The intended usage is as follows:
         - <i>create</i> an instance
         - <i>set</i> the basic statistical parameters by either
                - calling one of the update() member functions, or
                - using the set... methods
                .
         - <i>do</i> something with the basic statistical parameters, e.g.
                - using the get... methods, or
                - obtain samples from a normal approximation with these parameters
                - whatever member function you might want to add to this class ;-)
         .

         @ingroup Math
    */
    template <typename RealT = double>
    class
    BasicStatistics
    {

public:

      /// The real type specified as template argument.
      typedef RealT RealType;

      typedef std::vector<RealType> probability_container;
      typedef std::vector<RealType> coordinate_container;

      /// Default constructor.
      BasicStatistics() :
        mean_(0),
        variance_(0),
        sum_(0)
      {}

      /// Copy constructor.
      BasicStatistics(BasicStatistics const & arg) :
        mean_(arg.mean_),
        variance_(arg.variance_),
        sum_(arg.sum_)
      {}

      /// Assignment.
      BasicStatistics & operator=(BasicStatistics const & arg)
      {
        mean_      = arg.mean_;
        variance_  = arg.variance_;
        sum_       = arg.sum_;
        return *this;
      }

      /// Set sum, mean, and variance to zero.
      void clear()
      {
        mean_ = 0;
        variance_ = 0;
        sum_ = 0;
      }

      /// This does the actual calculation.
      /** You can call this as often as you like, using different input vectors. */
      template <typename ProbabilityIterator>
      void update(ProbabilityIterator probability_begin,
                  ProbabilityIterator const probability_end
                  )
      {
        clear();
        unsigned pos = 0;
        ProbabilityIterator iter = probability_begin;

        for (; iter != probability_end; ++iter, ++pos)
        {
          sum_  += *iter;
          mean_ += *iter * pos;
        }
        mean_ /= sum_;

        for (iter = probability_begin, pos = 0; iter != probability_end; ++iter, ++pos)
        {
          RealType diff = RealType(pos) - mean_;
          diff *= diff;
          variance_ += *iter * diff;
        }
        variance_ /= sum_;

        if (sum_ == 0 && (std::isnan(mean_) || std::isinf(mean_)) )
        {
          mean_ = 0;
          variance_ = 0;
        }
      }

      /// This does the actual calculation.
      /** You can call this as often as you like, using different input vectors. */
      template <typename ProbabilityIterator, typename CoordinateIterator>
      void update(ProbabilityIterator const probability_begin,
                  ProbabilityIterator const probability_end,
                  CoordinateIterator  const coordinate_begin
                  )
      {
        clear();
        ProbabilityIterator prob_iter = probability_begin;
        CoordinateIterator  coord_iter = coordinate_begin;

        for (; prob_iter != probability_end; ++prob_iter, ++coord_iter)
        {
          sum_  += *prob_iter;
          mean_ += *prob_iter * *coord_iter;
        }
        mean_ /= sum_;

        for (prob_iter = probability_begin, coord_iter = coordinate_begin;
             prob_iter != probability_end;
             ++prob_iter, ++coord_iter
             )
        {
          RealType diff = *coord_iter - mean_;
          diff *= diff;
          variance_ += *prob_iter * diff;
        }
        variance_ /= sum_;
        return;
      }

      /// Returns the mean.
      RealType mean() const { return mean_; }
      void setMean(RealType const & mean) { mean_ = mean; }

      /// Returns the variance.
      RealType variance() const { return variance_; }
      void setVariance(RealType const & variance) { variance_ = variance; }

      /// Returns the sum.
      RealType sum() const { return sum_; }
      void setSum(RealType const & sum) { sum_ = sum; }


      /**@brief Returns the density of the normal approximation at point,
           multiplied by sqrt( 2 * pi ).  This saves a division operation compared
           to normalDensity()
      */
      RealType normalDensity_sqrt2pi(RealType coordinate) const
      {
        coordinate -= mean();
        coordinate *= coordinate;
        return exp(-coordinate / RealType(2) / variance());
      }

      /// Returns sqrt( 2 * pi ), which is useful to normalize the result of normalDensity_sqrt2pi().
      static RealType sqrt2pi() { return 2.50662827463100050240; }

      /**@brief See normalDensity_sqrt2pi().  Returns the density of the normal
           distribution at point.
      */
      inline RealType normalDensity(RealType const coordinate) const
      {
        return normalDensity_sqrt2pi(coordinate) / sqrt2pi();
      }

      /**@brief The argument \c probability is filled with values according to
           the normal approximation.  Its \c size() is not changed.  The
           approximation takes place at coordinate positions 0, 1, ..., size()-1.
      */
      void normalApproximation(probability_container & probability)
      {
        normalApproximationHelper_(probability, probability.size());
        return;
      }

      /** The argument \c probability is filled with values according to the
              normal approximation.  Its size() is set to \c size.  The
              approximation takes place at coordinate positions 0, 1, ..., size-1.
      */
      void normalApproximation(probability_container & probability,
                               typename probability_container::size_type const size
                               )
      {
        probability.resize(size);
        normalApproximationHelper_(probability, size);
        return;
      }

      /**@brief The argument probability is filled with values according to the
           normal approximation.  The second argument coordinate contains the
           positions where the approximation takes place.  probability.size() is
           set to coordinate.size().
      */
      void normalApproximation(probability_container & probability,
                               coordinate_container  const & coordinate
                               )
      {
        probability.resize(coordinate.size());
        normalApproximationHelper_(probability, coordinate);
        return;
      }

      /**@brief A convenient overload for debugging purposes.

      @relatesalso BasicStatistics
      */
      friend std::ostream & operator<<(std::ostream & os, BasicStatistics & arg)
      {
        os << "BasicStatistics:  mean=" << arg.mean() << "  variance=" << arg.variance() << "  sum=" << arg.sum();
        return os;
      }

protected:

      /// @name Protected Members
      //@{

      RealType mean_;
      RealType variance_;
      RealType sum_;

private:
      //@}

      /// @name Private Methods
      //@{

      void normalApproximationHelper_(probability_container & probability,
                                      typename probability_container::size_type const size
                                      )
      {
        RealType gaussSum = 0;
        typename coordinate_container::size_type i;

        // precondition size == probability.size() is guaranteed by wrappers.
        for (i = 0; i < size; ++i)
        {
          gaussSum += normalDensity_sqrt2pi(RealType(i));
        }

        for (i = 0; i < size; ++i)
        {
          probability[i] = normalDensity_sqrt2pi(RealType(i)) / gaussSum * sum();
        }
        return;
      }

      void normalApproximationHelper_(probability_container & probability,
                                      coordinate_container const & coordinate
                                      )
      {
        RealType gaussSum = 0;
        typename coordinate_container::size_type i;
        typename coordinate_container::size_type const size = coordinate.size();

        for (i = 0; i < size; ++i)
        {
          gaussSum += normalDensity_sqrt2pi(coordinate[i]);
        }

        for (i = 0; i < size; ++i)
        {
          probability[i] = normalDensity_sqrt2pi(coordinate[i]) / gaussSum * sum();
        }
        return;
      }

      //@}

    };

  }   // namespace Math

} // namespace OpenMS
