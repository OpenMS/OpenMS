// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

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

    template <typename RealT = double>
    class AsymmetricStatistics :
      public BasicStatistics<RealT>
    {
      /// The real type and basic statistics specified as template argument.
      typedef BasicStatistics<RealT> Base;
      typedef typename Base::RealType RealType;

      using Base::clear;
      using Base::sum_;
      using Base::mean_;
      using Base::variance_;

public:

      /// Default constructor.
      AsymmetricStatistics() :
        BasicStatistics<>(),
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
      template <typename ProbabilityIterator, typename CoordinateIterator>
      void update(ProbabilityIterator const probability_begin,
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
        for (; prob_iter != probability_end; ++prob_iter, ++coord_iter)
        {
          RealType diff = *coord_iter - mean_;
          RealType diff_squared = diff * diff;

          if (diff_squared > variance_)
          {
            if (*coord_iter < mean_)
            {
              variance1_ += (*prob_iter * diff_squared);
              sum1  += *prob_iter;
            }
            else             // ( *coord_iter > mean_ )
            {
              variance2_ += (*prob_iter * diff_squared);
              sum2  += *prob_iter;
            }
          }
          else
          {
            RealType frac = (diff / stdev + 1.) / 2.;
            RealType prob_frac = frac * *prob_iter;
            variance2_ += prob_frac * diff_squared;
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

  }   // namespace Math

} // namespace OpenMS

