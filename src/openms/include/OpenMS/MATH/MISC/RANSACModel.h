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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{

  namespace Math
  {

    /**
      @brief Generic plug-in template base class using 'Curiously recurring template pattern' (CRTP)
             to allow for arbitrary RANSAC models (e.g. linear or quadratic fits).

      An actual model should derive from this class using CRTP (see RansacModelLinear for an example)
      and implement the *_impl() functions.

    */
    template<class ModelT = int> // dummy default to allow easy access to public typedefs
    class RansacModel
    {
  public:
      typedef std::pair<double, double> DPair;
      typedef std::vector<DPair> DVec;
      typedef DVec::const_iterator DVecIt;
      typedef std::vector<double> ModelParameters;


      /// fit a model and return its parameters
      ModelParameters rm_fit(const DVecIt& begin, const DVecIt& end) const
      {
        return static_cast<const ModelT*>(this)->rm_fit_impl(begin, end);
      }
      
      /**
        @brief Returns the R-squared of the data applied to the model (computed on-the-fly).
      
        Takes as input a standard vector of a standard pair of points in a 2D space.
        
        @param begin Iterator to first pair
        @param end Past-end iterator to last pair
        @return: R-squared value
      */
      double rm_rsq(const DVecIt& begin, const DVecIt& end) const
      {
        return static_cast<const ModelT*>(this)->rm_rsq_impl(begin, end);
      }

      /// calculates the residual sum of squares of the input points according to the model
      double rm_rss(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients) const
      {
        return static_cast<const ModelT*>(this)->rm_rss_impl(begin, end, coefficients);
      }

      /// calculates the squared residual of each input point vs. the model.
      /// All points with an error <= max_threshold, are accepted as inliers and returned.
      DVec rm_inliers(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients, double max_threshold) const
      {
        return static_cast<const ModelT*>(this)->rm_inliers_impl(begin, end, coefficients, max_threshold);
      }
    };

  }


}
