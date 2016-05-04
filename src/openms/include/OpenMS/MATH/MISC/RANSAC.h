// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#ifndef OPENMS_MATH_MISC_RANSAC_H
#define OPENMS_MATH_MISC_RANSAC_H

#include <OpenMS/config.h>

#include <OpenMS/MATH/MISC/RANSACModel.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/MATH/MISC/RANSACModelLinear.h>

#include <algorithm>    // std::random_shuffle
#include <limits>       // std::numeric_limits
#include <vector>       // std::vector

namespace OpenMS
{

  namespace Math
  {

    /**
      @brief This class provides a generic implementation of the RANSAC
       outlier detection algorithm. Is implemented and tested after the
       SciPy reference: http://wiki.scipy.org/Cookbook/RANSAC
    */
    template<typename TModelType = RansacModelLinear>
    class RANSAC
    {
public:

      /**
        @brief This function provides a generic implementation of the RANSAC
         outlier detection algorithm. Is implemented and tested after the
         SciPy reference: http://wiki.scipy.org/Cookbook/RANSAC

        @param pairs Input data (paired data of type <dim1, dim2>)
        @param n The minimum number of data points required to fit the model
        @param k The maximum number of iterations allowed in the algorithm 
        @param t Threshold value for determining when a data point fits a
         model. Corresponds to the maximal squared deviation in units of the
         _second_ dimension (dim2).
        @param d The number of close data values (according to 't') required to assert that a model fits well to data
        @param rng Custom RNG function (useful for testing with fixed seeds)

        @return A vector of pairs fitting the model well; data will be unsorted
      */
  
      static std::vector<std::pair<double, double> > ransac(
          const std::vector<std::pair<double, double> >& pairs, 
          size_t n, 
          size_t k, 
          double t, 
          size_t d, 
          int (*rng)(int) = NULL)
      {
        TModelType model;

        // implementation of the RANSAC algorithm according to http://wiki.scipy.org/Cookbook/RANSAC.

        if (pairs.size() <= n)
        {
          throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                        String("RANSAC: Number of pairs (") + String(pairs.size()) + ") must be larger than number of initial points (n=" + String(n) + ").");
        }

        std::vector< std::pair<double, double> > alsoinliers, betterdata, bestdata;

        double besterror = std::numeric_limits<double>::max();
    #ifdef DEBUG_MRMRTNORMALIZER
        std::pair<double, double > bestcoeff;
        double betterrsq = 0;
        double bestrsq = 0;
    #endif

        // mutable data. will be shuffled in every iteration
        std::vector<std::pair<double, double> > pairs_shuffled = pairs;

        for (size_t ransac_int=0; ransac_int<k; ransac_int++)
        {
          // check if the model already includes all points
          if (bestdata.size() == pairs.size()) break;

          if (rng != NULL)
          { // use portable RNG in test mode
            std::random_shuffle(pairs_shuffled.begin(), pairs_shuffled.end(), rng);
          } else {
            std::random_shuffle(pairs_shuffled.begin(), pairs_shuffled.end());
          }

          // test 'maybeinliers'
          typename TModelType::ModelParameters coeff = model.rm_fit(pairs_shuffled.begin(), pairs_shuffled.begin()+n);
          // apply model to remaining data; pick inliers
          alsoinliers = model.rm_inliers(pairs_shuffled.begin()+n, pairs_shuffled.end(), coeff, t);
          // ... and add data
          if (alsoinliers.size() > d)
          {
            betterdata.clear();
            std::copy( pairs_shuffled.begin(), pairs_shuffled.begin()+n, back_inserter(betterdata) );
            betterdata.insert( betterdata.end(), alsoinliers.begin(), alsoinliers.end() );
            typename TModelType::ModelParameters bettercoeff = model.rm_fit(betterdata.begin(), betterdata.end());
            double bettererror = model.rm_rss(betterdata.begin(), betterdata.end(), bettercoeff);
    #ifdef DEBUG_MRMRTNORMALIZER
            betterrsq = model.rm_rsq(betterdata);
    #endif

            if (bettererror < besterror)
            {
              besterror = bettererror;
              bestdata = betterdata;
    #ifdef DEBUG_MRMRTNORMALIZER
              bestcoeff = bettercoeff;
              bestrsq = betterrsq;
              std::cout << "RANSAC " << ransac_int << ": Points: " << betterdata.size() << " RSQ: " << bestrsq << " Error: " << besterror << " c0: " << bestcoeff.first << " c1: " << bestcoeff.second << std::endl;
    #endif
            }
          }
        }

    #ifdef DEBUG_MRMRTNORMALIZER
        std::cout << "=======STARTPOINTS=======" << std::endl;
        for (std::vector<std::pair<double, double> >::iterator it = bestdata.begin(); it != bestdata.end(); ++it)
        {
          std::cout << it->first << "\t" << it->second << std::endl;
        }
        std::cout << "=======ENDPOINTS=======" << std::endl;
    #endif

        return(bestdata);
      } // ransac()

    }; // class
  
  } // namespace Math


} // namespace OpenMS
#endif // OPENMS_MATH_MISC_RANSAC_H
