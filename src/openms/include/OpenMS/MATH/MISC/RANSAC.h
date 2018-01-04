// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <sstream>      // stringstream

namespace OpenMS
{

  namespace Math
  {
    /**
      @brief A simple struct to carry all the parameters required for a RANSAC run.
    */
    struct RANSACParam
    {
      /// Default constructor
      RANSACParam()
        : n(0), k(0), t(0), d(0), relative_d(false), rng(nullptr)
        {
        }
      /// Full constructor
      RANSACParam(size_t p_n, size_t p_k, double p_t, size_t p_d, bool p_relative_d = false, int (*p_rng)(int) = nullptr)
        : n(p_n), k(p_k), t(p_t), d(p_d), relative_d(p_relative_d), rng(p_rng)
      {
        if (relative_d)
        {
          if (d >= 100) throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("RANSAC: Relative 'd' >= 100% given. Use a lower value; the more outliers you expect, the lower it should be."));
        }
      }

      std::string toString() const
      {
        std::stringstream r;
        r << "RANSAC param:\n  n: " << n << "\n  k: " << k << " iterations\n  t: " << t << " threshold\n  d: " << d << " inliers\n\n";
        return r.str();
      }

      size_t n; ///< data points: The minimum number of data points required to fit the model
      size_t k; ///< iterations: The maximum number of iterations allowed in the algorithm 
      double t; ///< Threshold value: for determining when a data point fits a model. Corresponds to the maximal squared deviation in units of the _second_ dimension (dim2).
      size_t d; ///< The number of close data values (according to 't') required to assert that a model fits well to data
      bool relative_d; ///< Should 'd' be interpreted as percentages (0-100) of data input size.
      int (*rng)(int); ///< Optional RNG function (useful for testing with fixed seeds)
    };

    /**
      @brief This class provides a generic implementation of the RANSAC
       outlier detection algorithm. Is implemented and tested after the
       SciPy reference: http://wiki.scipy.org/Cookbook/RANSAC
    */
    template<typename TModelType = RansacModelLinear>
    class RANSAC
    {
public:

      /// alias for ransac() with full params
      static std::vector<std::pair<double, double> > ransac(
        const std::vector<std::pair<double, double> >& pairs, 
        const RANSACParam& p)
      {
        return ransac(pairs, p.n, p.k, p.t, p.d, p.relative_d, p.rng);
      }

      /**
        @brief This function provides a generic implementation of the RANSAC
         outlier detection algorithm. Is implemented and tested after the
         SciPy reference: http://wiki.scipy.org/Cookbook/RANSAC

        If possible, restrict 'n' to the minimal number of points which the model requires to make a
        fit, i.e. n=2 for linear, n=3 for quadratic. Any higher number will result in increasing the chance
        of including an outlier, hence a lost iteration.

        While iterating, this RANSAC implementation will consider any model which explains more data points
        than the currently best model as even better. If the data points are equal, RSS (residual sum of squared error)
        will be used.

        Making 'd' a relative measure (1-99%) is useful if you cannot predict how many points RANSAC will actually receive
        at runtime, but you have a rough idea how many percent will be outliers. E.g. if you expect 20% outliers,
        then setting d=60, relative_d=true (i.e. 60% inliers with some margin for error) is a good bet for a larger input set.
        (Consider that 2-3 data points will be used for the initial model already -- they cannot possibly become inliers).

        @param pairs Input data (paired data of type <dim1, dim2>)
        @param n The minimum number of data points required to fit the model
        @param k The maximum number of iterations allowed in the algorithm 
        @param t Threshold value for determining when a data point fits a
         model. Corresponds to the maximal squared deviation in units of the
         _second_ dimension (dim2).
        @param d The number of close data values (according to 't') required to assert that a model fits well to data
        @param relative_d Should 'd' be interpreted as percentages (0-100) of data input size
        @param rng Custom RNG function (useful for testing with fixed seeds)

        @return A vector of pairs fitting the model well; data will be unsorted
      */
  
      static std::vector<std::pair<double, double> > ransac(
          const std::vector<std::pair<double, double> >& pairs, 
          size_t n, 
          size_t k, 
          double t, 
          size_t d, 
          bool relative_d = false,
          int (*rng)(int) = nullptr)
      {
        // translate relative percentages into actual numbers
        if (relative_d)
        {
          if (d >= 100) throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("RANSAC: Relative 'd' >= 100% given. Use a lower value; the more outliers you expect, the lower it should be."));
          d = pairs.size() * d / 100;
        }

        // implementation of the RANSAC algorithm according to http://wiki.scipy.org/Cookbook/RANSAC.

        if (pairs.size() <= n)
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        String("RANSAC: Number of total data points (") + String(pairs.size()) + ") must be larger than number of initial points (n=" + String(n) + ").");
        }

        TModelType model;

        std::vector< std::pair<double, double> > alsoinliers, betterdata, bestdata;
        std::vector<std::pair<double, double> > pairs_shuffled = pairs;  // mutable data. will be shuffled in every iteration
        double besterror = std::numeric_limits<double>::max();
        typename TModelType::ModelParameters coeff;
    #ifdef DEBUG_RANSAC
        std::pair<double, double > bestcoeff;
        double betterrsq = 0;
        double bestrsq = 0;
    #endif

        for (size_t ransac_int=0; ransac_int<k; ransac_int++)
        {
          // check if the model already includes all points
          if (bestdata.size() == pairs.size()) break;

          if (rng != nullptr)
          { // use portable RNG in test mode
            std::random_shuffle(pairs_shuffled.begin(), pairs_shuffled.end(), rng);
          } else {
            std::random_shuffle(pairs_shuffled.begin(), pairs_shuffled.end());
          }

          // test 'maybeinliers'
          try
          { // fitting might throw UnableToFit if points are 'unfortunate'
            coeff = model.rm_fit(pairs_shuffled.begin(), pairs_shuffled.begin()+n);
          }
          catch (...)
          {
            continue;
          }
          // apply model to remaining data; pick inliers
          alsoinliers = model.rm_inliers(pairs_shuffled.begin()+n, pairs_shuffled.end(), coeff, t);
          // ... and add data
          if (alsoinliers.size() > d 
              || alsoinliers.size() >= (pairs_shuffled.size()-n)) // maximum number of inliers we can possibly have (i.e. remaining data)
          {
            betterdata.clear();
            std::copy( pairs_shuffled.begin(), pairs_shuffled.begin()+n, back_inserter(betterdata) );
            betterdata.insert( betterdata.end(), alsoinliers.begin(), alsoinliers.end() );
            typename TModelType::ModelParameters bettercoeff = model.rm_fit(betterdata.begin(), betterdata.end());
            double bettererror = model.rm_rss(betterdata.begin(), betterdata.end(), bettercoeff);
    #ifdef DEBUG_RANSAC
            betterrsq = model.rm_rsq(betterdata);
    #endif

            // If the current model explains more points, we assume its better (these points pass the error threshold 't', so they should be ok);
            // If the number of points is equal, we trust rss.
            // E.g. imagine gaining a zillion more points (which pass the threshold!) -- then rss will automatically be worse, no matter how good
            //      these points fit, since its a simple absolute SUM() of residual error over all points.
            if (betterdata.size() > bestdata.size() || (betterdata.size() == bestdata.size() && (bettererror < besterror)))
            {
              besterror = bettererror;
              bestdata = betterdata;
    #ifdef DEBUG_RANSAC
              bestcoeff = bettercoeff;
              bestrsq = betterrsq;
              std::cout << "RANSAC " << ransac_int << ": Points: " << betterdata.size() << " RSQ: " << bestrsq << " Error: " << besterror << " c0: " << bestcoeff.first << " c1: " << bestcoeff.second << std::endl;
    #endif
            }
          }
        }

    #ifdef DEBUG_RANSAC
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
