// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/ML/RANSAC/RANSACModel.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ML/RANSAC/RANSACModelLinear.h>
#include <OpenMS/MATH/MathFunctions.h>

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
        : n(0), k(0), t(0), d(0), relative_d(false)
        {
        }
      /// Full constructor
      RANSACParam(size_t p_n, size_t p_k, double p_t, size_t p_d, bool p_relative_d = false)
        : n(p_n), k(p_k), t(p_t), d(p_d), relative_d(p_relative_d)
      {
        if (relative_d)
        {
          if (d >= 100) throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("RANSAC: Relative 'd' >= 100% given. Use a lower value; the more outliers you expect, the lower it should be."));
        }
      }

      [[nodiscard]] std::string toString() const
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

      explicit RANSAC(uint64_t seed = time(nullptr)):
      shuffler_(seed)
      {}

      ~RANSAC() = default;


      /// set seed for random shuffle
      void setSeed(uint64_t seed)
      {
        shuffler_.seed(seed);
      }

      /// alias for ransac() with full params
      std::vector<std::pair<double, double> > ransac(
        const std::vector<std::pair<double, double> >& pairs,
        const RANSACParam& p)
      {
        return ransac(pairs, p.n, p.k, p.t, p.d, p.relative_d);
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

        @return A vector of pairs fitting the model well; data will be unsorted
      */
      std::vector<std::pair<double, double> > ransac(
          const std::vector<std::pair<double, double> >& pairs, 
          size_t n, 
          size_t k, 
          double t, 
          size_t d, 
          bool relative_d = false)
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

          // use portable RNG in test mode
          shuffler_.portable_random_shuffle(pairs_shuffled.begin(), pairs_shuffled.end());

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

    private:
      Math::RandomShuffler shuffler_{};
    }; // class
  
  } // namespace Math


} // namespace OpenMS
