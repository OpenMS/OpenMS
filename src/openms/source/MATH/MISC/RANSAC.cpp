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
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/RANSAC.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/CONCEPT/LogStream.h> // LOG_DEBUG

#include <numeric>
#include <boost/math/special_functions/erf.hpp>
#include <algorithm>

namespace OpenMS
{
  std::pair<double, double > Math::RANSAC::llsm_fit_(std::vector<std::pair<double, double> >& pairs)
  {
    std::vector<double> x, y;

    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      x.push_back(it->first);
      y.push_back(it->second);
    }

    // GSL implementation
    //double c0, c1, cov00, cov01, cov11, sumsq;
    //double* xa = &x[0];
    //double* ya = &y[0];
    //gsl_fit_linear(xa, 1, ya, 1, pairs.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    // OpenMS::MATH implementation
    Math::LinearRegression lin_reg;
    lin_reg.computeRegression(0.95, x.begin(), x.end(), y.begin());
    double c0, c1;
    c0 = lin_reg.getIntercept();
    c1 = lin_reg.getSlope();

    return(std::make_pair(c0,c1));
  }

  double Math::RANSAC::llsm_rsq(std::vector<std::pair<double, double> >& pairs)
  {
    std::vector<double> x, y;

    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      x.push_back(it->first);
      y.push_back(it->second);
    }

    // GSL implementation
    //double* xa = &x[0];
    //double* ya = &y[0];
    //double r = gsl_stats_correlation(xa, 1, ya, 1, pairs.size());
    //return(r);

    // OpenMS::MATH implementation
    Math::LinearRegression lin_reg;
    lin_reg.computeRegression(0.95, x.begin(), x.end(), y.begin());

    return lin_reg.getRSquared();
  }

  double Math::RANSAC::llsm_rss_(std::vector<std::pair<double, double> >& pairs, std::pair<double, double >& coefficients)
  {
    double rss = 0;

    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      rss += pow(it->second - (coefficients.first + ( coefficients.second * it->first)), 2);
    }

    return rss;
  }

  std::vector<std::pair<double, double> > Math::RANSAC::llsm_rss_inliers_(
      std::vector<std::pair<double, double> >& pairs,
      std::pair<double, double >& coefficients, double max_threshold)
  {
    std::vector<std::pair<double, double> > alsoinliers;

    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      if (pow(it->second - (coefficients.first + ( coefficients.second * it->first)), 2) < max_threshold)
      {
        alsoinliers.push_back(*it);
      }
    }

    return alsoinliers;
  }

  std::vector<std::pair<double, double> > Math::RANSAC::ransac(
      std::vector<std::pair<double, double> >& pairs, size_t n, size_t k, double t, size_t d, bool test)
  {
    // implementation of the RANSAC algorithm according to http://wiki.scipy.org/Cookbook/RANSAC.

    std::vector<std::pair<double, double> > maybeinliers, test_points, alsoinliers, betterdata, bestdata;

    double besterror = std::numeric_limits<double>::max();
    double bettererror;
#ifdef DEBUG_MRMRTNORMALIZER
    std::pair<double, double > bestcoeff;
    double betterrsq = 0;
    double bestrsq = 0;
#endif

    for (size_t ransac_int=0; ransac_int<k; ransac_int++)
    {
      std::vector<std::pair<double, double> > pairs_shuffled = pairs;

      if (!test)
      { // disables random selection in test mode
        std::random_shuffle(pairs_shuffled.begin(), pairs_shuffled.end());
      }

      maybeinliers.clear();
      test_points.clear();
      std::copy( pairs_shuffled.begin(), pairs_shuffled.begin()+n, std::back_inserter(maybeinliers) );
      std::copy( pairs_shuffled.begin()+n, pairs_shuffled.end(), std::back_inserter(test_points) );

      std::pair<double, double > coeff = Math::RANSAC::llsm_fit_(maybeinliers);

      alsoinliers = Math::RANSAC::llsm_rss_inliers_(test_points,coeff,t);

      if (alsoinliers.size() > d)
      {
        betterdata = maybeinliers;
        betterdata.insert( betterdata.end(), alsoinliers.begin(), alsoinliers.end() );
        std::pair<double, double > bettercoeff = Math::RANSAC::llsm_fit_(betterdata);
        bettererror = Math::RANSAC::llsm_rss_(betterdata,bettercoeff);
#ifdef DEBUG_MRMRTNORMALIZER
        betterrsq = Math::RANSAC::llsm_rsq(betterdata);
#endif

        if (bettererror < besterror)
        {
          besterror = bettererror;
#ifdef DEBUG_MRMRTNORMALIZER
          bestcoeff = bettercoeff;
#endif
          bestdata = betterdata;

#ifdef DEBUG_MRMRTNORMALIZER
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
  }

}
