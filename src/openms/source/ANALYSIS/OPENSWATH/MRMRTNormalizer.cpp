// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>

namespace OpenMS
{
  std::pair<double, double > MRMRTNormalizer::llsm_fit(std::vector<std::pair<double, double> >& pairs) {
    // standard least-squares fit to a straight line
    // takes as input a standard vector of a standard pair of points in a 2D space
    // and returns the coefficients of the linear regression Y(c,x) = c0 + c1 * x
    std::vector<double> x, y;
  
    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      x.push_back(it->first);
      y.push_back(it->second);
    }
  
    double c0, c1, cov00, cov01, cov11, sumsq;
    double* xa = &x[0];
    double* ya = &y[0];
    gsl_fit_linear(xa, 1, ya, 1, pairs.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  
    return(std::make_pair(c0,c1));
  }
  
  double MRMRTNormalizer::llsm_rsq(std::vector<std::pair<double, double> >& pairs) {
    // standard least-squares fit to a straight line
    // takes as input a standard vector of a standard pair of points in a 2D space
    // and returns the coefficients of the linear regression Y(c,x) = c0 + c1 * x
    std::vector<double> x, y;

    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      x.push_back(it->first);
      y.push_back(it->second);
    }

    double* xa = &x[0];
    double* ya = &y[0];
 
    double r = gsl_stats_correlation(xa, 1, ya, 1, pairs.size());
  
    return(pow(r,2));
  }
  
  double MRMRTNormalizer::llsm_rss(std::vector<std::pair<double, double> >& pairs, std::pair<double, double >& coefficients  ) {
    // calculates the residual sum of squares of the input points and the linear fit with coefficients c0 & c1.
    double rss = 0;
  
    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      rss += pow( it->second - (coefficients.first + ( coefficients.second * it->first)), 2);
    }
  
    return(rss);
  }
  
  std::vector<std::pair<double, double> > MRMRTNormalizer::llsm_rss_inliers(std::vector<std::pair<double, double> >&   pairs, std::pair<double, double >& coefficients, double max_threshold) {
    // calculates the residual sum of squares of the input points and the linear fit with coefficients c0 & c1.
    // further removes all points that have an error larger or equal than max_threshold.
  
    std::vector<std::pair<double, double> > alsoinliers;
  
    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      if(pow( it->second - (coefficients.first + ( coefficients.second * it->first)), 2) < max_threshold) {
        alsoinliers.push_back(*it);
      }
    }
  
    return alsoinliers;
  }
  
  std::vector<std::pair<double, double> > MRMRTNormalizer::rm_outliers_ransac(std::vector<std::pair<double, double> >&   pairs, double min_rsq, double min_coverage, int max_iterations, double max_rt_threshold) {
    // implementation of the RANSAC algorithm according to http://wiki.scipy.org/Cookbook/RANSAC.
    // adaptations for MRMRTNormalizer: estimation of parameters.
    double max_threshold = pow(max_rt_threshold,2);
    size_t coverage_num = (int)(min_coverage*pairs.size());
    size_t min_points = (size_t)(coverage_num/3);
  
    std::vector<std::pair<double, double> > maybeinliers, test_points, alsoinliers, betterdata, bestdata;
    std::pair<double, double > bestcoeff;
    double besterror = std::numeric_limits<double>::max();
    double bettererror;
    double rsq = 0;
  
    for (int ransac_int=0; ransac_int<max_iterations && rsq < min_rsq; ransac_int++) {
      std::vector<std::pair<double, double> > pairs_shuffled = pairs;
      std::random_shuffle(pairs_shuffled.begin(), pairs_shuffled.end());
  
      maybeinliers.clear();
      test_points.clear();
      std::copy( pairs_shuffled.begin(), pairs_shuffled.begin()+min_points, std::back_inserter(maybeinliers) );
      std::copy( pairs_shuffled.begin()+min_points, pairs_shuffled.end(), std::back_inserter(test_points) );
  
      std::pair<double, double > coeff = llsm_fit(maybeinliers);
  
      alsoinliers = llsm_rss_inliers(test_points,coeff,max_threshold);
  
      if (alsoinliers.size() > coverage_num) {
        betterdata = maybeinliers;
        betterdata.insert( betterdata.end(), alsoinliers.begin(), alsoinliers.end() );
        std::pair<double, double > bettercoeff = llsm_fit(betterdata);
        bettererror = llsm_rss(betterdata,bettercoeff);
  
        if (bettererror < besterror) {
          besterror = bettererror;
          bestcoeff = bettercoeff;
          bestdata = betterdata;
          rsq = llsm_rsq(betterdata);
  
          std::cout << besterror << " " << rsq << " " << ransac_int << std::endl;
        }
      }
    }
  
    return(bestdata);
  } 
 
  int MRMRTNormalizer::outlier_candidate(std::vector<double>& x, std::vector<double>& y)
  {
    // Returns candidate outlier: A linear regression and rsq is calculated for
    // the data points with one removed pair. The combination resulting in
    // highest rsq is considered corresponding to the outlier candidate. The
    // corresponding iterator position is then returned.
    std::vector<double> x_tmp, y_tmp, rsq_tmp;

    for (Size i = 0; i < x.size(); i++)
    {
      x_tmp = x;
      y_tmp = y;
      x_tmp.erase(x_tmp.begin() + i);
      y_tmp.erase(y_tmp.begin() + i);

      Math::LinearRegression lin_reg;
      lin_reg.computeRegression(0.95, x_tmp.begin(), x_tmp.end(), y_tmp.begin());

      rsq_tmp.push_back(lin_reg.getRSquared());
    }
    return max_element(rsq_tmp.begin(), rsq_tmp.end()) - rsq_tmp.begin();
  }

  std::vector<std::pair<double, double> > MRMRTNormalizer::rm_outliers(
      std::vector<std::pair<double, double> >& pairs, double rsq_limit, double coverage_limit, bool use_chauvenet)
  {
    if (pairs.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
        "Need at least 2 points for the regression.");
    }

    // Removes outliers from vector of pairs until upper rsq and lower coverage limits are reached.
    std::vector<double> x, y;
    double confidence_interval = 0.95;

    std::vector<std::pair<double, double> > pairs_corrected;

    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      x.push_back(it->first);
      y.push_back(it->second);
      LOG_DEBUG << "RT Normalization pairs: " << it->first << " : " << it->second << std::endl;
    }

    double rsq;
    rsq = 0;

    while (x.size() >= coverage_limit * pairs.size() && rsq < rsq_limit)
    {
      Math::LinearRegression lin_reg;
      lin_reg.computeRegression(confidence_interval, x.begin(), x.end(), y.begin());

      rsq = lin_reg.getRSquared();

      std::cout << "rsq: " << rsq << " points: " << x.size() << std::endl;

      if (rsq < rsq_limit)
      {
        std::vector<double> residuals;

        // calculate residuals
        for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
        {
          double intercept = lin_reg.getIntercept();
          double slope = (double)lin_reg.getSlope();
          residuals.push_back(abs(it->second - (intercept + it->first * slope)));
          LOG_DEBUG << " RT Normalization residual is " << residuals.back() << std::endl;
        }

        int pos;
        // get candidate outlier: removal of which datapoint results in best rsq?
        pos = outlier_candidate(x, y);

        // remove if residual is an outlier according to Chauvenet's criterion
        // or if testing is turned off
        LOG_DEBUG << " Got outlier candidate " << pos << "(" << x[pos] << " / " << y[pos] << std::endl;
        if (!use_chauvenet || chauvenet(residuals, pos))
        {
          x.erase(x.begin() + pos);
          y.erase(y.begin() + pos);
        }
        else
        {
          break;
        }
      }
      else
      {
        break;
      }
    }

    if (rsq < rsq_limit)
    {
      // If the rsq is below the limit, this is an indication that something went wrong!
      throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-LinearRegression-RTNormalizer", "WARNING: rsq: " + boost::lexical_cast<std::string>(rsq) + " is below limit of " + boost::lexical_cast<std::string>(rsq_limit) + ". Validate assays for RT-peptides and adjust the limit for rsq or coverage.");
    }

    for (Size i = 0; i < x.size(); i++)
    {
      pairs_corrected.push_back(std::make_pair(x[i], y[i]));
    }

    return pairs_corrected;
  }

  bool MRMRTNormalizer::chauvenet(std::vector<double>& residuals, int pos)
  {
    double criterion = 1.0 / (2 * residuals.size());
    double prob = MRMRTNormalizer::chauvenet_probability(residuals, pos);

    LOG_DEBUG << " Chauvinet testing " << prob << " < " << criterion << std::endl;
    if (prob < criterion)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  double MRMRTNormalizer::chauvenet_probability(std::vector<double>& residuals, int pos)
  {
    double mean = std::accumulate(residuals.begin(), residuals.end(), 0.0) / residuals.size();
    double stdev = std::sqrt(std::inner_product(residuals.begin(), residuals.end(), residuals.begin(), 0.0) / residuals.size() - mean * mean);

    double d = fabs(residuals[pos] - mean) / stdev;
    d /= pow(2.0, 0.5);
    double prob = boost::math::erfc(d);

    return prob;
  }

}
