// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

  int MRMRTNormalizer::outlier_candidate(std::vector<double> & x, std::vector<double> & y)
  {
    // Returns candidate outlier: A linear regression and rsq is calculated for the data points with one removed pair. The combination resulting in highest rsq is considered corresponding to the outlier candidate. The corresponding iterator position is then returned.
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

  std::vector<std::pair<double, double> > MRMRTNormalizer::rm_outliers(std::vector<std::pair<double, double> > & pairs, double rsq_limit, double coverage_limit)
  {
    // Removes outliers from vector of pairs until upper rsq and lower coverage limits are reached.
    std::vector<double> x, y;
    double confidence_interval = 0.95;

    std::vector<std::pair<double, double> > pairs_corrected;

    for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
    {
      x.push_back(it->first);
      y.push_back(it->second);
    }

    double rsq;
    rsq = 0;

    while (x.size() >= coverage_limit * pairs.size() && rsq < rsq_limit)
    {
      Math::LinearRegression lin_reg;
      lin_reg.computeRegression(confidence_interval, x.begin(), x.end(), y.begin());

      rsq = lin_reg.getRSquared();

      cout << "rsq: " << rsq << " points: " << x.size() << endl;

      if (rsq < rsq_limit)
      {
        double intercept;
        double slope;
        std::vector<double> residuals;

        // calculate residuals
        for (std::vector<std::pair<double, double> >::iterator it = pairs.begin(); it != pairs.end(); ++it)
        {
          intercept = lin_reg.getIntercept();
          slope = (double)lin_reg.getSlope();
          residuals.push_back(abs(it->second - (intercept + it->first * slope)));
        }

        int pos;
        // get candidate outlier: removal of which datapoint results in best rsq?
        pos = outlier_candidate(x, y);

        if (chauvenet(residuals, pos)) // test for outlier: is the residual an outlier according to Chauvenet's criterion?
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
      throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-LinearRegression-RTNormalizer","WARNING: rsq: " + boost::lexical_cast<std::string>(rsq) + " is below limit of " + boost::lexical_cast<std::string>(rsq_limit) + ". Validate assays for RT-peptides and adjust the limit for rsq or coverage.");
    }

    for (Size i = 0; i < x.size(); i++)
    {
      pairs_corrected.push_back(make_pair(x[i], y[i]));
    }

    return pairs_corrected;
  }

  bool MRMRTNormalizer::chauvenet(std::vector<double> & residuals, int pos)
  {
    double criterion = 1.0 / (2 * residuals.size());
    double prob = MRMRTNormalizer::chauvenet_probability(residuals,pos);

    if (prob < criterion)
    {
      return(true);
    }
    else
    {
      return(false);
    }
  }

  double MRMRTNormalizer::chauvenet_probability(std::vector<double> & residuals, int pos)
  {
    double mean = std::accumulate(residuals.begin(), residuals.end(), 0.0) / residuals.size();
    double stdev = std::sqrt(std::inner_product(residuals.begin(), residuals.end(), residuals.begin(), 0.0) / residuals.size() - mean * mean);

    double d = abs(residuals[pos]-mean)/stdev;
    d /= pow(2.0,0.5);
    double prob=boost::math::erfc(d);

    return(prob);
  }
}
