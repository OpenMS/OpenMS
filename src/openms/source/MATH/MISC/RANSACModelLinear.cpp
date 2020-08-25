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

#include <OpenMS/MATH/MISC/RANSACModelLinear.h>

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <numeric>


namespace OpenMS
{
  namespace Math
  {
    RansacModelLinear::ModelParameters RansacModelLinear::rm_fit_impl(const DVecIt& begin, const DVecIt& end)
    {
      std::vector<double> x, y;

      for (DVecIt it = begin; it != end; ++it)
      {
        x.push_back(it->first);
        y.push_back(it->second);
      }
      LinearRegression lin_reg;
      lin_reg.computeRegression(0.95, x.begin(), x.end(), y.begin(), false); // no goodness of fit computation
      ModelParameters p;
      p.push_back(lin_reg.getIntercept());
      p.push_back(lin_reg.getSlope());
      return p;
    }

    double RansacModelLinear::rm_rsq_impl(const DVecIt& begin, const DVecIt& end)
    {
      std::vector<double> x, y;

      for (DVecIt it = begin; it != end; ++it)
      {
        x.push_back(it->first);
        y.push_back(it->second);
      }

      LinearRegression lin_reg;
      lin_reg.computeRegression(0.95, x.begin(), x.end(), y.begin(), false);

      return lin_reg.getRSquared();
    }

    double RansacModelLinear::rm_rss_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients)
    {
      double rss = 0;

      for (DVecIt it = begin; it != end; ++it)
      {
        rss += pow(it->second - (coefficients[0] + ( coefficients[1] * it->first)), 2);
      }

      return rss;
    }

    RansacModelLinear::DVec RansacModelLinear::rm_inliers_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients, double max_threshold)
    {
      DVec alsoinliers;
      //std::cerr << "\n\nRANSAC dists: ";
      for (DVecIt it = begin; it != end; ++it)
      {
        double dist = pow(it->second - (coefficients[0] + ( coefficients[1] * it->first)), 2);
        //std::cerr << dist << ", ";
        if (dist < max_threshold)
        {
          alsoinliers.push_back(*it);
        }
      }
      //std::cerr << "\n\n";

      return alsoinliers;
    }


  } // Math


} // OpenMS
