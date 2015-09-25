// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#ifndef OPENMS_MATH_MISC_RANSAC_H
#define OPENMS_MATH_MISC_RANSAC_H

#include <OpenMS/config.h>

#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{

  namespace Math
  {

    /**
      @brief This class provides a generic implementation of the RANSAC
       outlier detection algorithm. Is implemented and tested after the
       SciPy reference: http://wiki.scipy.org/Cookbook/RANSAC
    */
    class OPENMS_DLLAPI RANSAC
    {

protected:

    static std::pair<double, double > llsm_fit_(std::vector<std::pair<double, double> >& pairs);
  
    /// interface for GSL or OpenMS::MATH linear regression implementation
    /// calculates the residual sum of squares of the input points and the linear fit with coefficients c0 & c1.
    static double llsm_rss_(std::vector<std::pair<double, double> >& pairs, std::pair<double, double >& coefficients);
  
    /// calculates the residual sum of squares of the input points and the linear fit with coefficients c0 & c1.
    /// further removes all points that have an error larger or equal than max_threshold.
    static std::vector<std::pair<double, double> > llsm_rss_inliers_(std::vector<std::pair<double, double> >& pairs,
        std::pair<double, double >& coefficients, double max_threshold);

public:

    /**
      @brief This function provides a generic implementation of the RANSAC
       outlier detection algorithm. Is implemented and tested after the
       SciPy reference: http://wiki.scipy.org/Cookbook/RANSAC

      @param pairs Input data (paired data of type <dim1, dim2>)
      @param n the minimum number of data points required to fit the model
      @param k the maximum number of iterations allowed in the algorithm 
      @param t a threshold value for determining when a data point fits a
       model. Corresponds to the maximal squared deviation in units of the
       _second_ dimension (dim2).
      @param d the number of close data values required to assert that a model fits well to data
      @param test disables the random component of the algorithm

      @return A vector of pairs
    */
    static std::vector<std::pair<double, double> > ransac(std::vector<std::pair<double, double> >& pairs, size_t n, size_t k, double t, size_t d, bool test = false);

    /**
      @brief Interface for GSL or OpenMS::MATH linear regression implementation
      standard least-squares fit to a straight line takes as input a standard
      vector of a standard pair of points in a 2D space and returns the
      R-squared value for the corresponding linear regression Y(c,x) = c0 + c1 * x
    */
    static double llsm_rsq(std::vector<std::pair<double, double> >& pairs);

    };

  }


}
#endif // OPENMS_MATH_MISC_RANSAC_H
