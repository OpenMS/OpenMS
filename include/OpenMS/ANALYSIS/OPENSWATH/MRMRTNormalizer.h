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

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMRTNORMALIZER_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMRTNORMALIZER_H

// TODO (georger) check if needed
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

#include <numeric>
#include <boost/math/special_functions/erf.hpp>

namespace OpenMS
{

  /**
    @brief The MRMRTNormalizer will find retention time peptides in data.

    This tool will take a description of RT peptides and their normalized
    retention time to write out a transformation file on how to transform the
    RT space into the normalized space.

    The principle is adapted from the following publication:
    Escher, C. et al. (2012), Using iRT, a normalized retention time for more 
    targeted measurement of peptides. Proteomics, 12: 1111-1121.
  */
  class OPENMS_DLLAPI MRMRTNormalizer
  {

public:

    /**
      @brief This function computes a candidate outlier peptide by iteratively
       leaving one peptide out to find the one which results in the maximum R^2
       of a first order linear regression of the remaining ones. The data points
       are submitted as two vectors of doubles (x- and y-coordinates).

      @return The position of the candidate outlier peptide as supplied by the
       vector is returned.

      @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    static int outlier_candidate(std::vector<double>& x, std::vector<double>& y);

    /**
      @brief This function removes potential outliers from a set of paired points.
       Two thresholds need to be defined, first a lower R^2 limit to accept the
       regression for the RT normalization and second, the lower limit of peptide
       coverage. The algorithms then selects candidate outlier peptides and applies
       the Chauvenet's criterion on the assumption that the residuals are normal
       distributed to determine whether the peptides can be removed. This is done
       iteratively until both limits are reached.

      @param pairs Input data
      @param rsq_limit Minimal R^2 required
      @param coverage_limit Minimal coverage required (the number of points
      falls below this fraction, the algorithm aborts)
      @param use_chauvenet Whether to use Chauvenet's criterion to remove
      outliers (if false, removes any outlier candidate regardless of the
      criterion)

      @return A vector of pairs is returned if the R^2 limit was reached without
       reaching the coverage limit. If the limits are reached, an exception is
       thrown.

      @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    static std::vector<std::pair<double, double> > rm_outliers(std::vector<std::pair<double, double> >& pairs,
                                                               double rsq_limit, 
                                                               double coverage_limit, 
                                                               bool use_chauvenet);

    /**
      @brief This function computes Chauvenet's criterion probability for a vector
       and a value whose position is submitted.

      @return Chauvenet's criterion probability
    */
    static double chauvenet_probability(std::vector<double>& residuals, int pos);


    /**
      @brief This function computes Chauvenet's criterion for a vector and a value
       whose position is submitted.

      @return TRUE, if Chauvenet's criterion is fulfilled and the outlier can be removed.
    */
    static bool chauvenet(std::vector<double>& residuals, int pos);
  };

}
#endif
