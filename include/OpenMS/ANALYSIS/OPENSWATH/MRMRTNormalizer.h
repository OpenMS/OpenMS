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
// $Maintainer: Hannes Roest, George Rosenberger $
// $Authors: Hannes Roest, George Rosenberger $
// --------------------------------------------------------------------------

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
  retention time to write out a transformation file on how to transoform the
  RT space into the normalized space.

  TODO : elaborate

  */
  class MRMRTNormalizer
  {

public:

    /**
      @brief This function computes the outlier

      TODO george

      @return If an error occured during the fit.

      @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    static int outlier_candidate(std::vector<double> & x, std::vector<double> & y);

    /**
      @brief This function removes potential outliers from a set of paired points

      TODO george

      @return If an error occured during the fit.

      @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    static std::vector<std::pair<double, double> > rm_outliers(std::vector<std::pair<double, double> > & pairs, double rsq_limit, double coverage_limit);

    /**
      @brief This function computes Chauvenet's criterion probability for a vector.

      @return If an error occured during the fit.
    */
    static double chauvenet_probability(std::vector<double> & residuals, int pos);


    /**
      @brief This function computes Chauvenet's criterion for a vector.

      @return If an error occured during the fit.
    */
    static bool chauvenet(std::vector<double> & residuals, int pos);
  };

}
