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
// $Maintainer: Witold Wolski $
// $Authors: Witold Wolski  $
// --------------------------------------------------------------------------

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/StatsHelpers.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace OpenSwath
{
  void normalize(
    const std::vector<double> & intensities,
    double normalizer,
    std::vector<double> & normalized_intensities)
  {
    //compute total intensities
    //normalize intensities
    normalized_intensities.resize(intensities.size());
    if (normalizer > 0)
    {
      std::transform(intensities.begin(), intensities.end(), normalized_intensities.begin(), std::bind2nd(std::divides<double>(), normalizer));
    }
  }


  double dotprodScoring(std::vector<double> intExp, std::vector<double> theorint)
  {
    for (unsigned int i = 0; i < intExp.size(); ++i)
    {
      intExp[i] = sqrt(intExp[i]);
      theorint[i] = sqrt(theorint[i]);
    }

    double intExptotal = norm(intExp.begin(), intExp.end());
    double intTheorTotal = norm(theorint.begin(), theorint.end());
    OpenSwath::normalize(intExp, intExptotal, intExp);
    OpenSwath::normalize(theorint, intTheorTotal, theorint);
    double score2 = OpenSwath::dotProd(intExp.begin(), intExp.end(), theorint.begin());
    return score2;
  }

  double manhattanScoring(std::vector<double> intExp, std::vector<double> theorint)
  {

    for (unsigned int i = 0; i < intExp.size(); ++i)
    {
      intExp[i] = sqrt(intExp[i]);
      theorint[i] = sqrt(theorint[i]);
      //std::transform(intExp.begin(), intExp.end(), intExp.begin(), sqrt);
      //std::transform(theorint.begin(), theorint.end(), theorint.begin(), sqrt);
    }

    double intExptotal = std::accumulate(intExp.begin(), intExp.end(), 0.0);
    double intTheorTotal = std::accumulate(theorint.begin(), theorint.end(), 0.0);
    OpenSwath::normalize(intExp, intExptotal, intExp);
    OpenSwath::normalize(theorint, intTheorTotal, theorint);
    double score2 = OpenSwath::manhattanDist(intExp.begin(), intExp.end(), theorint.begin());
    return score2;
  }

}
