//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>

#include <iomanip>

namespace OpenMS
{
  double QScore::getQScore(const PeakGroup *pg, const int abs_charge)
  {
    if (pg->empty())
    { // all zero
      return .0;
    }
    //const std::vector<double> weights({ 1.492, -2.0041, -14.3891, -0.9853, 0.4568, 0.063, 14.4072});
    const std::vector<double> weights({ -.0941, -1.9804, -12.7522, 0.2622, -1.2431, 0.0815, 13.5244});

    //ChargeCos         1.492
    //ChargeSNR       -2.0041
    //Cos            -14.3891
    //SNR             -0.9853
    //ChargeScore      0.4568
    //AvgPPMerror       0.063
    //Intercept       14.4072

    // ChargeCos           0.0941
    // ChargeSNR           1.9804
    // Cos                12.7522
    // SNR                -0.2622
    // ChargeScore         1.2431
    // AvgPPMerror        -0.0815
    // Intercept         -13.5244

    double score = weights.back();
    auto fv = toFeatureVector_(pg, abs_charge);

    for (Size i = 0; i < weights.size() - 1; i++)
    {
      score += fv[i] * weights[i];
    }
    double qscore = 1.0 / (1.0 + exp(score));

    return qscore;
  }

  std::vector<double> QScore::toFeatureVector_(const PeakGroup *pg, const int abs_charge)
  {
    std::vector<double> fvector(6);

    double a = pg->getChargeIsotopeCosine(abs_charge);
    double d = 1;
    int index = 0;
    fvector[index++] = (log2(a + d));
    a = pg->getChargeSNR(abs_charge);
    fvector[index++] = (log2(d + a / (1 + a)));
    a = pg->getIsotopeCosine();
    fvector[index++] = (log2(a + d));
    a = pg->getSNR();
    fvector[index++] = (log2(d + a / (1 + a)));
    a = pg->getChargeScore();
    fvector[index++] = (log2(a + d));
    a = pg->getAvgPPMError();
    fvector[index++] = (log2(abs(a) + d));
    return fvector;
  }

  void QScore::writeAttCsvFromDecoyHeader(std::fstream& f)
  {
    f << "MSLevel,ChargeCos,ChargeSNR,Cos,SNR,ChargeScore,AvgPPMerror,Class\n";
  }

  void QScore::writeAttCsvFromDecoy(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f)
  {
    int ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
    String cns[] = {"T", "D", "D", "D"};
    for(auto& pg:deconvolved_spectrum)
    {
      auto fv = toFeatureVector_(&pg, pg.getRepAbsCharge());
      f<< ms_level<<",";
      for (auto& item: fv)
      {
        f << item << ",";
      }
      f <<  cns[pg.getDecoyFlag()]<< "\n";
    }
  }
}
