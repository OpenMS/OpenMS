//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include "OpenMS/ANALYSIS/TOPDOWN/QScore.h"
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

namespace OpenMS
{

  double QScore::getQScore(const PeakGroup *pg, const int abs_charge)
  {
    if (pg == nullptr)
    { // all zero
      return -100;
    }
    const std::vector<double> weights_vh({1.3522, -1.0877, -16.4956, -2.036, -0.9439, 18.251});
    const std::vector<double> weights_h({-0.1738, -1.4679, -6.8065, -1.2966, -0.9714, 9.2687});
    const std::vector<double> weights_l({-3.203, -2.6899, 11.1909, -3.1146, -1.9595, -2.3368});

    //
    //ChargeCos      -0.1738 <<30kda
    //ChargeSNR      -1.4679
    //Cos            -6.8065
    //SNR            -1.2966
    //ChargeScore    -0.9714
    //Intercept       9.2687

    //ChargeCos        1.3522 ///>30kda
    //ChargeSNR       -1.0877
    //Cos            -16.4956
    //SNR              -2.036
    //ChargeScore     -0.9439
    //Intercept        18.251

    //ChargeCos          -3.203 // < abs_charge 6
    //ChargeSNR         -2.6899
    //Cos               11.1909
    //SNR               -3.1146
    //ChargeScore       -1.9595
    //Intercept         -2.3368
    const std::vector<double> &weights = (abs_charge > 6 ?
                                          (pg->getMonoMass() > 30000.0 ? weights_vh : weights_h) :
                                          weights_l);
    double score = weights[weights.size() - 1];
    auto fv = toFeatureVector_(pg, abs_charge);
    for (int i = 0; i < weights.size() - 1; i++)
    {
      score += fv[i] * weights[i];
    }
    return -score;
  }

  std::vector<double> QScore::toFeatureVector_(const PeakGroup *pg, const int abs_charge)
  {
    std::vector<double> fvector;

    double a = pg->getChargeIsotopeCosine(abs_charge);
    double d = 1;
    fvector.push_back(log2(a + d));
    a = pg->getChargeSNR(abs_charge);
    fvector.push_back(log2(d + a / (1 + a)));
    a = pg->getIsotopeCosine();
    fvector.push_back(log2(a + d));
    a = pg->getSNR();
    fvector.push_back(log2(d + a/(1+a)));
    a = pg->getChargeScore();
    fvector.push_back(log2(a + d));

    return fvector;
  }



  void QScore::writeAttHeader(std::fstream& f)
  {
    f
        << "ACC,RT,PrecursorMonoMass,PrecursorAvgMass,PrecursorMz,PrecursorCharge,ChargeCos,ChargeSNR,Cos,SNR,ChargeScore,Qscore,Class\n";
  }

  void QScore::writeAttTsv(const String &acc,
                           const double rt,
                           const double pmass,
                           const double pmz,
                           const PeakGroup pg,
                           const int charge,
                           const bool is_identified,
                           const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                           std::fstream &f)
  {
    auto avgpmass = avg.getAverageMassDelta(pmass) + pmass;
    if (pg.empty())
    {
      // return;
      f << acc << "," << rt << "," << (pmass <= .0 ? 0 : pmass) << "," << (pmass <= .0 ? 0 : avgpmass) << "," << pmz
        << ",0,";
      f << "0,0,0,0,0,-5,";
      f << (is_identified ? "T" : "F") << "\n";
    }
    else
    {
      auto fv = toFeatureVector_(&pg, charge);
      //if (pg.getChargeIsotopeCosine(charge) <= 0)
      //  return;

      double monomass = pmass <= .0? pg.getMonoMass() : pmass;
      double mass = pmass <= .0? avg.getAverageMassDelta(pg.getMonoMass()) + pg.getMonoMass() : avgpmass;
      f << acc << "," << rt << "," << monomass << "," << mass << "," << pmz << "," << pg.getRepAbsCharge() << ",";
      for (auto &item : fv)
      {
        f << item << ",";
      }
      f << pg.getQScore() << ",";
      f << (is_identified ? "T" : "F") << "\n";
    }
  }
}
