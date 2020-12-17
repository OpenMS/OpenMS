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

  double QScore::getQScore(PeakGroup *pg, int charge)
  {
    if (pg == nullptr)
    { // all zero
      return -100;
    }
    std::vector<double> weights({-0.7927, -1.9353, -9.4241, -1.5728, -2.4232, 14.0816 });
    //
    //ChargeCos      -0.7927
    //ChargeSNR      -1.9353
    //Cos            -9.4241
    //SNR            -1.5728
    //ChargeScore    -2.4232
    //Intercept      14.0816

    //ChargeCos      -5.6679
    //ChargeSNR      -0.6528
    //Cos             2.5078
    //SNR            -1.5879
    //ChargeScore    -3.6178
    //Intercept       7.5809

    double score = weights[weights.size() - 1];
    auto fv = toFeatureVector(pg, charge);
    for(int i=0;i<weights.size() - 1;i++){
      score += fv[i] * weights[i];
    }
    return -score;
  }

  std::vector<double> QScore::toFeatureVector(PeakGroup *pg, int charge)
  {
    std::vector<double> fvector;

    double a = pg->getChargeIsotopeCosine(charge);
    double d = 1;
    fvector.push_back(log2(a + d));
    a = pg->getChargeSNR(charge);
    fvector.push_back(log2(d + a/(1+a)));
    a = pg->getIsotopeCosine();
    fvector.push_back(log2(a + d));
    a = pg->getSNR();
    fvector.push_back(log2(d + a/(1+a)));
    a = pg->getChargeScore();
    fvector.push_back(log2(a + d));

    return fvector;
  }



  void QScore::writeAttHeader(std::fstream &f)
  {
    f<<"RT,PrecursorAvgMass,ChargeCos,ChargeSNR,Cos,SNR,ChargeScore,Qscore,Class\n";
  }

  void QScore::writeAttTsv(double rt, PeakGroup pg, int charge, bool isIdentified, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, std::fstream &f)
  {
    auto fv = toFeatureVector(&pg, charge);
    if (pg.getChargeIsotopeCosine(charge) <= 0) return;

    auto mass = avg.getAverageMassDelta(pg.getMonoMass()) + pg.getMonoMass();
    f << rt <<","<<mass <<",";
    for (auto& item : fv)
    {
      f<<item<<",";
    }
    f<<pg.getQScore()<<",";
    f<<(isIdentified?"T":"F")<<"\n";
  }
}