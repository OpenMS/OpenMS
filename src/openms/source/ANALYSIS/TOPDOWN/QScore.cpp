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

  double QScore::getQScore(const PeakGroup *pg, const int charge)
  {
    if (pg == nullptr)
    { // all zero
      return -100;
    }
    std::vector<double> weights({-0.0539, -2.6358, -8.0744, -1.7263, -0.8331, 11.8377});
    //
    //ChargeCos      -0.0539
    //ChargeSNR      -2.6358
    //Cos            -8.0744
    //SNR            -1.7263
    //ChargeScore    -0.8331
    //Intercept      11.8377

    double score = weights[weights.size() - 1];
    auto fv = toFeatureVector(pg, charge);
    for(int i=0;i<weights.size() - 1;i++){
      score += fv[i] * weights[i];
    }
    return -score;
  }

  std::vector<double> QScore::toFeatureVector(const PeakGroup *pg, const int charge)
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



  void QScore::writeAttHeader(std::fstream& f)
  {
    f<<"RT,PrecursorAvgMass,ChargeCos,ChargeSNR,Cos,SNR,ChargeScore,Qscore,Class\n";
  }

  void QScore::writeAttTsv(const double rt, const PeakGroup pg, const int charge, const bool is_identified,
                           const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, std::fstream& f)
  {
    auto fv = toFeatureVector(&pg, charge);
    if (pg.getChargeIsotopeCosine(charge) <= 0) return;

    double mass = avg.getAverageMassDelta(pg.getMonoMass()) + pg.getMonoMass();
    f << rt <<","<<mass <<",";
    for (auto& item : fv)
    {
      f<<item<<",";
    }
    f<<pg.getQScore()<<",";
    f<<(is_identified?"T":"F")<<"\n";
  }
}