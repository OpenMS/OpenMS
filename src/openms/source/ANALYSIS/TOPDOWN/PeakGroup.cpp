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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include "include/OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h"

namespace OpenMS
{
  void PeakGroup::clearChargeInfo()
  {
    std::vector<float>().swap(perChargeSNR);
    std::vector<float>().swap(perChargeCos);
    std::vector<float>().swap(perChargeInt);

    //for (auto& item : perChargeInfo)
    //{
    //  std::vector<float>().swap(item.second);
    //}
    // std::unordered_map<int, std::vector<float>>().swap(perChargeInfo);
  }

  bool PeakGroup::operator<(const PeakGroup &a) const
  {
    //if(this->spec->getRT() == a.spec->getRT()){
    if(this->monoisotopicMass == a.monoisotopicMass){
      return this->intensity < a.intensity;
    }
    return this->monoisotopicMass < a.monoisotopicMass;
    //}
    //return this->spec->getRT() < a.spec->getRT();
  }

  bool PeakGroup::operator>(const PeakGroup &a) const
  {
    if(this->monoisotopicMass == a.monoisotopicMass){
      return this->intensity > a.intensity;
    }
    return this->monoisotopicMass > a.monoisotopicMass;
  }

  bool PeakGroup::operator==(const PeakGroup &a) const
  {
    return
        this->monoisotopicMass == a.monoisotopicMass
        && this->intensity == a.intensity;
  }

  void PeakGroup::updateMassesAndIntensity(int offset,
                                           int maxIsoIndex)
  {
    //

    if (offset != 0)
    {
      std::vector<LogMzPeak> tmpPeaks;
      tmpPeaks.swap(*this);
      reserve(tmpPeaks.size());

      for (auto &p : tmpPeaks)
      {
        p.isotopeIndex -= offset;
        if (p.isotopeIndex < 0 || p.isotopeIndex >= maxIsoIndex)
        {
          continue;
        }
        push_back(p);
      }
    }

    intensity = .0;
    double nominator = .0;

    for (auto &p : *this)
    {
      double pi = p.intensity;
      intensity += pi;
      nominator += pi * (p.getUnchargedMass() - p.isotopeIndex * Constants::C13C12_MASSDIFF_U);
    }
    monoisotopicMass = nominator / intensity;
    // auto massDelta = averagines.getAverageMassDelta(monoisotopicMass);
    //avgMass = monoisotopicMass + massDelta;

  }


  float PeakGroup::getChargeSNR(int charge) const
  {
    if (maxCharge < charge)
    {
      return 0;
    }
    return perChargeSNR[charge];
  }

  float PeakGroup::getChargeIsotopeCosine(int charge) const
  {
    if (maxCharge < charge)
    {
      return 0;
    }
    return perChargeCos[charge];
  }

  float PeakGroup::getChargeIntensity(int charge) const
  {
    if (maxCharge < charge)
    {
      return 0;
    }
    return perChargeInt[charge];
  }

  void PeakGroup::setChargeSNR(int charge, float snr)
  {
    if (maxCharge < charge)
    {
      return;
    }
    perChargeSNR[charge] = snr;
  }

  void PeakGroup::setChargeIsotopeCosine(int charge, float cos)
  {
    if (maxCharge < charge)
    {
      return;
    }
    perChargeCos[charge] = cos;
  }

  void PeakGroup::setChargeIntensity(int charge, float i)
  {
    if (maxCharge < charge)
    {
      return;
    }
    perChargeInt[charge] = i;
  }

  PeakGroup::PeakGroup(int maxCharge) :
      maxCharge(maxCharge)
  {
    perChargeSNR = std::vector<float>(1 + maxCharge, .0);
    perChargeInt = std::vector<float>(1 + maxCharge, .0);
    perChargeCos = std::vector<float>(1 + maxCharge, .0);
  }

  void PeakGroup::setMaxQScoreMzRange(double min, double max)
  {
    maxQScoreMzStart = min;
    maxQScoreMzEnd = max;
  }

  std::tuple<double, double> PeakGroup::getMzxQScoreMzRange() const
  {
    return std::tuple<double, double>{maxQScoreMzStart, maxQScoreMzEnd};
  }

  std::tuple<int, int> PeakGroup::getChargeRange() const
  {
    return std::tuple<int, int>{minCharge, maxCharge};
  }

  void PeakGroup::setChargeRange(int min, int max)
  {
    minCharge = min;
    maxCharge = max;
  }

  void PeakGroup::setScanNumber(int sn)
  {
    scanNumber = sn;
  }

  int PeakGroup::getScanNumber() const
  {
    return scanNumber;
  }

  void PeakGroup::setMonoMass(double m)
  {
    monoisotopicMass = m;
  }

  double PeakGroup::getMonoMass() const
  {
    return monoisotopicMass;
  }

  void PeakGroup::setIntensity(double i)
  {
    intensity = i;
  }

  double PeakGroup::getIntensity() const
  {
    return intensity;
  }

  void PeakGroup::setIsotopeCosine(float cos)
  {
    isotopeCosineScore = cos;
  }

  void PeakGroup::setRepCharge(int c)
  {
    maxQScoreCharge = c;
  }

  float PeakGroup::getIsotopeCosine() const
  {
    return isotopeCosineScore;
  }

  int PeakGroup::getRepCharge() const
  {
    return maxQScoreCharge;
  }

  void PeakGroup::setQScore(float q)
  {
    qScore = q;
  }

  float PeakGroup::getQScore() const
  {
    return qScore;
  }

  void PeakGroup::setSNR(float snr)
  {
    totalSNR = snr;
  }

  float PeakGroup::getSNR() const
  {
    return totalSNR;
  }


}