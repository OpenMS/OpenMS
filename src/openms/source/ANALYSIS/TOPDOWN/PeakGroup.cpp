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

  PeakGroup::~PeakGroup()
  {
    std::vector<LogMzPeak>().swap(peaks);
    clearChargeInfo();
  }

  void PeakGroup::push_back(FLASHDeconvHelperStructs::LogMzPeak &p)
  {
    peaks.push_back(p);
  }

  void PeakGroup::reserve(Size n)
  {
    peaks.reserve(n);
  }

  void PeakGroup::clearChargeInfo()
  {
    for (auto& item : perChargeInfo)
    {
      std::vector<float>().swap(item.second);
    }
    std::unordered_map<int, std::vector<float>>().swap(perChargeInfo);
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

  void PeakGroup::updateMassesAndIntensity(FLASHDeconvHelperStructs::PrecalculatedAveragine &averagines,
                                           double chargeMass,
                                           int offset,
                                           int maxIsoIndex)
  {
    //

    if (offset != 0)
    {
      std::vector<LogMzPeak> tmpPeaks;
      tmpPeaks.swap(peaks);
      peaks.reserve(tmpPeaks.size());

      for (auto &p : tmpPeaks)
      {
        p.isotopeIndex -= offset;
        if (p.isotopeIndex < 0 || p.isotopeIndex >= maxIsoIndex)
        {
          continue;
        }
        peaks.push_back(p);
      }
    }

    intensity = .0;
    double nominator = .0;

    for (auto &p : peaks)
    {
      double pi = p.intensity;
      intensity += pi;
      nominator += pi * (p.getUnchargedMass(chargeMass) - p.isotopeIndex * Constants::C13C12_MASSDIFF_U);

    }
    monoisotopicMass = nominator / intensity;
    auto massDelta = averagines.getAverageMassDelta(monoisotopicMass);
    avgMass = monoisotopicMass + massDelta;

  }
}