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
  PeakGroup::PeakGroup(const int min_charge, const int max_charge) :
      min_charge_(min_charge),
      max_charge_(max_charge)
  {
    per_charge_snr_ = std::vector<float>(1 + max_charge, .0);
    per_charge_int_ = std::vector<float>(1 + max_charge, .0);
    per_charge_cos_ = std::vector<float>(1 + max_charge, .0);
  }

  PeakGroup::~PeakGroup()
  {
    clearChargeInfo();
  }

  void PeakGroup::clearChargeInfo()
  {
    std::vector<float>().swap(per_charge_snr_);
    std::vector<float>().swap(per_charge_cos_);
    std::vector<float>().swap(per_charge_int_);
  }

  bool PeakGroup::operator<(const PeakGroup& a) const
  {
    //if(this->spec->getRT() == a.spec->getRT()){
    if(this->monoisotopic_mass_ == a.monoisotopic_mass_){
      return this->intensity_ < a.intensity_;
    }
    return this->monoisotopic_mass_ < a.monoisotopic_mass_;
    //}
    //return this->spec->getRT() < a.spec->getRT();
  }

  bool PeakGroup::operator>(const PeakGroup& a) const
  {
    if(this->monoisotopic_mass_ == a.monoisotopic_mass_){
      return this->intensity_ > a.intensity_;
    }
    return this->monoisotopic_mass_ > a.monoisotopic_mass_;
  }

  bool PeakGroup::operator==(const PeakGroup& a) const
  {
    return
        this->monoisotopic_mass_ == a.monoisotopic_mass_
       && this->intensity_ == a.intensity_;
  }

  void PeakGroup::updateMassesAndIntensity(const int offset,
                                           const int max_isotope_index)
  {
    if (offset != 0)
    {
      std::vector<LogMzPeak> tmpPeaks;
      tmpPeaks.swap(*this);
      reserve(tmpPeaks.size());

      for (auto& p : tmpPeaks)
      {
        p.isotopeIndex -= offset;
        if (p.isotopeIndex < 0 || p.isotopeIndex >= max_isotope_index)
        {
          continue;
        }
        push_back(p);
      }
    }

    intensity_ = .0;
    double nominator = .0;

    for (auto& p : *this)
    {
      double pi = p.intensity;
      intensity_ += pi;
      nominator += pi * (p.getUnchargedMass() - p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U);
    }
    monoisotopic_mass_ = nominator / intensity_;
    // auto massDelta = averagines.getAverageMassDelta(monoisotopicMass);
    //avgMass = monoisotopicMass + massDelta;

  }

  void PeakGroup::setChargeSNR(const int charge, const float snr)
  {
    if (max_charge_ < charge)
    {
      return;
    }
    per_charge_snr_[charge] = snr;
  }

  void PeakGroup::setChargeIsotopeCosine(const int charge, const float cos)
  {
    if (max_charge_ < charge)
    {
      return;
    }
    per_charge_cos_[charge] = cos;
  }

  void PeakGroup::setChargeIntensity(const int charge, const float intensity)
  {
    if (max_charge_ < charge)
    {
      return;
    }
    per_charge_int_[charge] = intensity;
  }


  void PeakGroup::setMaxQScoreMzRange(const double min, const double max)
  {
    max_qscore_mz_start_ = min;
    max_qscore_mz_end_ = max;
  }

  void PeakGroup::setChargeRange(const int min, const int max)
  {
    min_charge_ = min;
    max_charge_ = max;
  }

  void PeakGroup::setScanNumber(const int sn)
  {
    scan_number_ = sn;
  }

  void PeakGroup::setIsotopeCosine(const float cos)
  {
    isotope_cosine_score_ = cos;
  }

  void PeakGroup::setRepCharge(const int c)
  {
    max_qscore_charge_ = c;
  }

  void PeakGroup::setChargeScore(const float score)
  {
    charge_score_ = score;
  }


  void PeakGroup::setSNR(const float snr)
  {
    total_snr_ = snr;
  }

  void PeakGroup::setQScore(const float q)
  {
    qscore_ = q;
  }

  std::tuple<double, double> PeakGroup::getMaxQScoreMzRange() const
  {
    return std::tuple<double, double>{max_qscore_mz_start_, max_qscore_mz_end_};
  }

  std::tuple<int, int> PeakGroup::getChargeRange() const
  {
    return std::tuple<int, int>{min_charge_, max_charge_};
  }

  int PeakGroup::getScanNumber() const
  {
    return scan_number_;
  }

  double PeakGroup::getMonoMass() const
  {
    return monoisotopic_mass_;
  }

  double PeakGroup::getIntensity() const
  {
    return intensity_;
  }
  float PeakGroup::getIsotopeCosine() const
  {
    return isotope_cosine_score_;
  }

  int PeakGroup::getRepCharge() const
  {
    return max_qscore_charge_;
  }

  float PeakGroup::getQScore() const
  {
    return qscore_;
  }


  float PeakGroup::getSNR() const
  {
    return total_snr_;
  }

  float PeakGroup::getChargeScore() const
  {
    return charge_score_;
  }


  float PeakGroup::getChargeSNR(const int charge) const
  {
    if (max_charge_ < charge)
    {
      return 0;
    }
    return per_charge_snr_[charge];
  }

  float PeakGroup::getChargeIsotopeCosine(const int charge) const
  {
    if (max_charge_ < charge)
    {
      return 0;
    }
    return per_charge_cos_[charge];
  }

  float PeakGroup::getChargeIntensity(const int charge) const
  {
    if (max_charge_ < charge)
    {
      return 0;
    }
    return per_charge_int_[charge];
  }
}
