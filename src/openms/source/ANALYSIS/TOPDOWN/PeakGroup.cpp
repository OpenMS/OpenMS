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

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

namespace OpenMS
{
  PeakGroup::PeakGroup(const int min_abs_charge, const int max_abs_charge, const bool is_positive) :
      min_abs_charge_(min_abs_charge),
      max_abs_charge_(max_abs_charge),
      is_positive_(is_positive)
  {
  }

  PeakGroup::~PeakGroup()
  {
  }

  bool PeakGroup::operator<(const PeakGroup& a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
    {
      return this->intensity_ < a.intensity_;
    }
    return this->monoisotopic_mass_ < a.monoisotopic_mass_;
    //}
    //return this->spec->getRT() < a.spec->getRT();
  }

  bool PeakGroup::operator>(const PeakGroup& a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
    {
      return this->intensity_ > a.intensity_;
    }
    return this->monoisotopic_mass_ > a.monoisotopic_mass_;
  }

  bool PeakGroup::operator==(const PeakGroup& a) const
  {
    return
        this->monoisotopic_mass_ == a.monoisotopic_mass_
&&        this->intensity_ == a.intensity_;
  }

  void PeakGroup::updateMassesAndIntensity(const int offset,
                                           const int max_isotope_index)
  {
    if (offset != 0)
    {
      std::vector<LogMzPeak> tmpPeaks;
      tmpPeaks.swap(logMzpeaks_);
      reserve(tmpPeaks.size());

      for (auto& p: tmpPeaks)
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

    for (auto& p: *this)
    {
      double pi = p.intensity + 1;
      intensity_ += pi;
      nominator += pi * (p.getUnchargedMass() - p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U);
    }
    monoisotopic_mass_ = nominator / intensity_;
    // auto massDelta = averagines.getAverageMassDelta(monoisotopicMass);
    //avgMass = monoisotopicMass + massDelta;
  }

  bool PeakGroup::isSignalMZ(const double mz, const double tol) const
  {
    for (auto& p: *this)
    {
      if (abs(p.mz - mz) < p.mz * tol * 1e-6)
      {
        return true;
      }
    }
    return false;
  }

  void PeakGroup::setSNR(const float snr)
  {
    snr_ = snr;
  }

  void PeakGroup::setChargeSNR(const int abs_charge, const float c_snr)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }
    if (per_charge_snr_.empty())
    {
      per_charge_snr_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    per_charge_snr_[abs_charge] = c_snr;
  }

  void PeakGroup::setTargeted()
  {
    is_targeted_ = true;
  }

  void PeakGroup::setChargePower(const int abs_charge, const double pwr)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }
    if (per_charge_pwr_.empty())
    {
      per_charge_pwr_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    per_charge_pwr_[abs_charge] = pwr;
  }

  void PeakGroup::setChargeSignalPower(const int abs_charge, const double pwr)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }
    if (per_charge_signal_pwr_.empty())
    {
      per_charge_signal_pwr_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    per_charge_signal_pwr_[abs_charge] = pwr;
  }

  void PeakGroup::setChargeIsotopeCosine(const int abs_charge, const float cos)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }
    if (per_charge_cos_.empty())
    {
      per_charge_cos_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    per_charge_cos_[abs_charge] = cos;
  }

  void PeakGroup::setChargeIntensity(const int abs_charge, const float intensity)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }
    if (per_charge_int_.empty())
    {
      per_charge_int_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    per_charge_int_[abs_charge] = intensity;
  }


  void PeakGroup::setMaxQScoreMzRange(const double min, const double max)
  {
    max_qscore_mz_start_ = min;
    max_qscore_mz_end_ = max;
  }

  void PeakGroup::setAbsChargeRange(const int min_abs_charge, const int max_abs_charge)
  {
    min_abs_charge_ = min_abs_charge;
    max_abs_charge_ = max_abs_charge;
  }

  void PeakGroup::setScanNumber(const int sn)
  {
    scan_number_ = sn;
  }

  void PeakGroup::setIsotopeCosine(const float cos)
  {
    isotope_cosine_score_ = cos;
  }

  void PeakGroup::setRepAbsCharge(const int max_qscore_charge)
  {
    max_qscore_abs_charge_ = max_qscore_charge;
  }

  void PeakGroup::setChargeScore(const float score)
  {
    charge_score_ = score;
  }

  void PeakGroup::setAvgPPMError(const float error)
  {
    avg_ppm_error_ = error;
  }

  void PeakGroup::setQScore(const float q)
  {
    qscore_ = q;
  }

  std::tuple<double, double> PeakGroup::getMaxQScoreMzRange() const
  {
    return std::tuple<double, double>{max_qscore_mz_start_, max_qscore_mz_end_};
  }

  std::tuple<double, double> PeakGroup::getMzRange(int abs_charge) const
  {
    double mz_start = -1;
    double mz_end = -10;
    if (!(abs_charge > max_abs_charge_ || abs_charge < min_abs_charge_))
    {
      for (auto& tmp_p: *this)
      {
        if (tmp_p.abs_charge != abs_charge)
        {
          continue;
        }
        if (mz_start < 0)
        {
          mz_start = tmp_p.mz;
        }
        else
        {
          mz_start = mz_start < tmp_p.mz ? mz_start : tmp_p.mz;
        }
        mz_end = mz_end > tmp_p.mz ? mz_end : tmp_p.mz;
      }
    }
    return std::tuple<double, double>{mz_start, mz_end};
  }

  std::tuple<int, int> PeakGroup::getAbsChargeRange() const
  {
    return std::tuple<int, int>{min_abs_charge_, max_abs_charge_};
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

  int PeakGroup::getRepAbsCharge() const
  {
    return max_qscore_abs_charge_;
  }

  float PeakGroup::getQScore() const
  {
    return qscore_;
  }

  bool PeakGroup::isTargeted() const
  {
    return is_targeted_;
  }

  void PeakGroup::updateSNR()
  {
    float cos_squred = isotope_cosine_score_ * isotope_cosine_score_;
    float signal = 0, noise = 0;
    per_charge_snr_ = std::vector<float>(1 + max_abs_charge_, .0);

    for (int c = min_abs_charge_; c <= std::min((int) per_charge_signal_pwr_.size(), max_abs_charge_); ++c)
    {
      float per_charge_cos_squared = per_charge_cos_[c] * per_charge_cos_[c];
      float signal_pwr =
          per_charge_signal_pwr_[c] < per_charge_pwr_[c] ? per_charge_signal_pwr_[c] : per_charge_pwr_[c];
      float nom = per_charge_cos_squared * signal_pwr;
      float denom = per_charge_pwr_[c] - signal_pwr
                    + (1 - per_charge_cos_squared) * signal_pwr;

      per_charge_snr_[c] = denom <= 0 ? .0 : nom / denom;

      signal += per_charge_signal_pwr_[c];
      noise += per_charge_pwr_[c] - signal_pwr;
    }

    float t_nom = cos_squred * signal;
    float t_denom = noise
                    + (1 - cos_squred) * signal;
    snr_ = t_denom <= 0 ? .0 : t_nom / t_denom;
  }

  float PeakGroup::getSNR() const
  {
    return snr_;
  }

  float PeakGroup::getChargeScore() const
  {
    return charge_score_;
  }

  float PeakGroup::getAvgPPMError() const
  {
    return avg_ppm_error_;
  }


  float PeakGroup::getChargeSNR(const int abs_charge) const
  {
    if (per_charge_snr_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_snr_[abs_charge];
  }

  float PeakGroup::getChargePower(const int abs_charge) const
  {
    if (per_charge_pwr_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_pwr_[abs_charge];
  }


  float PeakGroup::getChargeIsotopeCosine(const int abs_charge) const
  {
    if (per_charge_cos_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_cos_[abs_charge];
  }

  float PeakGroup::getChargeIntensity(const int abs_charge) const
  {
    if (per_charge_int_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_int_[abs_charge];
  }

  bool PeakGroup::isPositive() const
  {
    return is_positive_;
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator PeakGroup::begin() const noexcept
  {
    return logMzpeaks_.begin();
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator PeakGroup::end() const noexcept
  {
    return logMzpeaks_.end();
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator PeakGroup::begin() noexcept
  {
    return logMzpeaks_.begin();
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator PeakGroup::end() noexcept
  {
    return logMzpeaks_.end();
  }

  const FLASHDeconvHelperStructs::LogMzPeak& PeakGroup::operator[](const Size i) const
  {
    return logMzpeaks_[i];
  }

  void PeakGroup::push_back (const FLASHDeconvHelperStructs::LogMzPeak& pg)
  {
    logMzpeaks_.push_back(pg);
  }

  Size PeakGroup::size() const noexcept
  {
    return logMzpeaks_.size();
  }

  void PeakGroup::clear()
  {
    logMzpeaks_.clear();
  }

  void PeakGroup::reserve (Size n)
  {
    logMzpeaks_.reserve(n);
  }

  bool PeakGroup::empty() const
  {
    return logMzpeaks_.empty();
  }

  void PeakGroup::swap (std::vector<FLASHDeconvHelperStructs::LogMzPeak>& x)
  {
    logMzpeaks_.swap(x);
  }

  void PeakGroup::shrink_to_fit()
  {
    logMzpeaks_.shrink_to_fit();
  }

  void PeakGroup::sort()
  {
    std::sort(logMzpeaks_.begin(), logMzpeaks_.end());
  }
}
