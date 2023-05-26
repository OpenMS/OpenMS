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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qscore.h>

namespace OpenMS
{
  PeakGroup::PeakGroup(const int min_abs_charge, const int max_abs_charge, const bool is_positive) : min_abs_charge_(min_abs_charge), max_abs_charge_(max_abs_charge), is_positive_(is_positive)
  {
  }

  bool PeakGroup::operator<(const PeakGroup& a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
    {
      return this->intensity_ < a.intensity_;
    }
    return this->monoisotopic_mass_ < a.monoisotopic_mass_;
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
    return this->monoisotopic_mass_ == a.monoisotopic_mass_ && this->intensity_ == a.intensity_;
  }

  void PeakGroup::updateAvgPPMError_()
  {
    avg_ppm_error_ = 0;
    for (auto& p : *this)
    {
      avg_ppm_error_ += getAbsPPMError_(p);
    }
    avg_ppm_error_ /= (float)size();
  }

  void PeakGroup::updateAvgDaError_()
  {
    avg_da_error_ = .0f;
    for (auto& p : *this)
    {
      avg_da_error_ += getAbsDaError_(p);
    }
    avg_da_error_ /= size();
  }

  float PeakGroup::getAbsPPMError_(const LogMzPeak& p) const
  {
    float average_mass = (float)(monoisotopic_mass_ + p.isotopeIndex * iso_da_distance_);
    return (float)(abs(average_mass / (float)p.abs_charge + FLASHDeconvHelperStructs::getChargeMass(p.is_positive) - p.mz) / p.mz * 1e6);
  }

  float PeakGroup::getAbsDaError_(LogMzPeak& p) const
  {
    float average_mass = (float)(monoisotopic_mass_ + p.isotopeIndex * iso_da_distance_);
    return (float)(abs(average_mass - p.getUnchargedMass()));
  }

  void PeakGroup::updatePerChargeCos(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg)
  {
    auto iso_dist = avg.get(monoisotopic_mass_);
    int iso_size = (int)iso_dist.size();
    auto current_per_isotope_intensities = std::vector<float>(getIsotopeIntensities().size(), .0f);

    for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
    {
      std::fill(current_per_isotope_intensities.begin(), current_per_isotope_intensities.end(), .0f);
      int min_isotope_index = (int)current_per_isotope_intensities.size();
      int max_isotope_index = -1; // this is inclusive!!

      for (auto& peak : logMzpeaks_)
      {
        if (peak.abs_charge != abs_charge)
        {
          continue;
        }

        if (peak.isotopeIndex >= (int)current_per_isotope_intensities.size())
        {
          continue;
        }

        if (peak.isotopeIndex < 0)
        {
          continue;
        }
        current_per_isotope_intensities[peak.isotopeIndex] += peak.intensity;
        min_isotope_index = min_isotope_index < peak.isotopeIndex ? min_isotope_index : peak.isotopeIndex;
        max_isotope_index = max_isotope_index < peak.isotopeIndex ? peak.isotopeIndex : max_isotope_index;
      }

      float cos_score = FLASHDeconvAlgorithm::getCosine(current_per_isotope_intensities, min_isotope_index, max_isotope_index + 1, iso_dist, iso_size, 0, 0);
      setChargeIsotopeCosine(abs_charge, cos_score); //
    }
  }

  int PeakGroup::updateIsotopeCosineSNRAvgErrorAndQscore(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double min_cos)
  {
    qscore_ = 0;
    if (empty())
    {
      return 0;
    }

    updateChargeFitScoreAndChargeIntensities_();
    if (charge_score_ < .7f) //
    {
      return 0;
    }

    updateMonomassAndIsotopeIntensities(); //
    if (per_isotope_int_.empty() || max_abs_charge_ < min_abs_charge_)
    {
      return 0;
    }

    int h_offset;
    isotope_cosine_score_ = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(monoisotopic_mass_, per_isotope_int_, h_offset, avg, -1);

    if (isotope_cosine_score_ < min_cos)
    {
      return 0;
    }

    if (h_offset != 0)
    {
      return h_offset;
    }

    updatePerChargeCos(avg);
    updateAvgPPMError_();
    updateAvgDaError_();
    updateSNR_();
    for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
    {
      if (getChargeSNR(abs_charge) <= 0 || getChargeIsotopeCosine(abs_charge) <= 0)
      {
        continue;
      }

      float q_score = Qscore::getQscore(this, abs_charge);
      if (qscore_ > q_score)
      {
        continue;
      }
      max_qscore_abs_charge_ = abs_charge;
      qscore_ = q_score;
    }
    return h_offset;
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak> PeakGroup::recruitAllPeaksInSpectrum(const MSSpectrum& spec, const double tol, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                                                                        double mono_mass, const std::unordered_set<double>& excluded_peak_mzs, int charge_offset,
                                                                                        double charge_multiple, double mz_off)
  {
    std::vector<LogMzPeak> noisy_peaks;
    if (mono_mass < 0)
    {
      return noisy_peaks;
    }

    // Adjust charge and mass information in case this is for decoy peak group generation.
    // When a decoy is generated,  charge_multiple, charge_offset, or mz_off values are different from default values.
    mono_mass = mono_mass * charge_multiple * (charge_offset + getRepAbsCharge()) / getRepAbsCharge();
    max_abs_charge_ += charge_offset;
    min_abs_charge_ += charge_offset;
    max_abs_charge_ = (int)(max_abs_charge_ * charge_multiple);
    min_abs_charge_ = (int)(min_abs_charge_ * charge_multiple);
    min_abs_charge_ = std::max(min_abs_charge_, 1);
    max_abs_charge_ = std::max(max_abs_charge_, 1);
    monoisotopic_mass_ = mono_mass;

    int iso_margin = 3; // how many isotopes do we want to scan before the monoisotopic mass?
    int max_isotope = (int)avg.getLastIndex(mono_mass);
    int min_isotope = (int)(avg.getApexIndex(mono_mass) - avg.getLeftCountFromApex(mono_mass) - iso_margin);
    min_isotope = std::max(0, min_isotope);

    clear_(); // clear logMzPeaks

    reserve((max_isotope) * (max_abs_charge_ - min_abs_charge_ + 1) * 2);
    noisy_peaks.reserve(max_isotope * (max_abs_charge_ - min_abs_charge_ + 1) * 2);

    int new_max_abs_charge = -1;              // new max abs charge after peak recruiting. The original max_abs_charge is replaced by this after recruiting.
    int new_min_abs_charge = min_abs_charge_; // new min abs charge after peak recruiting. The original min_abs_charge is replaced by this after recruiting.
    int max_sig_charge = 0;
    float max_sig = 0;

    per_charge_noise_pwr_ = std::vector<float>(1 + max_abs_charge_, .0);
    per_charge_signal_pwr_ = std::vector<float>(1 + max_abs_charge_, .0);
    per_charge_int_ = std::vector<float>(1 + max_abs_charge_, .0);

    // scan from the largest to the smallest charges and recruit the raw peaks for this monoisotopic_mass
    for (int c = max_abs_charge_; c >= min_abs_charge_; c--)
    {
      if (c <= 0)
      {
        break;
      }
      float charge_noise_pwr = .0, charge_sig_pwr = .0;
      float charge_intensity = .0;
      double cmz = (mono_mass) / c + FLASHDeconvHelperStructs::getChargeMass(is_positive_) + mz_off;
      double left_mz = (mono_mass - iso_margin * iso_da_distance_) / c + FLASHDeconvHelperStructs::getChargeMass(is_positive_) + mz_off;
      Size index = spec.findNearest(left_mz * (1 - tol));
      double iso_delta = iso_da_distance_ / c;
      float max_charge_signal_intensity = .0;

      for (; index < spec.size(); index++)
      {
        float pint = spec[index].getIntensity();
        if (pint <= 0)
        {
          continue;
        }
        double pmz = spec[index].getMZ();
        int iso_index = (int)round((pmz - cmz) / iso_delta);
        if (iso_index > max_isotope)
        {
          break;
        }

        if (iso_index < min_isotope)
        {
          continue;
        }
        // if excluded_peak_mzs_ is not empty, these mzs should be ignored in this raw spectrum for this peak group! But they can be included in noisy peaks.
        bool excluded = excluded_peak_mzs.size() > 0 && excluded_peak_mzs.find(pmz) != excluded_peak_mzs.end();

        if (!excluded && abs(pmz - cmz - iso_index * iso_delta) <= pmz * tol)
        {
          auto p = LogMzPeak(spec[index], is_positive_);
          p.isotopeIndex = iso_index;
          p.abs_charge = c;
          push_back(p);

          charge_sig_pwr += pint * pint;
          charge_intensity += pint;
          if (max_charge_signal_intensity < p.intensity)
          {
            max_charge_signal_intensity = p.intensity;
          }
        }
        else
        {
          auto p = LogMzPeak(spec[index], is_positive_);
          p.isotopeIndex = iso_index;
          p.abs_charge = c;
          noisy_peaks.push_back(p);
        }
      }
      // update charge range and per charge information (i.e., power, noise power, intensity)
      if (charge_sig_pwr > 0)
      {
        if (new_max_abs_charge < 0)
        {
          new_max_abs_charge = c;
        }
        new_min_abs_charge = c;
        for (auto& p : noisy_peaks)
        {
          if (p.abs_charge != c)
          {
            continue;
          }
          if (p.isotopeIndex < min_isotope - 1 || p.isotopeIndex > max_isotope)
          {
            continue;
          }
          // filter out all noisy peaks whose intensities are lower than 10% the maximum signal peak intensity
          if (p.intensity < max_charge_signal_intensity / 10.0)
          {
            continue;
          }
          charge_noise_pwr += p.intensity * p.intensity;
        }

        setChargePowers_(c, charge_sig_pwr, charge_noise_pwr, charge_intensity);

        if (max_sig < charge_sig_pwr / (1 + charge_noise_pwr))
        {
          max_sig = charge_sig_pwr / (1 + charge_noise_pwr);
          max_sig_charge = c;
        }
      }

      if (index >= spec.size())
      {
        break;
      }
    }

    // determine the final charge ranges based on per charge power.
    // If more than two consecutive charges do not contain any signal peak, the charge range stops at that charge.
    if (max_sig_charge > 0)
    {
      int t_nmax_abs_charge = new_max_abs_charge;
      int t_nmin_abs_charge = new_min_abs_charge;
      new_max_abs_charge = new_min_abs_charge = max_sig_charge;
      int c_zero = 0;
      float threshold = std::min(max_sig / 10, 1.0f);
      for (int z = max_sig_charge; z <= t_nmax_abs_charge; z++)
      {
        if ((per_charge_signal_pwr_[z] / (1 + per_charge_noise_pwr_[z])) < threshold)
        {
          break;
        }

        if (per_charge_signal_pwr_[z] <= 0)
        {
          c_zero++;
        }
        else
        {
          c_zero = 0;
          new_max_abs_charge = z;
        }
        if (c_zero > 2)
        {
          break;
        }
      }
      c_zero = 0;
      for (int z = max_sig_charge; z >= t_nmin_abs_charge; z--)
      {
        if ((per_charge_signal_pwr_[z] / (1 + per_charge_noise_pwr_[z])) < threshold)
        {
          break;
        }

        if (per_charge_signal_pwr_[z] <= 0)
        {
          c_zero++;
        }
        else
        {
          c_zero = 0;
          new_min_abs_charge = z;
        }
        if (c_zero > 2)
        {
          break;
        }
      }
    }
    // if the updated charge range is different from the original one, signal and noisy peaks are again updated
    if (max_abs_charge_ != new_max_abs_charge || min_abs_charge_ != new_min_abs_charge)
    {
      std::vector<LogMzPeak> new_logMzpeaks;
      new_logMzpeaks.reserve(size());
      std::vector<LogMzPeak> new_noisy_peaks;
      new_noisy_peaks.reserve(noisy_peaks.size());

      for (const auto& p : logMzpeaks_)
      {
        if (p.abs_charge < new_min_abs_charge || p.abs_charge > new_max_abs_charge)
        {
          continue;
        }
        new_logMzpeaks.push_back(p);
      }
      for (const auto& p : noisy_peaks)
      {
        if (p.abs_charge < new_min_abs_charge || p.abs_charge > new_max_abs_charge)
        {
          continue;
        }
        new_noisy_peaks.push_back(p);
      }

      new_logMzpeaks.swap(logMzpeaks_);
      new_noisy_peaks.swap(noisy_peaks);
      max_abs_charge_ = new_max_abs_charge;
      min_abs_charge_ = new_min_abs_charge;
    }

    // if no peak has been found...
    if (min_abs_charge_ > max_abs_charge_)
    {
      clear_();
    }
    else
    {
      sort();
    }
    return noisy_peaks;
  }

  void PeakGroup::updateChargeFitScoreAndChargeIntensities_()
  {
    float max_per_charge_intensity = .0;
    float summed_intensity = .0;
    int max_index = -1;
    int first_index = -1;
    int last_index = -1;

    for (int c = min_abs_charge_; c <= max_abs_charge_; c++)
    {
      summed_intensity += per_charge_int_[c];
      if (per_charge_int_[c] > 0)
      {
        if (first_index < 0)
        {
          first_index = c;
        }
        last_index = c;
      }

      if (max_per_charge_intensity > per_charge_int_[c])
      {
        continue;
      }
      max_per_charge_intensity = per_charge_int_[c];
      max_index = c;
    }
    if (max_index < 0 || summed_intensity <= 0)
    {
      charge_score_ = 0;
      return;
    }

    first_index = first_index < 0 ? 0 : first_index;
    float p = .0f;
    for (int c = max_index; c < last_index; c++)
    {
      float diff = per_charge_int_[c + 1] - per_charge_int_[c];
      if (diff <= 0)
      {
        continue;
      }
      p += diff;
    }
    for (int c = max_index; c > first_index; c--)
    {
      float diff = per_charge_int_[c - 1] - per_charge_int_[c];
      if (diff <= 0)
      {
        continue;
      }
      p += diff;
    }
    charge_score_ = std::max(.0f, 1.0f - p / summed_intensity);
  }

  void PeakGroup::updateMonomassAndIsotopeIntensities()
  {
    if (logMzpeaks_.size() == 0)
    {
      return;
    }
    int max_isotope_index = 0;
    std::sort(logMzpeaks_.begin(), logMzpeaks_.end());
    for (auto& p : logMzpeaks_)
    {
      max_isotope_index = max_isotope_index < p.isotopeIndex ? p.isotopeIndex : max_isotope_index;
    }

    per_isotope_int_ = std::vector<float>(max_isotope_index + 1, .0f);
    intensity_ = .0;
    double nominator = .0;

    for (auto& p : logMzpeaks_)
    {
      float pi = p.intensity + 1;
      if (p.isotopeIndex < 0)
      {
        continue;
      }

      per_isotope_int_[p.isotopeIndex] += p.intensity;
      nominator += pi * (p.getUnchargedMass() - p.isotopeIndex * iso_da_distance_);
      intensity_ += pi;
    }

    monoisotopic_mass_ = nominator / intensity_;
  }

  bool PeakGroup::isSignalMZ(const double mz, const double tol) const
  {
    for (auto& p : logMzpeaks_)
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

  void PeakGroup::setChargePowers_(const int abs_charge, const float signal_pwr, const float noise_pwr, const float intensity)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }

    per_charge_int_[abs_charge] = intensity;
    per_charge_noise_pwr_[abs_charge] = noise_pwr;
    per_charge_signal_pwr_[abs_charge] = signal_pwr;
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

  void PeakGroup::setMonoisotopicMass(double mono_mass)
  {
    monoisotopic_mass_ = mono_mass;
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

  void PeakGroup::Qscore(const float qscore)
  {
    qscore_ = qscore;
  }

  std::tuple<double, double> PeakGroup::getRepMzRange() const
  {
    return getMzRange(getRepAbsCharge());
  }

  std::tuple<double, double> PeakGroup::getMzRange(int abs_charge) const
  {
    double mz_start = -1;
    double mz_end = -10;
    if (!(abs_charge > max_abs_charge_ || abs_charge < min_abs_charge_))
    {
      for (auto& tmp_p : logMzpeaks_)
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
    return std::tuple<double, double> {mz_start, mz_end};
  }

  std::tuple<int, int> PeakGroup::getAbsChargeRange() const
  {
    return std::tuple<int, int> {min_abs_charge_, max_abs_charge_};
  }

  const std::vector<float>& PeakGroup::getIsotopeIntensities() const
  {
    return per_isotope_int_;
  }

  int PeakGroup::getScanNumber() const
  {
    return scan_number_;
  }

  double PeakGroup::getMonoMass() const
  {
    return monoisotopic_mass_;
  }

  float PeakGroup::getIntensity() const
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

  float PeakGroup::getQscore() const
  {
    return qscore_;
  }

  bool PeakGroup::isTargeted() const
  {
    return is_targeted_;
  }

  void PeakGroup::updateSNR_()
  {
    float cos_squred = isotope_cosine_score_ * isotope_cosine_score_;
    float signal = 0, noise = 0;
    per_charge_snr_ = std::vector<float>(1 + max_abs_charge_, .0);

    for (size_t c = min_abs_charge_; c < std::min(per_charge_signal_pwr_.size(), size_t(1 + max_abs_charge_)); ++c)
    {
      float signal_pwr = per_charge_signal_pwr_[c];
      if (per_charge_cos_.size() > c)
      {
        float per_charge_cos_squared = per_charge_cos_[c] * per_charge_cos_[c];
        float nom = per_charge_cos_squared * signal_pwr;
        float denom = 1 + per_charge_noise_pwr_[c] + (1 - per_charge_cos_squared) * signal_pwr;

        per_charge_snr_[c] = denom <= 0 ? .0f : (nom / denom);
      }
      signal += per_charge_signal_pwr_[c];
      noise += per_charge_noise_pwr_[c];
    }

    per_charge_signal_pwr_.clear();
    per_charge_noise_pwr_.clear();
    float t_nom = cos_squred * signal;
    float t_denom = 1 + noise + (1 - cos_squred) * signal;

    snr_ = t_denom <= 0 ? .0f : (t_nom / t_denom);
  }

  float PeakGroup::getQvalue(PeakGroup::TargetDummyType flag) const
  {
    if (flag == PeakGroup::TargetDummyType::target)
    {
      return std::min(1.0f, getQvalue(PeakGroup::TargetDummyType::charge_dummy) + getQvalue(PeakGroup::TargetDummyType::noise_dummy) + getQvalue(PeakGroup::TargetDummyType::isotope_dummy));
    }
    if (qvalue_.find(flag) == qvalue_.end())
    {
      return 1.0f;
    }
    return qvalue_.at(flag);
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

  float PeakGroup::getAvgDaError() const
  {
    return avg_da_error_;
  }

  float PeakGroup::getChargeSNR(const int abs_charge) const
  {
    if (abs_charge < 0 || (int)per_charge_snr_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_snr_[abs_charge];
  }

  float PeakGroup::getChargeIsotopeCosine(const int abs_charge) const
  {
    if (abs_charge < 0 || (int)per_charge_cos_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_cos_[abs_charge];
  }

  float PeakGroup::getChargeIntensity(const int abs_charge) const
  {
    if (abs_charge < 0 || (int)per_charge_int_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_int_[abs_charge];
  }

  bool PeakGroup::isPositive() const
  {
    return is_positive_;
  }

  PeakGroup::TargetDummyType PeakGroup::getTargetDummyType() const
  {
    return target_dummy_type_;
  }

  void PeakGroup::setTargetDummyType(PeakGroup::TargetDummyType index)
  {
    target_dummy_type_ = index;
  }

  void PeakGroup::setIsotopeDaDistance(const double d)
  {
    iso_da_distance_ = d;
  }

  double PeakGroup::getIsotopeDaDistance() const
  {
    return iso_da_distance_;
  }

  void PeakGroup::setIndex(const uint i)
  {
    index_ = i;
  }

  uint PeakGroup::getIndex() const
  {
    return index_;
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

  void PeakGroup::push_back(const FLASHDeconvHelperStructs::LogMzPeak& pg)
  {
    logMzpeaks_.push_back(pg);
  }

  Size PeakGroup::size() const noexcept
  {
    return logMzpeaks_.size();
  }

  void PeakGroup::clear_()
  {
    std::vector<LogMzPeak>().swap(logMzpeaks_);
  }

  void PeakGroup::reserve(Size n)
  {
    logMzpeaks_.reserve(n);
  }

  bool PeakGroup::empty() const
  {
    return logMzpeaks_.empty();
  }

  void PeakGroup::swap(std::vector<FLASHDeconvHelperStructs::LogMzPeak>& x)
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

  void PeakGroup::setQvalue(float q, PeakGroup::TargetDummyType flag)
  {
    qvalue_[flag] = q;
  }
} // namespace OpenMS
