// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qscore.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

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
    for (const auto& p : *this)
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

  float PeakGroup::getAbsDaError_(const LogMzPeak& p) const
  {
    float average_mass = (float)(monoisotopic_mass_ + p.isotopeIndex * iso_da_distance_);
    return (float)(abs(average_mass - p.getUnchargedMass()));
  }

  int PeakGroup::getMinNegativeIsotopeIndex() const
  {
    return min_negative_isotope_index_;
  }

  void PeakGroup::updatePerChargeCos_(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg)
  {
    auto iso_dist = avg.get(monoisotopic_mass_);
    auto current_per_isotope_intensities = std::vector<float>(getIsotopeIntensities().size() + min_negative_isotope_index_, .0f);

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

      float cos_score =
        SpectralDeconvolution::getCosine(current_per_isotope_intensities, min_isotope_index, max_isotope_index + 1, iso_dist, 0, 0, target_decoy_type_ == PeakGroup::TargetDecoyType::noise_decoy);
      setChargeIsotopeCosine(abs_charge, cos_score); //
    }
  }

  int PeakGroup::updateQscore(std::vector<LogMzPeak>& noisy_peaks, const MSSpectrum& spec, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double min_cos, bool is_low_charge, int allowed_iso_error)
  {
    qscore_ = 0;

    if (empty())
    {
      return 0;
    }

    updatePerChargeInformation_(noisy_peaks);
    updateChargeRange_(noisy_peaks);
    updateChargeFitScoreAndChargeIntensities_(is_low_charge);

    if (charge_score_ < .6f) //
    {
      return 0;
    }

    updateMonoMassAndIsotopeIntensities(); //
    if (per_isotope_int_.empty() || max_abs_charge_ < min_abs_charge_)
    {
      return 0;
    }

    int h_offset;
    isotope_cosine_score_ =
      SpectralDeconvolution::getIsotopeCosineAndDetermineIsotopeIndex(monoisotopic_mass_, per_isotope_int_, h_offset, avg, -min_negative_isotope_index_, // change if to select cosine calculation and if to get second best hits
                                                                      target_decoy_type_ == PeakGroup::TargetDecoyType::isotope_decoy? 0 : -1, allowed_iso_error, target_decoy_type_ == PeakGroup::TargetDecoyType::isotope_decoy ? PeakGroup::TargetDecoyType::target : target_decoy_type_);

    if (h_offset != 0) return h_offset;

    if (isotope_cosine_score_ < min_cos)
    {
      return h_offset;
    }

    if (max_abs_charge_ - min_abs_charge_ < max_abs_charge_ / 20) // if charge range is too small ...
    {
      return h_offset;
    }

    updatePerChargeCos_(avg);
    updateAvgPPMError_();
    updateAvgDaError_();
    updateSNR_();

    for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
    {
      if (getChargeSNR(abs_charge) <= 0 || getChargeIsotopeCosine(abs_charge) <= 0)
      {
        continue;
      }

      double q_score = Qscore::getQscore(this, spec);
      if (qscore_ < q_score)
      {
        qscore_ = q_score;
      }

      if (getChargeSNR(abs_charge) > getChargeSNR(max_snr_abs_charge_))
      {
        max_snr_abs_charge_ = abs_charge;
      }
    }
    return h_offset;
  }

  float PeakGroup::getNoisePeakPower_(const std::vector<FLASHDeconvHelperStructs::LogMzPeak>& noisy_peaks, const std::vector<FLASHDeconvHelperStructs::LogMzPeak>& signal_peaks, const int z) const
  {
    if (noisy_peaks.empty())
      return 0;
    const Size max_noisy_peaks = 50; // too many noise peaks will slow down the process
    const Size max_bin_number = 29;  // 24 bin + 5 extra bin
    float threshold = 0;

    std::vector<FLASHDeconvHelperStructs::LogMzPeak> all_peaks;
    std::set<double> signal_mzs;

    all_peaks.reserve(max_noisy_peaks + signal_peaks.size());

    int noise_peak_count = 0;
    for (const auto& noisy_peak : noisy_peaks)
    {
      if (noisy_peak.abs_charge != z)
        continue;
      noise_peak_count++;
    }
    // get intensity threshold
    if (noise_peak_count > max_noisy_peaks)
    {
      std::vector<float> intensities;
      intensities.reserve(noise_peak_count);
      for (const auto& noisy_peak : noisy_peaks)
      {
        if (noisy_peak.abs_charge != z) continue;
        intensities.push_back(noisy_peak.intensity);
      }

      std::sort(intensities.begin(), intensities.end());
      threshold = intensities[intensities.size() - max_noisy_peaks];
    }

    // filter peaks and check which mzs are signal and which are noise.
    for (const auto& p : noisy_peaks)
    {
      if (p.abs_charge != z || p.intensity < threshold)
      {
        continue;
      }
      all_peaks.push_back(p);
    }

    for (const auto& p : signal_peaks)
    {
      if (p.abs_charge != z || p.intensity < threshold)
      {
        continue;
      }
      all_peaks.push_back(p);
      signal_mzs.insert(p.mz);
    }

    std::sort(all_peaks.begin(), all_peaks.end());
    boost::dynamic_bitset<> is_signal_bitset(all_peaks.size());
    // to speed up, use bitset
    for (Size i = 0; i < all_peaks.size(); i++)
    {
      if (signal_mzs.find(all_peaks[i].mz) != signal_mzs.end())
      {
        is_signal_bitset[i] = true;
      }
    }

    float charge_noise_pwr = 0;

    std::vector<std::vector<Size>> per_bin_edges(max_bin_number);
    std::vector<int> per_bin_start_index(max_bin_number, -2); // -2 means bin is empty. -1 means bin is used. zero or positive = edge index
    std::map<float, Size> max_intensity_sum_to_bin;

    for (Size k = 0; k < max_bin_number; k++)
    {
      per_bin_edges[k] = std::vector<Size>(all_peaks.size(), 0);
    }

    // first collect all possible edges. An edge means mass difference between two peaks.
    for (Size i = 0; i < all_peaks.size(); i++)
    {
      auto p1 = all_peaks[i];
      bool p1_signal = is_signal_bitset[i];
      for (Size j = i + 1; j < all_peaks.size(); j++)
      {
        auto p2 = all_peaks[j];
        double normalized_dist = (p2.mz - p1.mz) * z / iso_da_distance_;
        if (normalized_dist > .9 && normalized_dist < 1.1) // if the distance is too close to the isotope distance, skip
        {
          continue;
        }
        Size bin = (Size)round(normalized_dist * (max_bin_number - 5));
        if (bin == 0)
        {
          continue;
        }
        if (bin >= max_bin_number)
        {
          break;
        }

        if (p1_signal && is_signal_bitset[j] && normalized_dist >= .75) // if both are signals and they are different from each other by ~ one isotope distance, do not connect
        {
          continue;
        }
        per_bin_edges[bin][i] = j;
        per_bin_start_index[bin] = -1;
      }
    }

    // then from each bin find the highest intensity path consisting of the same mass differences.
    for (Size k = 0; k < max_bin_number; k++)
    {
      if (per_bin_start_index[k] == -2)
      {
        continue;
      }
      auto edges = per_bin_edges[k];
      float max_sum_intensity = 0;
      for (Size i = 0; i < edges.size(); i++)
      {
        if (edges[i] == 0)
        {
          continue;
        }
        float intensity = is_signal_bitset[i] ? 0 : all_peaks[i].intensity;
        float sum_intensity = intensity;
        Size j = edges[i];

        int cntr = 0; // how many edges?
        while (j < edges.size())
        {
          cntr++;
          j = edges[j];
          if (j <= 0)
          {
            break;
          }
          sum_intensity += intensity;
          if (!is_signal_bitset[j])
            intensity = all_peaks[j].intensity;
        }

        if (cntr > 2 && max_sum_intensity < sum_intensity) // at least two edges should be there.
        {
          max_sum_intensity = sum_intensity;
          per_bin_start_index[k] = (int)i;
        }
      }
      max_intensity_sum_to_bin[max_sum_intensity] = k;
    }

    auto used = boost::dynamic_bitset<>(all_peaks.size());

    // Now from the highest intensity path to the lowest, sum up intensities excluding already used peaks or signal peaks.
    for (auto it = max_intensity_sum_to_bin.rbegin(); it != max_intensity_sum_to_bin.rend(); ++it)
    {
      Size bin = it->second;
      int index = per_bin_start_index[bin];

      if (index < 0)
      {
        continue;
      }

      auto edges = per_bin_edges[bin];
      float intensity = is_signal_bitset[index] ? 0 : all_peaks[index].intensity;
      float sum_intensity = .0, sum_squared_intensity = .0;
      int skiped_peak_cntr = 0;

      if (!used[index])
      {
        sum_intensity += intensity;
      }
      else
      {
        sum_squared_intensity += intensity * intensity;
        skiped_peak_cntr++;
      }
      used[index] = true;

      Size j = edges[index];

      while (j < edges.size())
      {
        j = edges[j];
        if (j <= 0)
        {
          break;
        }

        if (!used[j])
        {
          sum_intensity += intensity;
        }
        else
        {
          sum_squared_intensity += intensity * intensity;
          skiped_peak_cntr++;
        }

        used[j] = true;
        if (!is_signal_bitset[j])
          intensity = all_peaks[j].intensity;
      }

      if (skiped_peak_cntr < 2)
      {
        charge_noise_pwr += sum_intensity * sum_intensity;
      }
      else
      {
        charge_noise_pwr += sum_squared_intensity;
      }
    }

    // if still peaks are remaining, add their individual power.
    for (Size i = 0; i < all_peaks.size(); i++)
    {
      if (used[i] || is_signal_bitset[i])
      {
        continue;
      }
      charge_noise_pwr += all_peaks[i].intensity * all_peaks[i].intensity;
    }
    return charge_noise_pwr;
  }

  void PeakGroup::updatePerChargeInformation_(const std::vector<LogMzPeak>& noisy_peaks)
  {
    per_charge_noise_pwr_ = std::vector<float>(1 + max_abs_charge_, .0);
    per_charge_sum_signal_squared_ = std::vector<float>(1 + max_abs_charge_, .0);
    per_charge_int_ = std::vector<float>(1 + max_abs_charge_, .0);

    // calculate per charge intensity, and per charge sum of signal intensity squared
    for (const auto& p : logMzpeaks_)
    {
      per_charge_int_[p.abs_charge] += p.intensity;
      per_charge_sum_signal_squared_[p.abs_charge] += p.intensity * p.intensity;
    }

    // for each charge calculate signal and noise power
    for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
    {
      per_charge_noise_pwr_[z] = getNoisePeakPower_(noisy_peaks, logMzpeaks_, z);
    }
  }

  void PeakGroup::updateChargeRange_(std::vector<LogMzPeak>& noisy_peaks)
  {
    int max_sig_charge = 0;
    float max_sig = 0;

    // first, find the maximum snr charge.
    for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
    {
      double tmp_snr = per_charge_int_[z] * per_charge_int_[z] / (1 + per_charge_noise_pwr_[z]);
      if (max_sig < tmp_snr)
      {
        max_sig = tmp_snr;
        max_sig_charge = z;
      }
    }

    // determine the final charge ranges based on per charge power.
    // If more than two consecutive charges do not contain any signal peak, the charge range stops at that charge.

    int new_max_abs_charge;
    int new_min_abs_charge;

    new_max_abs_charge = new_min_abs_charge = max_sig_charge;
    float threshold = std::min(max_sig / 10, 1.0f);
    for (int z = max_sig_charge; z <= max_abs_charge_; z++)
    {
      float per_charge_signal_power = per_charge_int_[z] * per_charge_int_[z];
      if ((per_charge_signal_power / (1 + per_charge_noise_pwr_[z])) < threshold)
      {
        break;
      }
      new_max_abs_charge = z;
    }

    for (int z = max_sig_charge; z >= min_abs_charge_; z--)
    {
      float per_charge_signal_power = per_charge_int_[z] * per_charge_int_[z];
      if ((per_charge_signal_power / (1 + per_charge_noise_pwr_[z])) < threshold)
      {
        break;
      }
      new_min_abs_charge = z;
    }

    // if the updated charge range is different from the original one, signal and noisy peaks are again updated
    if (max_abs_charge_ != new_max_abs_charge || min_abs_charge_ != new_min_abs_charge)
    {
      std::vector<LogMzPeak> new_logMzpeaks;
      new_logMzpeaks.reserve(size());
      std::vector<LogMzPeak> new_noisy_peaks;
      new_noisy_peaks.reserve(noisy_peaks.size());

      // now only take the signal and noise peaks within the new charge range.
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
    if (min_abs_charge_ > max_abs_charge_)
    {
      clear_();
    }
    else
    {
      sort();
    }
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak> PeakGroup::recruitAllPeaksInSpectrum(const MSSpectrum& spec, const double tol, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                                                                        double mono_mass)
  {
    const double mul_tol = .8; // not all peaks within tolerance are considered as signal peaks.
    std::vector<LogMzPeak> noisy_peaks;
    if (mono_mass < 0)
    {
      return noisy_peaks;
    }

    monoisotopic_mass_ = mono_mass;

    // int iso_margin = 3; // how many isotopes do we want to scan before the monoisotopic mass?
    int max_isotope = (int)avg.getLastIndex(mono_mass);
    int min_isotope = (int)(avg.getApexIndex(mono_mass) - avg.getLeftCountFromApex(mono_mass) + min_negative_isotope_index_);
    min_isotope = std::max(min_negative_isotope_index_, min_isotope);

    clear_(); // clear logMzPeaks
    negative_iso_peaks_.clear();

    reserve((max_isotope) * (max_abs_charge_ - min_abs_charge_ + 1) * 2);
    noisy_peaks.reserve(max_isotope * (max_abs_charge_ - min_abs_charge_ + 1) * 2);

    // scan from the largest to the smallest charges and recruit the raw peaks for this monoisotopic_mass
    for (int c = max_abs_charge_; c >= min_abs_charge_; c--)
    {
      if (c <= 0)
      {
        break;
      }

      double cmz = (mono_mass) / c + FLASHDeconvHelperStructs::getChargeMass(is_positive_);
      double left_mz = (mono_mass - (1 - min_negative_isotope_index_) * iso_da_distance_) / c + FLASHDeconvHelperStructs::getChargeMass(is_positive_);
      Size index = spec.findNearest(left_mz * (1 - tol * mul_tol));
      double iso_delta = iso_da_distance_ / c;

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

        if (abs(pmz - cmz - iso_index * iso_delta) <= pmz * tol * mul_tol)
        {
          auto p = LogMzPeak(spec[index], is_positive_);
          p.isotopeIndex = iso_index;
          p.abs_charge = c;
          if (iso_index < 0)
          {
            negative_iso_peaks_.push_back(p);
          }
          else
          {
            push_back(p);
          }
        }
        else if (iso_index >= 0)
        {
          auto p = LogMzPeak(spec[index], is_positive_);
          p.isotopeIndex = iso_index;
          p.abs_charge = c;
          noisy_peaks.push_back(p);
        }
      }

      if (index >= spec.size())
      {
        break;
      }
    }

    return noisy_peaks;
  }

  void PeakGroup::updateChargeFitScoreAndChargeIntensities_(bool is_low_charge)
  {
    if (max_abs_charge_ == min_abs_charge_)
    {
      charge_score_ = 1;
      return;
    }
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
    if (max_index < 0)
    {
      charge_score_ = 0;
      return;
    }

    first_index = first_index < 0 ? 0 : first_index;
    float p = .0f;
    const float factor = .3;
    for (int c = max_index; c < last_index; c++)
    {
      float diff = per_charge_int_[c + 1] - per_charge_int_[c];
      if (diff > 0)
        p += diff;
      else if (!is_low_charge && diff < -(per_charge_int_[c]) * factor)
        p -= diff + (per_charge_int_[c]) * factor;
    }
    for (int c = max_index; c > first_index; c--)
    {
      float diff = per_charge_int_[c - 1] - per_charge_int_[c];
      if (diff > 0)
        p += diff;
      else if (!is_low_charge && diff < -(per_charge_int_[c]) * factor)
        p -= diff + (per_charge_int_[c]) * factor;
    }
    charge_score_ = std::max(.0f, 1.0f - p / summed_intensity);
  }

  void PeakGroup::updateMonoMassAndIsotopeIntensities()
  {
    if (logMzpeaks_.size() == 0)
    {
      return;
    }
    int max_isotope_index = 0;
    std::sort(logMzpeaks_.begin(), logMzpeaks_.end());
    for (const auto& p : logMzpeaks_)
    {
      max_isotope_index = max_isotope_index < p.isotopeIndex ? p.isotopeIndex : max_isotope_index;
    }

    per_isotope_int_ = std::vector<float>(max_isotope_index + 1 - min_negative_isotope_index_, .0f);
    intensity_ = .0;
    double nominator = .0;

    for (const auto& p : logMzpeaks_)
    {
      float pi = p.intensity;
      if (p.isotopeIndex < 0)
      {
        continue;
      }
      per_isotope_int_[p.isotopeIndex - min_negative_isotope_index_] += pi;
      nominator += pi * (p.getUnchargedMass() - p.isotopeIndex * iso_da_distance_);
      intensity_ += pi;
    }
    for (const auto& p : negative_iso_peaks_)
    {
      if (p.isotopeIndex - min_negative_isotope_index_ < 0)
      {
        continue;
      }
      per_isotope_int_[p.isotopeIndex - min_negative_isotope_index_] += p.intensity;
    }

    monoisotopic_mass_ = nominator / intensity_;
  }

  bool PeakGroup::isSignalMZ(const double mz, const double tol) const
  {
    for (const auto& p : logMzpeaks_)
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

  void PeakGroup::setRepAbsCharge(const int max_snr_abs_charge)
  {
    max_snr_abs_charge_ = max_snr_abs_charge;
  }

  void PeakGroup::setChargeScore(const float score)
  {
    charge_score_ = score;
  }

  void PeakGroup::setAvgPPMError(const float error)
  {
    avg_ppm_error_ = error;
  }

  void PeakGroup::setQscore(const double qscore)
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
      for (const auto& tmp_p : logMzpeaks_)
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
    return max_snr_abs_charge_;
  }

  double PeakGroup::getQscore() const
  {
    return qscore_;
  }

  double PeakGroup::getQscore2D() const
  {
    // if (qscore2D_ < 0)
    //   return qscore_;
    return std::max(qscore_, qscore2D_);
  }

  void PeakGroup::setFeatureIndex(uint findex)
  {
    findex_ = findex;
  }

  bool PeakGroup::isTargeted() const
  {
    return is_targeted_;
  }

  void PeakGroup::updateSNR_()
  {
    float cos_squared = isotope_cosine_score_ * isotope_cosine_score_;
    float signal = 0, noise = 0, sum_signal_squared = 0;
    per_charge_snr_ = std::vector<float>(1 + max_abs_charge_, .0);

    for (size_t c = min_abs_charge_; c < std::min(per_charge_sum_signal_squared_.size(), size_t(1 + max_abs_charge_)); ++c)
    {
      if (per_charge_cos_.size() > c)
      {
        float per_charge_cos_squared = per_charge_cos_[c] * per_charge_cos_[c];
        float nom = per_charge_cos_squared * per_charge_int_[c] * per_charge_int_[c];
        float denom = 1 + per_charge_noise_pwr_[c] + (1 - per_charge_cos_squared) * per_charge_sum_signal_squared_[c];

        per_charge_snr_[c] = denom <= 0 ? .0f : (nom / denom);
      }
      sum_signal_squared += per_charge_sum_signal_squared_[c];
      signal += per_charge_int_[c] * per_charge_int_[c];
      noise += per_charge_noise_pwr_[c];
    }

    per_charge_sum_signal_squared_.clear();
    per_charge_noise_pwr_.clear();
    float t_nom = cos_squared * signal;
    float t_denom = 1 + noise + (1 - cos_squared) * sum_signal_squared;

    snr_ = t_denom <= 0 ? .0f : (t_nom / t_denom);
  }

  float PeakGroup::getQvalue(PeakGroup::TargetDecoyType flag) const
  {
    if (flag == PeakGroup::TargetDecoyType::target)
    {
      return std::min(1.0f, getQvalue(PeakGroup::TargetDecoyType::charge_decoy) + getQvalue(PeakGroup::TargetDecoyType::noise_decoy) + getQvalue(PeakGroup::TargetDecoyType::isotope_decoy));
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
    if (abs_charge < 0 || per_charge_int_.empty() || (int)per_charge_int_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_int_[abs_charge];
  }

  bool PeakGroup::isPositive() const
  {
    return is_positive_;
  }

  PeakGroup::TargetDecoyType PeakGroup::getTargetDecoyType() const
  {
    return target_decoy_type_;
  }

  void PeakGroup::setTargetDecoyType(PeakGroup::TargetDecoyType index)
  {
    target_decoy_type_ = index;
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

  uint PeakGroup::getFeatureIndex() const
  {
    return findex_;
  }

  void PeakGroup::setQscore2D(double fqscore)
  {
    qscore2D_ = fqscore;
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

  FLASHDeconvHelperStructs::LogMzPeak& PeakGroup::back()
  {
    return logMzpeaks_.back();
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

  void PeakGroup::sort()
  {
    std::sort(logMzpeaks_.begin(), logMzpeaks_.end());
  }

  void PeakGroup::setQvalue(double q, PeakGroup::TargetDecoyType flag)
  {
    qvalue_[flag] = std::min(1.0, q);
  }

  void PeakGroup::calculateDLMatrices(const MSSpectrum& spec, double tol, const PrecalculatedAveragine& avg)
  {
    dl_matrices_.clear();
    std::vector<LogMzPeak> noisy_peaks = recruitAllPeaksInSpectrum(spec, tol, avg, getMonoMass());

    int center_z = std::distance(per_charge_snr_.begin(), std::max_element(per_charge_snr_.begin(), per_charge_snr_.end()));

    if (center_z - charge_range_for_DL_ / 2 > min_abs_charge_ || center_z + charge_range_for_DL_ / 2 < max_abs_charge_)
    {
      if (center_z - charge_range_for_DL_ / 2 < min_abs_charge_)
      {
        center_z -= std::min(center_z - charge_range_for_DL_ / 2 - min_abs_charge_, center_z + charge_range_for_DL_ / 2 - max_abs_charge_);
      }
      else if (center_z + charge_range_for_DL_ / 2 > max_abs_charge_)
      {
        center_z += std::min(min_abs_charge_ - center_z + charge_range_for_DL_ / 2, max_abs_charge_ - center_z - charge_range_for_DL_ / 2);
      }
    }

    int min_z = center_z - charge_range_for_DL_ / 2;
    int max_z = center_z + charge_range_for_DL_ / 2;

    auto iso = avg.get(getMonoMass());
    int apex_iso = avg.getApexIndex(getMonoMass());
    int min_iso = apex_iso - (iso_range_for_DL_ / 2);
    int max_iso = apex_iso + (iso_range_for_DL_ / 2); // inclusive

    max_iso = std::min(max_iso, (int)iso.size() - 1);

    Matrix<float> sig, distortion, noise;

    sig.resize(charge_range_for_DL_, iso_range_for_DL_ / bin_width_DL_, .0);
    distortion.resize(charge_range_for_DL_, iso_range_for_DL_ / bin_width_DL_, .0);
    noise.resize(charge_range_for_DL_, iso_range_for_DL_ / bin_width_DL_, .0);

    std::vector<float> iso_vector(max_iso - min_iso + 1, 0);

    float sum = .0;
    for (int i = min_iso; i <= max_iso; i++)
    {
      if (i < 0)
        continue;
      sum += iso[i].getIntensity() * iso[i].getIntensity();
    }

    for (int i = min_iso; i <= max_iso; i++)
    {
      if (i < 0)
        continue;
      iso_vector[i - min_iso] = iso[i].getIntensity() / sqrt(sum);
    }
    std::vector<LogMzPeak> new_noisy_peaks;
    new_noisy_peaks.reserve(noisy_peaks.size());

    for (int z = min_z; z <= max_z; z++)
    {
      if (z <= 0)
        continue;
      std::vector<float> observed_iso_vector(iso_vector.size(), 0);
      std::vector<float> diff_iso_vector(iso_vector);
      float mul_factor = 0;

      for (const auto& p : logMzpeaks_)
      {
        if (p.abs_charge != z)
          continue;
        int index = p.isotopeIndex;

        if (index < min_iso || index > max_iso)
          continue;
        observed_iso_vector[index - min_iso] += p.intensity;
      }

      for (const auto& p : noisy_peaks)
      {
        if (p.abs_charge != z)
          continue;
        if (p.isotopeIndex < 0)
          continue;
        if (std::abs(monoisotopic_mass_ - p.getUnchargedMass() + p.isotopeIndex * iso_da_distance_) < bin_width_DL_)
        {
          int index = p.isotopeIndex;
          if (index < min_iso || index > max_iso)
            continue;
          observed_iso_vector[index - min_iso] += p.intensity;
        }
        else
        {
          new_noisy_peaks.push_back(p);
        }
      }

      for (int i = min_iso; i <= max_iso; i++)
      {
        mul_factor += diff_iso_vector[i - min_iso] * observed_iso_vector[i - min_iso];
      }
      for (int i = min_iso; i <= max_iso; i++)
      {
        diff_iso_vector[i - min_iso] *= mul_factor;
      }
      for (int i = min_iso; i <= max_iso; i++)
      {
        diff_iso_vector[i - min_iso] -= observed_iso_vector[i - min_iso];
      }
      for (int i = min_iso; i <= max_iso; i++)
      {
        int iso_bin = (int)round((i - min_iso) / bin_width_DL_);
        if (z - min_z >= sig.rows())
          break;
        if (iso_bin >= sig.cols())
          break;
        if (iso_bin < 0 || z - min_z < 0)
          continue;
        sig.setValue(z - min_z, iso_bin, observed_iso_vector[i - min_iso]);
        distortion.setValue(z - min_z, iso_bin, (diff_iso_vector[i - min_iso]));
      }
    }

    sum = .0;
    for (auto v : sig.asVector())
    {
      sum += v * v;
    }

    float normalization_factor = sqrt(sum);

    noisy_peaks = new_noisy_peaks;

    for (const auto& p : noisy_peaks)
    {
      if (p.isotopeIndex < 0)
        continue;
      int iso_bin = (int)round(((p.getUnchargedMass() - monoisotopic_mass_) / iso_da_distance_ - min_iso) / bin_width_DL_);
      int z = p.abs_charge;
      if (z - min_z >= noise.rows() || iso_bin >= noise.cols() || iso_bin < 0 || z - min_z < 0)
        continue;
      noise.setValue(z - min_z, iso_bin, (noise.getValue(z - min_z, iso_bin) + p.intensity));
    }

    for (const auto& p : negative_iso_peaks_)
    {
      int iso_bin = (int)round(((p.getUnchargedMass() - monoisotopic_mass_) / iso_da_distance_ - min_iso) / bin_width_DL_);
      int z = p.abs_charge;
      if (z - min_z >= noise.rows() || iso_bin >= noise.cols() || iso_bin < 0 || z - min_z < 0)
        continue;
      distortion.setValue(z - min_z, iso_bin, distortion.getValue(z - min_z, iso_bin) + p.intensity);
    }

    if (normalization_factor > 0)
    {
      for (Size r = 0; r < sig.rows(); r++)
      {
        for (Size c = 0; c < sig.cols(); c++)
        {
          sig.setValue(r, c, sig.getValue(r, c) / normalization_factor);
          distortion.setValue(r, c, distortion.getValue(r, c) / normalization_factor);
          noise.setValue(r, c, noise.getValue(r, c) / normalization_factor);
        }
      }
    }
    dl_matrices_.push_back(sig);
    dl_matrices_.push_back(distortion);
    dl_matrices_.push_back(noise);
  }

  Matrix<float> PeakGroup::getDLMatrix(int index) const
  {
    assert(index >= 0 && index < dl_matrices_.size());
    return dl_matrices_[index];
  }
} // namespace OpenMS
