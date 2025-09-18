// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

  int PeakGroup::getMinNegativeIsotopeIndex() const
  {
    return min_negative_isotope_index_;
  }

  void PeakGroup::updatePerChargeCos_(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg)
  {
    auto iso_dist = avg.get(monoisotopic_mass_);
    int iso_size = (int)iso_dist.size();
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

      float cos_score = FLASHDeconvAlgorithm::getCosine(current_per_isotope_intensities, min_isotope_index, max_isotope_index + 1, iso_dist, iso_size, 0, 0);
      setChargeIsotopeCosine(abs_charge, cos_score); //
    }
  }

  int PeakGroup::updateQscore(std::vector<LogMzPeak>& noisy_peaks, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double min_cos, int allowed_iso_error)
  {
    qscore_ = 0;
    updatePerChargeInformation_(noisy_peaks);
    updateChargeRange_(noisy_peaks);

    if (empty())
    {
      return 0;
    }

    updateChargeFitScoreAndChargeIntensities_();
    if (charge_score_ < .7f) //
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
      FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(monoisotopic_mass_, per_isotope_int_, h_offset, avg, -min_negative_isotope_index_, -1, allowed_iso_error, target_dummy_type_);

    if (isotope_cosine_score_ < min_cos)
    {
      return 0;
    }

    if (max_abs_charge_ - min_abs_charge_ < max_abs_charge_ / 20) // if charge range is too small ...
    {
      return 0;
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

      float q_score = Qscore::getQscore(this);
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

  float PeakGroup::getNoisePeakPower_(const std::vector<FLASHDeconvHelperStructs::LogMzPeak>& noisy_peaks, const std::vector<FLASHDeconvHelperStructs::LogMzPeak>& signal_peaks) const
  {
    const Size max_noisy_peaks = 50; // too many noise peaks will slow down the process
    const Size max_bin_number = 29;  // 24 bin + 5 extra bin
    float threshold = 0;

    std::vector<FLASHDeconvHelperStructs::LogMzPeak> all_peaks;
    std::set<double> signal_mzs;

    all_peaks.reserve(noisy_peaks.size() + signal_peaks.size());
    int z = 0;

    // get intensity threshold
    if (noisy_peaks.size() + signal_peaks.size() > max_noisy_peaks)
    {
      std::vector<float> intensities(noisy_peaks.size() + signal_peaks.size(), .0f);
      for (const auto& noisy_peak : noisy_peaks)
      {
        intensities.push_back(noisy_peak.intensity);
      }

      std::sort(intensities.begin(), intensities.end());
      threshold = intensities[intensities.size() - max_noisy_peaks];
    }

    // filter peaks and check which mzs are signal and which are noise.
    for (auto& p : noisy_peaks)
    {
      if (p.intensity < threshold)
      {
        continue;
      }
      all_peaks.push_back(p);
    }

    for (auto& p : signal_peaks)
    {
      z = p.abs_charge;
      if (p.intensity < threshold)
      {
        continue;
      }
      all_peaks.push_back(p);
      signal_mzs.insert(p.mz);
    }

    if (z == 0)
    {
      return 0;
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
    for (auto& p : logMzpeaks_)
    {
      per_charge_int_[p.abs_charge] += p.intensity;
      per_charge_sum_signal_squared_[p.abs_charge] += p.intensity * p.intensity;
    }

    std::vector<LogMzPeak> charge_noisy_peaks;
    std::vector<LogMzPeak> charge_signal_peaks;

    // for each charge calculate signal and noise power
    for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
    {
      charge_noisy_peaks.clear();
      charge_signal_peaks.clear();
      charge_noisy_peaks.reserve(noisy_peaks.size());
      charge_signal_peaks.reserve(size());
      for (auto& p : noisy_peaks)
      {
        if (p.abs_charge != z)
        {
          continue;
        }
        charge_noisy_peaks.push_back(p);
      }

      for (auto& p : logMzpeaks_)
      {
        if (p.abs_charge != z)
        {
          continue;
        }
        charge_signal_peaks.push_back(p);
      }
      per_charge_noise_pwr_[z] = getNoisePeakPower_(charge_noisy_peaks, charge_signal_peaks);
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
                                                                                        double mono_mass, const std::unordered_set<double>& excluded_peak_mzs)
  {
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
      Size index = spec.findNearest(left_mz * (1 - tol));
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
        // if excluded_peak_mzs_ is not empty, these mzs should be ignored in this raw spectrum for this peak group! But they can be included in noisy peaks.
        bool excluded = excluded_peak_mzs.size() > 0 && excluded_peak_mzs.find(pmz) != excluded_peak_mzs.end();

        if (!excluded && abs(pmz - cmz - iso_index * iso_delta) <= pmz * tol)
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

  void PeakGroup::updateChargeFitScoreAndChargeIntensities_()
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
    float min_intensity = -1;
    for (int c = min_abs_charge_; c <= max_abs_charge_; c++)
    {
      if (min_intensity < 0 || min_intensity > per_charge_int_[c])
      {
        min_intensity = per_charge_int_[c];
      }
    }
    for (int c = min_abs_charge_; c <= max_abs_charge_; c++)
    {
      summed_intensity += per_charge_int_[c] - min_intensity;
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

  void PeakGroup::updateMonoMassAndIsotopeIntensities()
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

    per_isotope_int_ = std::vector<float>(max_isotope_index + 1 - min_negative_isotope_index_, .0f);
    intensity_ = .0;
    double nominator = .0;

    for (auto& p : logMzpeaks_)
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
    for (auto& p : negative_iso_peaks_)
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
    return max_snr_abs_charge_;
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

  void PeakGroup::sort()
  {
    std::sort(logMzpeaks_.begin(), logMzpeaks_.end());
  }

  void PeakGroup::setQvalue(float q, PeakGroup::TargetDummyType flag)
  {
    qvalue_[flag] = q;
  }
} // namespace OpenMS
