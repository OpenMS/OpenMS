// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
    auto mass = (float)(monoisotopic_mass_ + p.isotopeIndex * iso_da_distance_);
    return (float)(abs(mass / (float)p.abs_charge + FLASHDeconvHelperStructs::getChargeMass(p.is_positive) - p.mz) / p.mz * 1e6);
  }

  float PeakGroup::getAbsDaError_(const LogMzPeak& p) const
  {
    auto mass = (float)(monoisotopic_mass_ + p.isotopeIndex * iso_da_distance_);
    return (float)(abs(mass - p.getUnchargedMass()));
  }

  int PeakGroup::getMinNegativeIsotopeIndex() const
  {
    return min_negative_isotope_index_;
  }

  void PeakGroup::updatePerChargeCos_(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg)
  {
    auto iso_dist = avg.get(monoisotopic_mass_);
    auto current_per_isotope_intensities = std::vector<float>(getIsotopeIntensities().size() + min_negative_isotope_index_, .0f);

    if (min_abs_charge_ == max_abs_charge_)
      setChargeIsotopeCosine(min_abs_charge_, getIsotopeCosine());
    else
    {
      for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
      {
        std::fill(current_per_isotope_intensities.begin(), current_per_isotope_intensities.end(), .0f);
        int min_isotope_index = (int)current_per_isotope_intensities.size();
        int max_isotope_index = -1; // this is inclusive!!

        for (const auto& peak : logMzpeaks_)
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

        float cos_score = SpectralDeconvolution::getCosine(current_per_isotope_intensities, min_isotope_index, max_isotope_index, iso_dist, 0, SpectralDeconvolution::min_iso_size,
                                                           target_decoy_type_ == PeakGroup::TargetDecoyType::noise_decoy);
        setChargeIsotopeCosine(abs_charge, cos_score); //
      }
    }
  }

  int PeakGroup::updateQscore(const std::vector<LogMzPeak>& noisy_peaks, const MSSpectrum& spec, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double min_cos, double tol,
                              bool is_low_charge, int allowed_iso_error, bool is_last)
  {
    qscore_ = 0;

    if (empty())
    {
      return 0;
    }
    updatePerChargeInformation_(noisy_peaks, tol, is_last);
    updateChargeRange_();
    updateChargeFitScoreAndChargeIntensities_(is_low_charge);
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
    int window_width = is_last ? 0 : -1;
    //auto target_decoy_type = is_last && target_decoy_type_ == PeakGroup::TargetDecoyType::isotope_decoy ? PeakGroup::TargetDecoyType::target : target_decoy_type_;
    allowed_iso_error = is_last? -1 : allowed_iso_error;
    isotope_cosine_score_ = SpectralDeconvolution::getIsotopeCosineAndIsoOffset(monoisotopic_mass_, per_isotope_int_, h_offset, avg,
                                                                                -min_negative_isotope_index_, // change if to select cosine calculation and if to get second best hits
                                                                                window_width, allowed_iso_error, target_decoy_type_);
    if (h_offset != 0)
      return h_offset;

    if (isotope_cosine_score_ < min_cos)
    {
      return 0;
    }
    updatePerChargeCos_(avg);
    updateAvgPPMError_();
    updateAvgDaError_();
    updateSNR_((float)avg.getSNRMultiplicationFactor(monoisotopic_mass_));

    for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
    {
      if (getChargeSNR(abs_charge) > getChargeSNR(max_snr_abs_charge_))
      {
        max_snr_abs_charge_ = abs_charge;
      }
    }

    qscore_ = Qscore::getQscore(this, spec);

    return h_offset;
  }

  float PeakGroup::getNoisePeakPower_(const std::vector<FLASHDeconvHelperStructs::LogMzPeak>& noisy_peaks, const int z, const double tol) const
  {
    if (noisy_peaks.empty())
      return 0;
    const Size max_noisy_peak_number = 40; // too many noise peaks will slow down the process
    const Size bin_number_margin = 8;
    const Size max_bin_number = bin_number_margin + 12; // 12 bin + 8 extra bin
    float threshold = -1;
    std::vector<std::pair<Peak1D, bool>> all_peaks; // peak + is signal?
    all_peaks.reserve(max_noisy_peak_number + logMzpeaks_.size());

    int noise_peak_count = 0;
    for (const auto& noisy_peak : noisy_peaks)
    {
      if (z > 0 && noisy_peak.abs_charge != z)
        continue;
      if (noisy_peak.abs_charge < min_abs_charge_ || noisy_peak.abs_charge > max_abs_charge_)
        continue;
      noise_peak_count++;
    }
    if (noise_peak_count == 0)
      return 0;
    // get intensity threshold
    if (noise_peak_count > (int)max_noisy_peak_number)
    {
      std::vector<float> intensities;
      intensities.reserve(noise_peak_count);
      for (const auto& noisy_peak : noisy_peaks)
      {
        if (z > 0 && noisy_peak.abs_charge != z)
          continue;
        if (noisy_peak.abs_charge < min_abs_charge_ || noisy_peak.abs_charge > max_abs_charge_)
          continue;
        intensities.push_back(noisy_peak.intensity);
      }

      std::sort(intensities.rbegin(), intensities.rend());
      threshold = intensities[max_noisy_peak_number];
    }

    // filter peaks and check which mzs are signal and which are noise.
    for (const auto& noisy_peak : noisy_peaks)
    {
      if ((z > 0 && noisy_peak.abs_charge != z) || noisy_peak.intensity < threshold)
        continue;
      if (noisy_peak.abs_charge < min_abs_charge_ || noisy_peak.abs_charge > max_abs_charge_)
        continue;
      all_peaks.emplace_back(Peak1D(noisy_peak.getUnchargedMass(), noisy_peak.intensity), false);
    }

    if (all_peaks.empty())
      return 0;

    for (const auto& peak : logMzpeaks_)
    {
      if ((z > 0 && peak.abs_charge != z) || peak.intensity < threshold)
        continue;
      if (peak.abs_charge < min_abs_charge_ || peak.abs_charge > max_abs_charge_)
        continue;
      all_peaks.emplace_back(Peak1D(peak.getUnchargedMass(), peak.intensity), true);
    }

    std::sort(all_peaks.begin(), all_peaks.end(), [](std::pair<Peak1D, bool>& left, std::pair<Peak1D, bool>& right) { return left.first.getMZ() < right.first.getMZ(); }); //

    float charge_noise_pwr = 0;

    std::vector<std::vector<Size>> per_bin_edges(max_bin_number);
    std::vector<int> per_bin_start_index(max_bin_number, -2); // -2 means bin is empty. -1 means bin is used. zero or positive = edge index
    std::map<float, Size> max_intensity_sum_to_bin;

    for (Size k = 0; k < max_bin_number; k++)
    {
      per_bin_edges[k] = std::vector<Size>(all_peaks.size(), 0);
    }
    // first collect all possible edges. An edge means mass difference between two peaks.
    const std::vector<double> div_factors {1.0, 2.0, 3.0}; // allow two skips for each bin
    for (Size i = 0; i < all_peaks.size(); i++)
    {
      const auto& [p1, p1_signal] = all_peaks[i];
      const auto p1_mass = p1.getMZ();
      std::vector<double> per_bin_error(max_bin_number, -1.0);
      for (Size j = i + 1; j < all_peaks.size(); j++)
      {
        const auto& [p2, p2_signal] = all_peaks[j];
        const double normalized_dist = (p2.getMZ() - p1_mass) / iso_da_distance_;

        if (p1_signal && p2_signal && normalized_dist >= .75) // if both are signals, and they are different from each other by more than .75 isotope distance, do not connect. Otherwise, connect as
                                                              // they may a part of consecutive other noisy peaks.
        {
          continue;
        }

        for (double d : div_factors)
        {
          double distance = normalized_dist / d * (max_bin_number - bin_number_margin);
          Size bin = (Size)round(distance);
          if (bin == 0)
          {
            continue;
          }
          if (bin >= max_bin_number)
          {
            break;
          }

          per_bin_start_index[bin] = -1;
          double current_error = d * max_bin_number + std::abs((double)bin - distance); // larger when d gets larger. For the same d, comparable.
          if (per_bin_error[bin] >= 0 && per_bin_error[bin] < current_error)
          {
            continue;
          }
          per_bin_edges[bin][i] = j;
          per_bin_error[bin] = current_error;
        }
      }
    }

    // then from each bin find the highest intensity path consisting of the same mass differences.
    for (Size k = 0; k < max_bin_number; k++)
    {
      if (per_bin_start_index[k] == -2)
      {
        continue;
      }
      const auto& edges = per_bin_edges[k];
      float max_sum_intensity = 0;
      for (Size i = 0; i < edges.size(); i++)
      {
        if (edges[i] == 0)
        {
          break;
        }
        const auto& [p1, p1_signal] = all_peaks[i];
        float intensity = // p1_signal ? 0 :
          p1.getIntensity();
        float sum_intensity = intensity;

        Size j = edges[i];

        while (j < edges.size())
        {
          const auto& [p2, p2_signal] = all_peaks[j];
          intensity = p2.getIntensity();
          sum_intensity += intensity;

          j = edges[j];
          if (j == 0)
          {
            break;
          }
        }

        if (max_sum_intensity < sum_intensity) // at least two edges should be there.
        {
          max_sum_intensity = sum_intensity;
          per_bin_start_index[k] = (int)i; //
        }
      }
      max_intensity_sum_to_bin[max_sum_intensity] = k; // how to deal with profile peaks?
    }

    auto unused = boost::dynamic_bitset<>(all_peaks.size());
    unused.flip();
    // Now from the highest intensity path to the lowest, sum up intensities excluding already used peaks or signal peaks.
    for (auto it = max_intensity_sum_to_bin.rbegin(); it != max_intensity_sum_to_bin.rend(); ++it)
    {
      Size bin = it->second;
      int index = per_bin_start_index[bin];
      if (index < 0)
      {
        continue;
      }

      const auto& edges = per_bin_edges[bin];
      const double ori_mass = all_peaks[index].first.getMZ();
      const int ori_index = index;
      float sum_intensity = .0;

      while (index < (int)all_peaks.size())
      {
        const auto& [p, p_signal] = all_peaks[index];
        if (p.getMZ() - ori_mass > tol / 2.0 * p.getMZ())
          break;
        float intensity = // p_signal ? 0 :
          p.getIntensity();

        Size j = edges[index];
        if (j == 0)
        {
          break;
        }

        if (unused[index])
        {
          sum_intensity += intensity;
          unused[index] = false;
        }
        else
        {
          break;
        }

        while (j < edges.size())
        {
          const auto& [p2, p2_signal] = all_peaks[j];
          // if (!p2_signal)
          intensity = p2.getIntensity();

          if (unused[j])
          {
            sum_intensity += intensity;
            unused[j] = false;
          }
          else
          {
            break;
          }
          j = edges[j];
          if (j == 0)
          {
            break;
          }
        }
        index++;
      }

      index = ori_index - 1;
      while (index >= 0)
      {
        const auto& [p, p_signal] = all_peaks[index];
        if (ori_mass - p.getMZ() > tol / 2.0 * ori_mass)
          break;
        float intensity = // p_signal ? 0 :
          p.getIntensity();

        Size j = edges[index];
        if (j == 0)
        {
          break;
        }

        if (unused[index])
        {
          sum_intensity += intensity;
          unused[index] = false;
        }
        else
        {
          break;
        }

        while (j < edges.size())
        {
          const auto& [p2, p2_signal] = all_peaks[j];
          // if (!p2_signal)
          intensity = p2.getIntensity();

          if (unused[j])
          {
            sum_intensity += intensity;
            unused[j] = false;
          }
          else
          {
            break;
          }
          j = edges[j];
          if (j == 0)
          {
            break;
          }
        }
        index--;
      }
      charge_noise_pwr += sum_intensity * sum_intensity;
    }

    // if still peaks are remaining, add their individual power.
    Size index = unused.find_first();
    double prev_mass = .0;
    float tmp_int_sum = 0;
    while (index != unused.npos)
    {
      const auto& [p, p_signal] = all_peaks[index];
      if (!p_signal)
      {
        if (p.getMZ() - prev_mass > p.getMZ() * tol)
        {
          charge_noise_pwr += tmp_int_sum * tmp_int_sum;
          prev_mass = p.getMZ();
          tmp_int_sum = 0;
        }
        tmp_int_sum += p.getIntensity();
      }
      index = unused.find_next(index);
    }
    charge_noise_pwr += tmp_int_sum * tmp_int_sum;

    return charge_noise_pwr;
  }

  void PeakGroup::updatePerChargeInformation_(const std::vector<LogMzPeak>& noisy_peaks, const double tol, const bool is_last)
  {
    per_charge_sum_signal_squared_ = std::vector<float>(1 + max_abs_charge_, .0f);
    per_charge_int_ = std::vector<float>(1 + max_abs_charge_, .0);
    int max_iso = 0;

    // calculate per charge intensity, and per charge sum of signal intensity squared
    for (const auto& p : logMzpeaks_)
    {
      per_charge_int_[p.abs_charge] += p.intensity;
      max_iso = std::max(max_iso, p.isotopeIndex);
    }
    Matrix<float> per_charge_isotope_int(1 + max_abs_charge_, 1 + max_iso, .0f);
    for (const auto& p : logMzpeaks_)
    {
      float prev_v = per_charge_isotope_int.getValue(p.abs_charge, p.isotopeIndex);
      per_charge_isotope_int.setValue(p.abs_charge, p.isotopeIndex, prev_v + p.intensity);
    }

    for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
    {
      for (Size i = 0; i < per_charge_isotope_int.cols(); i++)
      {
        float v = per_charge_isotope_int.getValue(z,i);
        per_charge_sum_signal_squared_[0] += v * v;
        per_charge_sum_signal_squared_[z] += v * v;
      }
    }

    // for each charge calculate signal and noise power
    per_charge_noise_pwr_ = std::vector<float>(1 + max_abs_charge_, .0f);

    if (is_last)
    {
      for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
      {
        per_charge_noise_pwr_[z] = getNoisePeakPower_(noisy_peaks, z, tol);
      }
      per_charge_noise_pwr_[0] = getNoisePeakPower_(noisy_peaks, 0, tol);
    }
    else
    {
      for (const auto& p : noisy_peaks)
      {
        float pwr = p.intensity * p.intensity;
        per_charge_noise_pwr_[0] += pwr;
        if (p.abs_charge < min_abs_charge_ || p.abs_charge > max_abs_charge_)
          continue;
        per_charge_noise_pwr_[p.abs_charge] += pwr;
      }
    }
  }

  void PeakGroup::updateChargeRange_()
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

      // now only take the signal and noise peaks within the new charge range.
      for (const auto& p : logMzpeaks_)
      {
        if (p.abs_charge < new_min_abs_charge || p.abs_charge > new_max_abs_charge)
        {
          continue;
        }
        new_logMzpeaks.push_back(p);
      }

      new_logMzpeaks.swap(logMzpeaks_);
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
    clear_(); // clear logMzPeaks

    if (mono_mass < 0)
    {
      return noisy_peaks;
    }

    if (max_abs_charge_ - min_abs_charge_ < max_abs_charge_ / 20) // if charge range is too small ...
    {
      return noisy_peaks;
    }

    monoisotopic_mass_ = mono_mass;

    // int iso_margin = 3; // how many isotopes do we want to scan before the monoisotopic mass?
    int max_isotope = (int)avg.getLastIndex(mono_mass);
    int min_isotope = (int)(avg.getApexIndex(mono_mass) - avg.getLeftCountFromApex(mono_mass) + min_negative_isotope_index_);
    int max_signal_isotope = 0;
    min_isotope = std::max(min_negative_isotope_index_, min_isotope);

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
            max_signal_isotope = std::max(max_signal_isotope, iso_index);
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

    std::vector<LogMzPeak> _noisy_peaks;
    _noisy_peaks.reserve(noisy_peaks.size());
    for (const auto& p : noisy_peaks)
    {
      if (p.isotopeIndex > max_signal_isotope)
        continue;
      _noisy_peaks.push_back(p);
    }

    return _noisy_peaks;
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

  /*
  float PeakGroup::getSumIntensityForMonoMass(float mono_mass, int hc, float iso_da, double tol) const
  {
    float sum = .0;
    for (auto& p : logMzpeaks_)
    {
      float mass_diff = (p.mass - mono_mass) * hc;
      float mass_error = abs(mass_diff - round(mass_diff/iso_da) * iso_da);
      if (mass_error < mono_mass * tol  * hc)
      {
        sum+=p.intensity;
      }
    }

    for (auto& p : noisy_peaks_)
    {
      float mass_diff = (p.mass - mono_mass) * hc;
      float mass_error = abs(mass_diff - round(mass_diff/iso_da) * iso_da);
      if (mass_error < mono_mass * tol  * hc)
      {
        sum+=p.intensity;
      }
    }
    return sum;
  }
*/

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

  float PeakGroup::getPeakOccupancy() const
  {
    int min_i = -1, max_i = 0;
    for (const auto& p : *this)
    {
      int i = p.isotopeIndex;
      max_i = std::max(max_i, i);
      if (min_i < 0)
        min_i = i;

      min_i = std::min(min_i, i);
    }

    auto used = std::vector<bool>((max_abs_charge_ - min_abs_charge_ + 1) * (max_i - min_i + 1), false);
    for (const auto& p : *this)
    {
      used[(p.abs_charge - min_abs_charge_ + 1) * (p.isotopeIndex - min_i + 1) - 1] = true;
    }
    int count = 0;
    for (const auto& b : used)
      if (b)
        count++;

    return (float)count / used.size();
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

  void PeakGroup::updateSNR_(float mul_factor)
  {
    per_charge_snr_ = std::vector<float>(1 + max_abs_charge_, .0);
    float total_nom = 1e-6;
    float total_denom = 1e-6;
    for (size_t c = min_abs_charge_; (int)c < 1 + max_abs_charge_; ++c)
    {
      if (per_charge_cos_.size() > c)
      {
        float per_charge_cos_squared = per_charge_cos_[c] * per_charge_cos_[c];
        float sig_pwr = per_charge_sum_signal_squared_[c] * per_charge_cos_squared;
        float nom = 1e-6f + mul_factor * sig_pwr;
        float denom = 1e-6f + per_charge_noise_pwr_[c] + (1 - per_charge_cos_squared) * per_charge_sum_signal_squared_[c];

        per_charge_snr_[c] = denom <= 0 ? .0f : (nom / denom);

        total_denom += (1 - per_charge_cos_squared) * per_charge_sum_signal_squared_[c];
        total_nom += sqrt(sig_pwr);
      }
    }

    snr_ = mul_factor * total_nom * total_nom / (total_denom + per_charge_noise_pwr_[0]);

    per_charge_sum_signal_squared_.clear();
    per_charge_noise_pwr_.clear();
  }

  float PeakGroup::getQvalue() const
  {
    return qvalue_;
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

  void PeakGroup::setQvalue(double q)
  {
    qvalue_ = q;
  }
} // namespace OpenMS
