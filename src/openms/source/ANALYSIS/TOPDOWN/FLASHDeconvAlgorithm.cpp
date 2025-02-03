// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  /// harmonic charge factors that will be considered for harmonic mass reduction.
  inline const std::vector<int> harmonic_charges_ {2, 3, 5, 7, 11};
  /// high and low charges are differently deconvolved. This value determines the (inclusive) threshold for low charge.
  inline const int low_charge_ = 10; // 10 inclusive
  inline const double tol_div_factor = 2.5; // use narrow tolerance for deconvolution and at the end use the input tolerance to filter out overlapping masses.

  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm() : DefaultParamHandler("FLASHDeconvAlgorithm")
  {
    defaults_.setValue("tol", DoubleList {10.0, 10.0}, "ppm tolerance for MS1, 2, ... (e.g., -tol 10.0 5.0 to specify 10.0 and 5.0 ppm for MS1 and MS2, respectively)");

    defaults_.setValue("min_mass", 50.0, "Minimum mass (Da)");
    defaults_.setValue("max_mass", 100000.0, "Maximum mass (Da)");

    defaults_.setValue("min_charge", 2, "Minimum charge state for MS1 spectra (can be negative for negative mode)");
    defaults_.setValue("max_charge", 100, "Maximum charge state for MS1 spectra (can be negative for negative mode)");

    defaults_.setValue("min_mz", -1.0, "If set to positive value, minimum m/z to deconvolve.");
    defaults_.setValue("max_mz", -1.0, "If set to positive value, maximum m/z to deconvolve.");
    defaults_.setValue("min_rt", -1.0, "If set to positive value, minimum RT to deconvolve.");
    defaults_.setValue("max_rt", -1.0, "If set to positive value, maximum RT to deconvolve.");

    defaults_.setValue("isolation_window", 5.0, "Default isolation window with. If the input mzML file does not contain isolation window width information, this width will be used.");
    defaults_.addTag("isolation_window", "advanced");

    defaults_.setValue("min_isotope_cosine", DoubleList {.85, .85},
                       "Cosine threshold between avg. and observed isotope pattern for MS1, 2, ... (e.g., -min_isotope_cosine_ 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)");

    defaults_.setValue("allowed_isotope_error", 1,
                       "Allowed isotope index error for decoy and qvalue report. If it is set to 1, for example, +-1 isotope errors are "
                       "not counted as false. Beta version.");
    defaults_.addTag("allowed_isotope_error", "advanced");

    defaults_.setValue("min_intensity", 0.0, "Intensity threshold");
    defaultsToParam_();
  }

  // Calculate the nominal (integer) mass from double mass. Multiplying 0.999497 to the original mass and then rounding reduce the error between the original and nominal masses.
  int FLASHDeconvAlgorithm::getNominalMass(const double mass)
  {
    return (int)(mass * 0.999497 + .5);
  }

  void FLASHDeconvAlgorithm::addMZsToExcludsionList(const DeconvolvedSpectrum& dspec, std::unordered_set<double>& excluded_mzs)
  {
    for (auto& pg : dspec)
    {
      for (auto& p : pg)
      {
        excluded_mzs.insert(p.mz);
      }
    }
  }

  // The main function called from outside. precursor_map_for_FLASHIda is used to read FLASHIda information
  void FLASHDeconvAlgorithm::performSpectrumDeconvolution(const MSSpectrum& spec, const std::vector<DeconvolvedSpectrum>& survey_scans, const int scan_number,
                                                          const std::map<int, std::vector<std::vector<float>>>& precursor_map_for_FLASHIda)
  {
    // First prepare for decoy runs.
    iso_da_distance_ =
      target_dummy_type_ != PeakGroup::noise_dummy ?
        Constants::ISOTOPE_MASSDIFF_55K_U :
        Constants::ISOTOPE_MASSDIFF_55K_U * sqrt(7.0) / 2.0; // sqrt(7.0)/2.0 Da is used instead of C13 - C12 to make sure masses detected with this nonsensical mass difference are not true.
    previously_deconved_mono_masses_for_dummy_.clear();
    previously_deconved_mass_bins_for_dummy_.reset();
    excluded_peak_mzs_.clear();

    if (target_dummy_type_ == PeakGroup::charge_dummy) // charge decoy
    {
      for (auto& pg : *target_dspec_for_dummy_calcualtion_)
      {
        int min_iso = -1, max_iso = 0;
        for (auto& p : pg)
        {
          previously_deconved_mono_masses_for_dummy_.push_back(p.getUnchargedMass());
          min_iso = min_iso < 0 ? p.isotopeIndex : std::min(min_iso, p.isotopeIndex);
          max_iso = std::max(max_iso, p.isotopeIndex);
        }
      }
    }
    if (target_dummy_type_ == PeakGroup::noise_dummy) // noise decoy
    {
      addMZsToExcludsionList(*target_dspec_for_dummy_calcualtion_, excluded_peak_mzs_);
    }

    ms_level_ = spec.getMSLevel();
    deconvolved_spectrum_ = DeconvolvedSpectrum(scan_number);
    deconvolved_spectrum_.setOriginalSpectrum(spec);

    // for MSn (n>1) register precursor peak and peak group.
    registerPrecursor_(survey_scans, precursor_map_for_FLASHIda);

    // rt range of analysis
    if (min_rt_ > 0 && spec.getRT() < min_rt_)
    {
      return;
    }
    if (max_rt_ > 0 && spec.getRT() > max_rt_)
    {
      return;
    }
    // based on MS level, adjust charge and mass ranges. Precursor charge and mass determine those.
    current_max_charge_ = deconvolved_spectrum_.getCurrentMaxAbsCharge(max_abs_charge_); //
    current_max_mass_ = deconvolved_spectrum_.getCurrentMaxMass(max_mass_);
    current_min_mass_ = deconvolved_spectrum_.getCurrentMinMass(min_mass_);

    // set universal pattern filter and harmonic pattern filters
    setFilters_();
    // LogMzPeaks are generated from raw peaks
    updateLogMzPeaks_();
    if (log_mz_peaks_.empty())
    {
      return;
    }

    // This is the main FLASHDeconv function in which deconvolution is performed.
    generatePeakGroupsFromSpectrum_();
  }

  void FLASHDeconvAlgorithm::updateMembers_()
  {
    min_mz_ = param_.getValue("min_mz");
    max_mz_ = param_.getValue("max_mz");

    min_rt_ = param_.getValue("min_rt");
    max_rt_ = param_.getValue("max_rt");

    min_abs_charge_ = param_.getValue("min_charge");
    max_abs_charge_ = param_.getValue("max_charge");
    is_positive_ = min_abs_charge_ > 0;

    min_abs_charge_ = abs(min_abs_charge_);
    max_abs_charge_ = abs(max_abs_charge_);

    if (min_abs_charge_ > max_abs_charge_)
    {
      int tmp = min_abs_charge_;
      min_abs_charge_ = max_abs_charge_;
      max_abs_charge_ = tmp;
    }

    max_mass_ = param_.getValue("max_mass");
    min_mass_ = param_.getValue("min_mass");

    isolation_window_size_ = param_.getValue("isolation_window");

    intensity_threshold_ = param_.getValue("min_intensity");

    bin_mul_factors_.clear();
    tolerance_ = param_.getValue("tol");

    for (double& j : tolerance_)
    {
      j *= 1e-6;
      j /= tol_div_factor; // finder bins are far better.
      bin_mul_factors_.push_back(1.0 / j);
    }

    min_isotope_cosine_ = param_.getValue("min_isotope_cosine");
    allowed_iso_error_ = param_.getValue("allowed_isotope_error");
  }

  const FLASHDeconvHelperStructs::PrecalculatedAveragine& FLASHDeconvAlgorithm::getAveragine()
  {
    return avg_;
  }

  void FLASHDeconvAlgorithm::calculateAveragine(const bool use_RNA_averagine)
  {
    CoarseIsotopePatternGenerator generator(300);

    auto iso = use_RNA_averagine ? generator.estimateFromRNAWeight(max_mass_) : generator.estimateFromPeptideWeight(max_mass_);
    iso.trimRight(0.0001 * iso.getMostAbundant().getIntensity());

    generator.setMaxIsotope(iso.size());
    avg_ = FLASHDeconvHelperStructs::PrecalculatedAveragine(50, max_mass_, 25, generator, use_RNA_averagine);
    avg_.setMaxIsotopeIndex((int)(iso.size() - 1));
  }

  // generate filters
  void FLASHDeconvAlgorithm::setFilters_()
  {
    filter_.clear();

    int charge_range = current_max_charge_;
    for (int i = 0; i < charge_range; i++)
    {
      filter_.push_back(-log(i + 1)); //+
    }

    harmonic_filter_matrix_.getEigenMatrix().resize(harmonic_charges_.size(), charge_range);
    harmonic_filter_matrix_.getEigenMatrix().setZero();

    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      int hc = harmonic_charges_[k];
      int n = hc / 2;

      for (int i = 0; i < charge_range; i++)
      {
        double a = i > 0 ? exp(-filter_[i - 1]) : 0;
        double b = exp(-filter_[i]);
        harmonic_filter_matrix_(k, i) = -log(b - (b - a) * n / hc);
      }
    }
  }

  // Generate uncharged log mz transformated peaks
  void FLASHDeconvAlgorithm::updateLogMzPeaks_()
  {
    log_mz_peaks_.clear();
    double threshold = intensity_threshold_;
    log_mz_peaks_.reserve(deconvolved_spectrum_.getOriginalSpectrum().size());

    // threshold = threshold < min_intensity * 2 ? min_intensity * 2 : threshold;
    for (auto& peak : deconvolved_spectrum_.getOriginalSpectrum())
    {
      if (peak.getIntensity() <= threshold) //
      {
        continue;
      }
      if (min_mz_ > 0 && peak.getMZ() < min_mz_)
      {
        continue;
      }
      if (max_mz_ > 0 && peak.getMZ() > max_mz_)
      {
        break;
      }
      if (!excluded_peak_mzs_.empty() && excluded_peak_mzs_.find(peak.getMZ()) != excluded_peak_mzs_.end())
      {
        continue;
      }

      LogMzPeak log_mz_peak(peak, is_positive_);
      log_mz_peaks_.push_back(log_mz_peak);
    }
  }

  // from bin to raw value
  double FLASHDeconvAlgorithm::getBinValue_(const Size bin, const double min_value, const double bin_mul_factor)
  {
    return min_value + (double)bin / bin_mul_factor;
  }

  // from value to bin number
  Size FLASHDeconvAlgorithm::getBinNumber_(const double value, const double min_value, const double bin_mul_factor)
  {
    if (value < min_value)
    {
      return 0;
    }
    return (Size)((value - min_value) * bin_mul_factor + .5);
  }

  // From log mz to mz bins.
  void FLASHDeconvAlgorithm::updateMzBins_(const Size bin_number, std::vector<float>& mz_bin_intensities)
  {
    mz_bins_ = boost::dynamic_bitset<>(bin_number);
    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];

    for (auto& p : log_mz_peaks_)
    {
      Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_mul_factor);
      if (bi >= bin_number)
      {
        continue;
      }
      mz_bins_.set(bi);

      mz_bin_intensities[bi] += p.intensity;
    }
  }

  // Find candidate mass bins from the current spectrum. The runtime of FLASHDeconv is determined by this function.
  void FLASHDeconvAlgorithm::updateCandidateMassBins_(std::vector<float>& mass_intensities, const std::vector<float>& mz_intensities)
  {
    Size mz_bin_index = mz_bins_.find_first();
    auto mz_bin_index_reverse = std::vector<Size>();
    mz_bin_index_reverse.reserve(mz_bins_.count());
    // invert mz bins so charges are counted from small to large given a mass
    while (mz_bin_index != mz_bins_.npos)
    {
      mz_bin_index_reverse.push_back(mz_bin_index);
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }
    int charge_range = current_max_charge_;
    size_t h_charge_size = harmonic_charges_.size();
    long bin_end = (long)mass_bins_.size();

    auto support_peak_count = std::vector<unsigned short>(mass_bins_.size(), 0); // per mass bin how many peaks are present

    // to calculate continuous charges, the previous charge value per mass should be stored
    auto prev_charges = std::vector<unsigned short>(mass_bins_.size(), charge_range + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prev_intensities = std::vector<float>(mass_bins_.size(), 1.0f);

    mass_intensities = std::vector<float>(mass_bins_.size(), .0f);

    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];
    std::vector<float> sub_max_h_intensity(h_charge_size, .0f);

    for (int i = (int)mz_bin_index_reverse.size() - 1; i >= 0; i--)
    {
      mz_bin_index = mz_bin_index_reverse[i];
      float intensity = mz_intensities[mz_bin_index];
      double mz = -1.0, log_mz = 0;

      log_mz = getBinValue_(mz_bin_index, mz_bin_min_value_, bin_mul_factor); // uncharged log mz
      mz = exp(log_mz);                                                       // uncharged mz
      // scan through charges
      for (int j = 0; j < charge_range; j++) // take all charge one ?
      {
        // mass is given by shifting by bin_offsets_[j]
        long mass_bin_index = (long)mz_bin_index + bin_offsets_[j];

        if (mass_bin_index < 0)
        {
          continue;
        }
        if (mass_bin_index >= bin_end)
        {
          break;
        }

        if (!previously_deconved_mono_masses_for_dummy_.empty() && previously_deconved_mass_bins_for_dummy_[mass_bin_index])
        {
          continue;
        }

        auto& spc = support_peak_count[mass_bin_index];
        int abs_charge = (j + 1);
        float& prev_intensity = prev_intensities[mass_bin_index];
        auto& prev_charge = prev_charges[mass_bin_index];
        bool charge_not_continous = prev_charge - j != -1 && (prev_charge <= charge_range);
        bool pass_first_check = false;

        // intensity ratio between consecutive charges should not exceed the factor.
        float factor = abs_charge <= low_charge_ ? 10.0f : (5.0f + 5.0f * low_charge_ / (float)abs_charge);
        // intensity ratio between consecutive charges for possible harmonic should be within this factor

        float hfactor = factor / 2.0f;
        // intensity of previous charge
        // intensity ratio between current and previous charges
        float intensity_ratio = intensity / prev_intensity;
        intensity_ratio = intensity_ratio < 1 ? 1.0f / intensity_ratio : intensity_ratio;
        float support_peak_intensity = 0;
        // check if peaks of continuous charges are present

        // if charge not continuous or intensity ratio is too high reset continuousChargePeakPairCount
        if (charge_not_continous || intensity_ratio > factor)
        {
          spc = 0;
        }
        else
        {
          pass_first_check = true;
          if (spc == 0 && abs_charge > low_charge_)
          {
            support_peak_intensity = prev_intensity;
          }
        }

        // for low charges, check isotope peak presence.
        float max_h_intensity = 0;

        if (!pass_first_check && abs_charge <= low_charge_)
        { // for low charges
          std::fill(sub_max_h_intensity.begin(), sub_max_h_intensity.end(), .0f);
          for (int d = 1; d >= -1; d -= 2)
          {
            bool iso_exist = false;
            double diff = d * iso_da_distance_ / abs_charge / mz;
            int next_iso_bin = 0;
            for (int t = -1; t < 2; t++)
            {
              int nib = (int)getBinNumber_(log_mz + diff, mz_bin_min_value_, bin_mul_factor) + t;
              if (nib != (int)mz_bin_index && nib > 0 && nib < (int)mz_bins_.size() && mz_bins_[nib])
              {
                iso_exist = true;
                pass_first_check = true;

                if (next_iso_bin == 0 || mz_intensities[next_iso_bin] < mz_intensities[nib])
                {
                  next_iso_bin = nib;
                }
              }
            }

            // harmonic check
            if (iso_exist)
            {
              double h_threshold = (intensity + mz_intensities[next_iso_bin]);

              for (size_t k = 0; k < h_charge_size; k++)
              {
                int hc = harmonic_charges_[k];
                int harmonic_cntr = 0;
                if (hc * abs_charge > current_max_charge_)
                {
                  break;
                }

                double hdiff = diff / hc;

                // check if there are harmonic peaks between the current peak and the next isotope peak.
                for (int t = -1; t < 2; t++)
                {
                  int next_harmonic_iso_bin = (int)getBinNumber_(log_mz + hdiff, mz_bin_min_value_, bin_mul_factor) + t;

                  // no perfect filtration. Just obvious ones are filtered out by checking if a peak is in the harmonic position and the intensity ratio is within two folds from the current peak
                  // (specified by mz_bin_index)
                  if (next_harmonic_iso_bin != (int)mz_bin_index && next_harmonic_iso_bin >= 0 && next_harmonic_iso_bin < (int)mz_bins_.size() && mz_bins_[next_harmonic_iso_bin] &&
                      mz_intensities[next_harmonic_iso_bin] > h_threshold / 2 && mz_intensities[next_harmonic_iso_bin] < h_threshold * 2)
                  {
                    harmonic_cntr++;
                    sub_max_h_intensity[k] = sub_max_h_intensity[k] < mz_intensities[next_harmonic_iso_bin] ? mz_intensities[next_harmonic_iso_bin] : sub_max_h_intensity[k];
                  }
                }

                if (harmonic_cntr > 0)
                {
                  pass_first_check = false;
                }
              }
              if (pass_first_check)
              {
                support_peak_intensity += mz_intensities[next_iso_bin];
              }
            }
          }
          max_h_intensity = *std::max_element(sub_max_h_intensity.begin(), sub_max_h_intensity.end());
          pass_first_check &= max_h_intensity <= 0;
        }

        if (pass_first_check)
        {
          if (prev_charge - j == -1) // check harmonic artifacts for high charge ranges
          {
            float max_intensity = intensity;
            float min_intensity = prev_intensity;
            if (prev_intensity <= 1.0)
            {
              max_intensity = intensity;
              min_intensity = intensity;
            }
            else if (min_intensity > max_intensity)
            {
              float tmpi = min_intensity;
              min_intensity = max_intensity;
              max_intensity = tmpi;
            }

            float high_threshold = max_intensity * hfactor;
            float low_threshold = min_intensity / hfactor;

            bool is_harmonic = false;

            // check if harmonic peaks are present with different harmonic multiple factors (2, 3, 5, 7  defined in harmonic_charges_).
            for (size_t k = 0; k < h_charge_size; k++)
            {
              for (int t = -1; t < 2; t++)
              {
                long hmz_bin_index = mass_bin_index - harmonic_bin_offset_matrix_(k, j) + t;

                if (hmz_bin_index > 0 && hmz_bin_index != (long)mz_bin_index && hmz_bin_index < (int)mz_bins_.size() && mz_bins_[hmz_bin_index])
                {
                  float harmonic_intensity = mz_intensities[hmz_bin_index];
                  if (harmonic_intensity > low_threshold && harmonic_intensity < high_threshold)
                  {
                    max_h_intensity = std::max(max_h_intensity, harmonic_intensity);
                    is_harmonic = true;
                  }
                }
              }
            }

            if (!is_harmonic) // if it is not harmonic
            {
              mass_intensities[mass_bin_index] += intensity + support_peak_intensity;
              spc++;
              if (spc >= min_support_peak_count_ || spc >= abs_charge / 2)
              {
                mass_bins_[mass_bin_index] = true;
              }
            }
            else // if harmonic
            {
              mass_intensities[mass_bin_index] -= max_h_intensity;
              if(spc > 0)
              {
                spc--;
              }
            }
          }
          else if (abs_charge <= low_charge_) // for low charge, include the mass if isotope is present
          {
            mass_intensities[mass_bin_index] += intensity + support_peak_intensity;
            spc++;
            mass_bins_[mass_bin_index] = true;
          }
        }
        prev_intensity = intensity;
        prev_charge = j;
      }
    }
  }

  // Subfunction of updateMassBins_. If a peak corresponds to multiple masses, only one mass is selected for the peak based on intensities.
  // mass level harmonic check is also performed in this function
  // it also outputs the charge range of each mass bin
  Matrix<int> FLASHDeconvAlgorithm::filterMassBins_(const std::vector<float>& mass_intensities)
  {
    int charge_range = current_max_charge_;
    Matrix<int> abs_charge_ranges(2, mass_bins_.size(), INT_MAX);
    for (Size i = 0; i < mass_bins_.size(); i++)
    {
      abs_charge_ranges(1, (int)i) = INT_MIN;
    }
    Size mz_bin_index = mz_bins_.find_first();
    long bin_size = (long)mass_bins_.size();

    auto to_skip = mass_bins_.flip();
    mass_bins_.reset();

    const int select_top_N = 3; // select top N charges per peak. We allow up to 3 just to consider frequent coelution
    std::vector<long> max_indices(select_top_N, -1);
    std::vector<int> max_intensity_abs_charge_ranges(select_top_N, -1);

    while (mz_bin_index != mz_bins_.npos)
    {
      std::fill(max_indices.begin(), max_indices.end(), -1);
      std::fill(max_intensity_abs_charge_ranges.begin(), max_intensity_abs_charge_ranges.end(), -1);

      float max_intensity = -1e11f;

      for (int j = 0; j < charge_range; j++)
      {
        long mass_bin_index = (long)mz_bin_index + bin_offsets_[j];

        if (mass_bin_index < 0)
        {
          continue;
        }
        if (mass_bin_index >= bin_size)
        {
          break;
        }

        if (!previously_deconved_mono_masses_for_dummy_.empty() && previously_deconved_mass_bins_for_dummy_[mass_bin_index])
        {
          continue;
        }

        if (!target_mono_masses_.empty() && target_mass_bins_[mass_bin_index])
        {
          float t = mass_intensities[mass_bin_index];

          if (t <= 0)
          { // no signal
            continue;
          }
          max_intensity = 1e38f;

          // store best values after shift by 1.
          for (int i = select_top_N - 1; i > 0; i--)
          {
            max_indices[i] = max_indices[i - 1];
            max_intensity_abs_charge_ranges[i] = max_intensity_abs_charge_ranges[i - 1];
          }

          max_indices[0] = mass_bin_index;
          max_intensity_abs_charge_ranges[0] = j;
        }
        else
        {
          if (to_skip[mass_bin_index])
          {
            continue;
          }

          float t = mass_intensities[mass_bin_index];
          if (t <= 0) // no signal
          {
            continue;
          }

          if (max_intensity < t)
          {
            max_intensity = t;
            // store best values after shift by 1.
            for (int i = select_top_N - 1; i > 0; i--)
            {
              max_indices[i] = max_indices[i - 1];
              max_intensity_abs_charge_ranges[i] = max_intensity_abs_charge_ranges[i - 1];
            }
            max_indices[0] = mass_bin_index;
            max_intensity_abs_charge_ranges[0] = j;
          }
        }
      }
      for (int i = 0; i < select_top_N; i++)
      {
        long max_index = max_indices[i];
        int max_intensity_abs_charge_range = max_intensity_abs_charge_ranges[i];
        if (max_index >= 0 && max_index < bin_size)
        {
          abs_charge_ranges(0, max_index) = std::min(abs_charge_ranges(0, max_index), max_intensity_abs_charge_range);
          abs_charge_ranges(1, max_index) = std::max(abs_charge_ranges(1, max_index), max_intensity_abs_charge_range);
          mass_bins_[max_index] = true;
        }
      }
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }

    return abs_charge_ranges;
  }

  // update mass bins which will be used to select peaks in the input spectrum...
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins_(const std::vector<float>& mz_intensities)
  {
    std::vector<float> mass_intensities;
    updateCandidateMassBins_(mass_intensities, mz_intensities);

    auto per_mass_abs_charge_ranges = filterMassBins_(mass_intensities);

    return per_mass_abs_charge_ranges;
  }

  // With mass_bins_ from updateMassBins_ function, select peaks from the same mass in the original input spectrum
  void FLASHDeconvAlgorithm::getCandidatePeakGroups_(const Matrix<int>& per_mass_abs_charge_ranges)
  {
    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];
    double tol = tolerance_[ms_level_ - 1];
    int charge_range = current_max_charge_;
    Size mass_bin_size = mass_bins_.size();
    int log_mz_peak_size = (int)log_mz_peaks_.size();
    // this stores which peak is now being considered per charge. Per charge, peak is considered from left (lowest m/z) to right (highest m/z).
    auto current_peak_index = std::vector<int>(charge_range, 0);
    deconvolved_spectrum_.reserve(mass_bins_.count());
    Size mass_bin_index = mass_bins_.find_first();
    auto peak_bin_numbers = std::vector<Size>(log_mz_peak_size);
    const double max_mass_dalton_tolerance = .16; // this is additional mass tolerance in Da to get more high signal-to-ratio peaks in this candidate peakgroup finding
    // per peak, store the m/z bin number for fast processing
    for (int i = 0; i < log_mz_peak_size; i++)
    {
      peak_bin_numbers[i] = getBinNumber_(log_mz_peaks_[i].logMz, mz_bin_min_value_, bin_mul_factor);
    }
    std::vector<double> total_harmonic_intensity(harmonic_charges_.size(), .0);
    std::vector<int> h_prev_iso(harmonic_charges_.size(), 0);
    std::vector<float> h_max_isotope_intensity(harmonic_charges_.size(), .0f);

    // main iteration. per_mass_abs_charge_ranges gives the range of charges for each mass bin
    while (mass_bin_index != mass_bins_.npos)
    {
      double log_m = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_mul_factor);
      double mass = exp(log_m);

      PeakGroup pg(1, per_mass_abs_charge_ranges(1, mass_bin_index) + 1, // make an empty peakGroup (mass)
                   is_positive_);

      pg.reserve(charge_range * 12);
      pg.setIsotopeDaDistance(iso_da_distance_);
      // the range of isotope span. For a given peak the peaks within the span are searched.
      Size right_index = avg_.getRightCountFromApex(mass);
      Size left_index = avg_.getLeftCountFromApex(mass);

      double total_signal_intensity = 0;
      std::fill(total_harmonic_intensity.begin(), total_harmonic_intensity.end(), .0);

      // scan through charge - from mass to m/z
      for (size_t j = per_mass_abs_charge_ranges(0, mass_bin_index); j <= (size_t)per_mass_abs_charge_ranges(1, mass_bin_index); j++)
      {
        int max_peak_index = -1;
        size_t abs_charge = j + 1;
        int bin_offset = bin_offsets_[j];

        if (mass_bin_index < (size_t)bin_offset)
        {
          continue;
        }

        Size b_index = mass_bin_index - bin_offset; // m/z bin
        int& cpi = current_peak_index[j];           // in this charge which peak is to be considered?
        double max_intensity = -1;

        while (cpi < log_mz_peak_size - 1) // scan through peaks from cpi
        {
          if (peak_bin_numbers[cpi] == b_index) // if the peak of consideration matches to this mass with charge abs_charge
          {
            double intensity = log_mz_peaks_[cpi].intensity;
            if (intensity > max_intensity) // compare with other matching peaks and select the most intense peak (in max_peak_index)
            {
              max_intensity = intensity;
              max_peak_index = cpi;
            }
          }
          else if (peak_bin_numbers[cpi] > b_index)
          {
            break;
          }
          cpi++;
        }

        if (max_peak_index < 0)
        {
          continue;
        }

        // Search for local max.
        if (max_peak_index > 0 && max_peak_index <= log_mz_peak_size && peak_bin_numbers[max_peak_index - 1] == b_index - 1 && log_mz_peaks_[max_peak_index - 1].intensity > max_intensity)
        {
          continue;
        }
        if (max_peak_index >= 0 && max_peak_index < log_mz_peak_size - 1 && peak_bin_numbers[max_peak_index + 1] == b_index + 1 && log_mz_peaks_[max_peak_index + 1].intensity > max_intensity)
        {
          continue;
        }

        // now we have a matching peak for this mass of charge  abs_charge. From here, isotope peaks are collected
        const double mz = log_mz_peaks_[max_peak_index].mz; // charged mz
        const double iso_delta = iso_da_distance_ / (double)abs_charge;
        double mz_delta = std::min(max_mass_dalton_tolerance / (double)abs_charge, 2.0 * tol * mz);

        double max_mz = mz;
        float max_peak_intensity = log_mz_peaks_[max_peak_index].intensity;
        float max_isotope_intensity = .0;
        int prev_iso = -1000;
        std::fill(h_prev_iso.begin(), h_prev_iso.end(), 0);
        std::fill(h_max_isotope_intensity.begin(), h_max_isotope_intensity.end(), .0);

        for (int peak_index = max_peak_index; peak_index < log_mz_peak_size; peak_index++)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          const float intensity = log_mz_peaks_[peak_index].intensity;
          double mz_diff = observed_mz - mz;
          int tmp_i = (int)round(mz_diff / iso_delta);

          if (observed_mz - max_mz > (double)right_index * iso_delta + mz_delta)
          {
            break;
          }

          if (abs(mz_diff - tmp_i * iso_delta) < mz_delta) // if peak is signal
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size && !(bin < previously_deconved_mass_bins_for_dummy_.size() && previously_deconved_mass_bins_for_dummy_[bin]))
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = (int)abs_charge;
              p.isotopeIndex = tmp_i;
              pg.push_back(p);
              if (max_peak_intensity < intensity)
              {
                max_peak_intensity = intensity;
              }
              if (prev_iso != tmp_i)
              {
                total_signal_intensity += max_isotope_intensity;
                max_isotope_intensity = .0;
              }
              max_isotope_intensity = std::max(max_isotope_intensity, intensity);
              prev_iso = tmp_i;
            }
          }
          else
          {
            for (Size l = 0; l < harmonic_charges_.size(); l++)
            {
              int hc = harmonic_charges_[l];
              if ((int)(hc * abs_charge) > current_max_charge_)
              {
                break;
              }
              double hiso_delta = iso_delta / hc;
              int tmp_hi = (int)round(mz_diff / hiso_delta);

              if ((double)tmp_hi / hc < tmp_i + max_mass_dalton_tolerance)
              {
                continue;
              }

              if ((double)tmp_hi / hc >= tmp_i + 1 - max_mass_dalton_tolerance)
              {
                break;
              }

              double err = abs(mz_diff - tmp_hi * hiso_delta);
              if (err < mz_delta)
              {
                if (h_prev_iso[l] != tmp_hi / hc)
                {
                  total_harmonic_intensity[l] += std::min(max_peak_intensity, h_max_isotope_intensity[l]);
                  h_max_isotope_intensity[l] = .0;
                }
                h_max_isotope_intensity[l] = std::max(h_max_isotope_intensity[l], intensity);
                h_prev_iso[l] = tmp_hi / hc;
              }
            }
          }
        }

        total_signal_intensity += max_isotope_intensity;
        for (Size l = 0; l < harmonic_charges_.size(); l++)
        {
          total_harmonic_intensity[l] += h_max_isotope_intensity[l];
        }

        max_isotope_intensity = .0;
        prev_iso = -1000;
        std::fill(h_prev_iso.begin(), h_prev_iso.end(), 0);
        std::fill(h_max_isotope_intensity.begin(), h_max_isotope_intensity.end(), .0);

        for (int peak_index = max_peak_index - 1; peak_index >= 0; peak_index--)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          const float intensity = log_mz_peaks_[peak_index].intensity;
          double mz_diff = observed_mz - mz;
          int tmp_i = (int)round(mz_diff / iso_delta);

          if (max_mz - observed_mz > (float)left_index * iso_delta + mz_delta)
          {
            break;
          }
          if (abs(mz_diff - tmp_i * iso_delta) < mz_delta)
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size && !(bin < previously_deconved_mass_bins_for_dummy_.size() && previously_deconved_mass_bins_for_dummy_[bin]))
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = (int)abs_charge;
              p.isotopeIndex = tmp_i;
              pg.push_back(p);
              // total_signal_intensity += intensity;
              if (max_peak_intensity < intensity)
              {
                max_peak_intensity = intensity;
              }
              if (prev_iso != tmp_i)
              {
                total_signal_intensity += max_isotope_intensity;
                max_isotope_intensity = .0;
              }
              max_isotope_intensity = std::max(max_isotope_intensity, intensity);
              prev_iso = tmp_i;
            }
          }
          else
          {
            for (Size l = 0; l < harmonic_charges_.size(); l++)
            {
              int hc = harmonic_charges_[l];
              if ((int)(hc * abs_charge) > current_max_charge_)
              {
                break;
              }
              double hiso_delta = iso_delta / hc;
              int tmp_hi = (int)round(mz_diff / hiso_delta);

              if ((double)tmp_hi / hc > tmp_i - max_mass_dalton_tolerance)
              {
                continue;
              }

              if ((double)tmp_hi / hc <= tmp_i - 1 + max_mass_dalton_tolerance)
              {
                break;
              }

              double err = abs(mz_diff - tmp_hi * hiso_delta);

              if (err < mz_delta)
              {
                if (h_prev_iso[l] != tmp_hi / hc)
                {
                  total_harmonic_intensity[l] += h_max_isotope_intensity[l];
                  h_max_isotope_intensity[l] = .0;
                }
                h_max_isotope_intensity[l] = std::max(h_max_isotope_intensity[l], intensity);
                h_prev_iso[l] = tmp_hi / hc;
              }
            }
          }
        }

        total_signal_intensity += max_isotope_intensity;
        for (Size l = 0; l < harmonic_charges_.size(); l++)
        {
          total_harmonic_intensity[l] += std::min(max_peak_intensity, h_max_isotope_intensity[l]);
        }
      }

      if (total_signal_intensity > *std::max_element(total_harmonic_intensity.begin(), total_harmonic_intensity.end())) //
      {
        double max_intensity = -1.0;
        double t_mass = .0;
        auto new_peaks = std::vector<LogMzPeak>();
        new_peaks.reserve(pg.size());
        for (auto& p : pg)
        {
          if (max_intensity < p.intensity)
          {
            max_intensity = p.intensity;
            t_mass = p.getUnchargedMass();
          }
        }
        double iso_tolerance = tol * t_mass;
        int min_off = 10000;
        int max_off = -1;
        int max_charge = -1;

        int apex_index = (int)avg_.getApexIndex(t_mass);
        for (auto& p : pg)
        {
          p.isotopeIndex = (int)round((p.getUnchargedMass() - t_mass) / iso_da_distance_);

          if (abs(t_mass - p.getUnchargedMass() + iso_da_distance_ * p.isotopeIndex) > iso_tolerance)
          {
            continue;
          }
          p.isotopeIndex += apex_index;
          new_peaks.push_back(p);
          min_off = min_off > p.isotopeIndex ? p.isotopeIndex : min_off;
          max_off = max_off < p.isotopeIndex ? p.isotopeIndex : max_off;
          max_charge = max_charge < p.abs_charge ? p.abs_charge : max_charge;
        }
        if (min_off != max_off)
        {
          pg.swap(new_peaks);
          pg.updateMonoMassAndIsotopeIntensities();

          if (pg.getMonoMass() < current_min_mass_ || pg.getMonoMass() > current_max_mass_)
          {
            continue;
          }
          pg.setScanNumber(deconvolved_spectrum_.getScanNumber());
          deconvolved_spectrum_.push_back(pg); //
        }
      }
      mass_bin_index = mass_bins_.find_next(mass_bin_index);
    }
  }

  DeconvolvedSpectrum& FLASHDeconvAlgorithm::getDeconvolvedSpectrum()
  {
    return deconvolved_spectrum_;
  }

  void FLASHDeconvAlgorithm::setTargetDummyType(PeakGroup::TargetDummyType target_dummy_type, DeconvolvedSpectrum& target_dspec_for_dummy_calcualtion)
  {
    target_dummy_type_ = target_dummy_type;
    target_dspec_for_dummy_calcualtion_ = &target_dspec_for_dummy_calcualtion;
  }

  // spectral deconvolution main function
  void FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_()
  {
    deconvolved_spectrum_.clear();
    int current_charge_range = current_max_charge_;
    int tmp_peak_cntr = current_charge_range - min_support_peak_count_;

    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;
    double mass_bin_max_value = std::min(log_mz_peaks_[log_mz_peaks_.size() - 1].logMz - filter_[tmp_peak_cntr], log(current_max_mass_ + (double)avg_.getRightCountFromApex(current_max_mass_) + 1.0));

    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];

    mass_bin_min_value_ = log(std::max(1.0, current_min_mass_ - avg_.getAverageMassDelta(current_min_mass_)));
    mz_bin_min_value_ = log_mz_peaks_[0].logMz;

    double mz_bin_max_value = log_mz_peaks_[log_mz_peaks_.size() - 1].logMz;
    Size mass_bin_number = getBinNumber_(mass_bin_max_value, mass_bin_min_value_, bin_mul_factor) + 1;
    bin_offsets_.clear();

    for (int i = 0; i < current_charge_range; i++)
    {
      bin_offsets_.push_back((int)round((mz_bin_min_value_ - filter_[i] - mass_bin_min_value_) * bin_mul_factor));
    }

    harmonic_bin_offset_matrix_.getEigenMatrix().resize(harmonic_charges_.size(), current_charge_range);
    
    
    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      for (int i = 0; i < current_charge_range; i++)
      {
        harmonic_bin_offset_matrix_(k, i) = (int)round((mz_bin_min_value_ - harmonic_filter_matrix_(k, i) - mass_bin_min_value_) * bin_mul_factor);
      }
    }

    Size mz_bin_number = getBinNumber_(mz_bin_max_value, mz_bin_min_value_, bin_mul_factor) + 1;
    auto mz_bin_intensities = std::vector<float>(mz_bin_number, .0f);
    updateMzBins_(mz_bin_number, mz_bin_intensities);
    mass_bins_ = boost::dynamic_bitset<>(mass_bin_number);

    if (target_dummy_type_ == PeakGroup::charge_dummy && !previously_deconved_mono_masses_for_dummy_.empty())
    {
      std::sort(previously_deconved_mono_masses_for_dummy_.begin(), previously_deconved_mono_masses_for_dummy_.end());
      previously_deconved_mass_bins_for_dummy_ = boost::dynamic_bitset<>(mass_bins_.size());
      // always positive
      unsigned bin_offset = (unsigned)round(tol_div_factor);
      for (double m : previously_deconved_mono_masses_for_dummy_)
      {
        if (m <= 0)
        {
          continue;
        }
        Size j = getBinNumber_(log(m), mass_bin_min_value_, bin_mul_factors_[ms_level_ - 1]);
        if (j >= bin_offset && j < previously_deconved_mass_bins_for_dummy_.size() - bin_offset - 1)
        {
          for (int k = -bin_offset; k <= (int) bin_offset; k++)
            previously_deconved_mass_bins_for_dummy_[j + k] = true;
        }
      }
    }

    if (!target_mono_masses_.empty())
    {
      target_mass_bins_.reset();
      target_mass_bins_ = boost::dynamic_bitset<>(mass_bins_.size());
      for (double& tm : target_mono_masses_)
      {
        for (int off = -1; off < 2; off++)
        {
          double m = tm + off * iso_da_distance_;
          double mass_delta = avg_.getMostAbundantMassDelta(m);

          Size j = getBinNumber_(log(m + mass_delta), mass_bin_min_value_, bin_mul_factors_[ms_level_ - 1]);
          if (j < 1)
          {
            continue;
          }

          if (j >= target_mass_bins_.size() - 2)
          {
            break;
          }

          target_mass_bins_[j - 1] = true;
          target_mass_bins_[j] = true;
          target_mass_bins_[j + 1] = true;
        }
      }
    }

    if (target_dummy_type_ != PeakGroup::isotope_dummy)
    {
      auto per_mass_abs_charge_ranges = updateMassBins_(mz_bin_intensities);
      getCandidatePeakGroups_(per_mass_abs_charge_ranges);
    }
    else
    {
      deconvolved_spectrum_ = *target_dspec_for_dummy_calcualtion_;
    }
    scoreAndFilterPeakGroups_();
  }

  void FLASHDeconvAlgorithm::scoreAndFilterPeakGroups_()
  {
    std::vector<PeakGroup> filtered_peak_groups;
    filtered_peak_groups.reserve(deconvolved_spectrum_.size());

    double tol = tolerance_[ms_level_ - 1];
#pragma omp parallel default(none) shared(tol, filtered_peak_groups)
    {
      std::vector<PeakGroup> filtered_peak_groups_private;
      filtered_peak_groups_private.reserve(deconvolved_spectrum_.size());
#pragma omp for nowait schedule(static)
      for (int i = 0; i < (int)deconvolved_spectrum_.size(); i++)
      {
        int offset = 0;
        auto peak_group = deconvolved_spectrum_[i];
        peak_group.setTargetDummyType(target_dummy_type_);
        float prev_cos = peak_group.getIsotopeCosine();
        float cos = getIsotopeCosineAndDetermineIsotopeIndex(peak_group.getMonoMass(), peak_group.getIsotopeIntensities(), offset, avg_, -peak_group.getMinNegativeIsotopeIndex(), -1, allowed_iso_error_,
                                                             target_dummy_type_);
        auto prev_mono_mass = peak_group.getMonoMass() + offset * iso_da_distance_;

        peak_group.setIsotopeCosine(cos);

        // first filtration to remove false positives before further processing.
        if (cos < std::min(.5, min_isotope_cosine_[ms_level_ - 1]) - .1)
        {
          continue;
        }

        for (int k = 0; k < 10; k++)
        {
          auto noisy_peaks = peak_group.recruitAllPeaksInSpectrum(deconvolved_spectrum_.getOriginalSpectrum(), tol, avg_, peak_group.getMonoMass() + offset * iso_da_distance_, excluded_peak_mzs_);
          // min cosine is checked in here. mono mass is also updated one last time. SNR, per charge SNR, and avg errors are updated here.
          offset = peak_group.updateQscore(noisy_peaks, avg_, min_isotope_cosine_[ms_level_ - 1], allowed_iso_error_);
          if (offset == 0)
          {
            break;
          }
        }

        if (peak_group.empty() || peak_group.getMonoMass() < current_min_mass_ || peak_group.getMonoMass() > current_max_mass_)
        {
          continue;
        }

        if (std::abs(prev_mono_mass - peak_group.getMonoMass()) > 3) // if they are off by more than 3, they are different envelopes.
        {
          continue;
        }

        auto [z1, z2] = peak_group.getAbsChargeRange();

        if (z1 > low_charge_ && (z2 - z1) < min_support_peak_count_)
        {
          continue;
        }

        if (target_dummy_type_ == PeakGroup::TargetDummyType::isotope_dummy)
        {
          if (peak_group.getIsotopeCosine() < prev_cos * .98) // if target cosine and isotope dummy cosine are too different, we do not take this dummy.
            continue;
        }

        if (!target_mono_masses_.empty())
        {
          double delta = peak_group.getMonoMass() * tolerance_[ms_level_ - 1] * 2;
          auto upper = std::upper_bound(target_mono_masses_.begin(), target_mono_masses_.end(), peak_group.getMonoMass() + delta);

          while (!peak_group.isTargeted())
          {
            if (upper != target_mono_masses_.end())
            {
              if (std::abs(*upper - peak_group.getMonoMass()) < delta)
              {
                peak_group.setTargeted();
              }
              if (peak_group.getMonoMass() - *upper > delta)
              {
                break;
              }
            }
            if (upper == target_mono_masses_.begin())
            {
              break;
            }
            --upper;
          }
        }

        float snr_threshold = .5;
        if (!peak_group.isTargeted() && (peak_group.getQscore() <= 0 || peak_group.getSNR() < snr_threshold)) // snr check prevents harmonics or noise.
        {
          continue;
        }

        if (target_dummy_type_ == PeakGroup::charge_dummy && !previously_deconved_mono_masses_for_dummy_.empty())
        {
          bool exclude = false;
          double delta = peak_group.getMonoMass() * tolerance_[ms_level_ - 1];
          auto upper = std::upper_bound(previously_deconved_mono_masses_for_dummy_.begin(), previously_deconved_mono_masses_for_dummy_.end(), peak_group.getMonoMass() + delta);

          while (!exclude)
          {
            if (upper != previously_deconved_mono_masses_for_dummy_.end())
            {
              if (std::abs(*upper - peak_group.getMonoMass()) < delta)
              {
                exclude = true;
              }
              if (peak_group.getMonoMass() - *upper > delta)
              {
                break;
              }
            }
            if (upper == previously_deconved_mono_masses_for_dummy_.begin())
            {
              break;
            }
            --upper;
          }
          if (exclude)
          {
            continue;
          }
        }

        if (peak_group.getQscore() <= 0)
        {
          continue;
        }
        filtered_peak_groups_private.push_back(peak_group);
      }

#ifdef _OPENMP
  #pragma omp for schedule(static) ordered
      for (int i = 0; i < omp_get_num_threads(); i++)
      {
  #pragma omp ordered
        filtered_peak_groups.insert(filtered_peak_groups.end(), filtered_peak_groups_private.begin(), filtered_peak_groups_private.end());
      }
#else
      filtered_peak_groups = filtered_peak_groups_private;
#endif
    }

    deconvolved_spectrum_.setPeakGroups(filtered_peak_groups);
    deconvolved_spectrum_.sort();

    removeChargeErrorPeakGroups_(deconvolved_spectrum_);
    removeOverlappingPeakGroups_(deconvolved_spectrum_, tol * tol_div_factor * 1.5);
  }

  float FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(const double mono_mass, const std::vector<float>& per_isotope_intensities, int& offset, const PrecalculatedAveragine& avg,
                                                                       int iso_int_shift, int window_width, int allowed_iso_error_for_second_best_cos, PeakGroup::TargetDummyType target_dummy_type)
  {
    offset = 0;
    if ((int)per_isotope_intensities.size() < min_iso_size_ + iso_int_shift)
    {
      return .0;
    }
    auto iso = avg.get(mono_mass);

    int iso_size = (int)iso.size();
    int right = (int)avg.getApexIndex(mono_mass) / 4 + 1;
    int left = right;

    if (window_width >= 0)
    {
      right = std::min(right, window_width);
      left = std::min(left, window_width);
    }

    float max_cos = -1000;
    float second_max_cos = -1000;
    int second_max_offset = -1000;
    int max_isotope_index = (int)per_isotope_intensities.size(); // exclusive
    int min_isotope_index = -1;                                  // inclusive

    left -= iso_int_shift;
    right += iso_int_shift;

    for (int i = 0; i < max_isotope_index; i++)
    {
      if (per_isotope_intensities[i] <= 0)
      {
        continue;
      }

      if (min_isotope_index < 0)
      {
        min_isotope_index = i;
      }
    }
    if (max_isotope_index - min_isotope_index < min_iso_size_)
    {
      return .0;
    }

    for (int tmp_offset = -left; tmp_offset <= right; tmp_offset++)
    {
      float tmp_cos = getCosine(per_isotope_intensities, min_isotope_index, max_isotope_index, iso,
                                iso_size, // apex_index,
                                tmp_offset, min_iso_size_);

      if (max_cos < tmp_cos)
      {
        max_cos = tmp_cos;
        offset = tmp_offset;
      }
    }

    if (target_dummy_type == PeakGroup::TargetDummyType::isotope_dummy)
    {
      for (int tmp_offset = offset - 3; tmp_offset <= offset + 3; tmp_offset++)
      {
        if (abs(offset - tmp_offset) <= allowed_iso_error_for_second_best_cos) //
        {
          continue;
        }
        if (tmp_offset < -left || tmp_offset > right)
        {
          continue;
        }
        float tmp_cos = getCosine(per_isotope_intensities, min_isotope_index, max_isotope_index, iso,
                                  iso_size, // apex_index,
                                  tmp_offset, min_iso_size_);

        if (second_max_cos < tmp_cos && tmp_cos < max_cos)
        {
          second_max_cos = tmp_cos;
          second_max_offset = tmp_offset;
        }
      }
      max_cos = second_max_cos;
      offset = second_max_offset;
    }

    offset -= iso_int_shift;
    return max_cos;
  }

  float FLASHDeconvAlgorithm::getCosine(const std::vector<float>& a, int a_start, int a_end, const IsotopeDistribution& b, int b_size, int offset, int min_iso_size)
  {
    float n = .0, a_norm = .0;
    a_start = std::max(0, a_start);
    a_end = std::min((int)a.size(), a_end);

    if (a_end - a_start < min_iso_size)
    {
      return 0;
    }

    int max_intensity_index = 0;
    float max_intensity = 0;

    for (int j = a_start; j < a_end; j++)
    {
      int i = j - offset;
      a_norm += a[j] * a[j];

      if (max_intensity < a[j])
      {
        max_intensity = a[j];
        max_intensity_index = j;
      }

      if (i < 0 && a[j] > 0)
      {
        // n -= a[j] * b[0].getIntensity();
      }
      else if (i >= b_size || i < 0 || b[i].getIntensity() <= 0)
      {
        continue;
      }
      else
      {
        n += a[j] * b[i].getIntensity(); //
      }
    }

    // two consecutive isotopes around the max intensity isotope
    if (min_iso_size > 0)
    {
      if (max_intensity_index == a_end - 1)
      {
        if (max_intensity_index > 0 && a[max_intensity_index - 1] == 0)
        {
          return 0;
        }
      }
      else if (max_intensity_index == a_start)
      {
        if (max_intensity_index + 1 < (int)a.size() && a[max_intensity_index + 1] == 0)
        {
          return 0;
        }
      }
      else if (max_intensity_index > 0 && max_intensity_index + 1 < (int)a.size())
      {
        if (a[max_intensity_index + 1] == 0 && a[max_intensity_index - 1] == 0)
        {
          return 0;
        }
      }
    }

    if (a_norm <= 0)
    {
      return 0;
    }

    return n / sqrt(a_norm);
  }

  void FLASHDeconvAlgorithm::removeChargeErrorPeakGroups_(DeconvolvedSpectrum& dspec) const
  {
    std::map<double, std::set<int>> peak_to_pgs;
    std::map<double, float> mz_to_intensities;
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(dspec.size());
    std::vector<float> overlap_intensity(dspec.size(), .0f);

    for (Size i = 0; i < dspec.size(); i++)
    {
      auto& pg = dspec[i];
      for (auto& p : pg)
      {
        peak_to_pgs[p.mz].insert((int)i);
        mz_to_intensities[p.mz] = p.intensity;
      }
    }

    for (auto& e : peak_to_pgs)
    {
      auto& pg_is = e.second;
      double pmz = e.first;
      float pint = mz_to_intensities[pmz];

      if (pg_is.size() == 1)
      {
        continue;
      }

      for (auto i : pg_is)
      {
        bool is_overlap = false;
        double mass1 = dspec[i].getMonoMass();
        int repz1 = (int)round(mass1 / (pmz - FLASHDeconvHelperStructs::getChargeMass(is_positive_)));
        for (auto j : pg_is)
        {
          if (i == j)
          {
            continue;
          }
          double mass2 = dspec[j].getMonoMass();
          int repz2 = (int)round(mass2 / (pmz - FLASHDeconvHelperStructs::getChargeMass(is_positive_)));
          if (repz1 == repz2)
          {
            continue;
          }
          if (dspec[i].getChargeSNR(repz1) > dspec[j].getChargeSNR(repz2) * 2.0) // if ith is way better than jth, jth is overlapped not ith
          {
            continue;
          }
          is_overlap = true;
          break;
        }
        if (is_overlap)
          overlap_intensity[i] += pint;
      }
    }

    for (Size i = 0; i < dspec.size(); i++)
    {
      if (dspec[i].getTargetDummyType() != target_dummy_type_)
      {
        continue;
      }
      if (!dspec[i].isTargeted() &&
          overlap_intensity[i] >= dspec[i].getIntensity() * .5) // If the overlapped intensity takes more than 50% total intensity then it is a peakgroup with a charge error. the smaller, the harsher
      {
        continue;
      }
      if ((ms_level_ == 1 && dspec[i].getRepAbsCharge() < min_abs_charge_) || dspec[i].getRepAbsCharge() > max_abs_charge_)
      {
        continue;
      }
      filtered_pg_vec.push_back(dspec[i]);
    }
    dspec.setPeakGroups(filtered_pg_vec);
  }

  void FLASHDeconvAlgorithm::removeOverlappingPeakGroups_(DeconvolvedSpectrum& dspec, double tol)
  {
    if (dspec.empty())
    {
      return;
    }
    std::vector<PeakGroup> filtered_pg_vec; //
    filtered_pg_vec.reserve(dspec.size());
    double start_mass = dspec[0].getMonoMass();
    float local_max_SNR = 0;
    Size local_max_index = 0;

    for (Size i = 0; i < dspec.size(); i++)
    {
      double mass = dspec[i].getMonoMass();
      if (mass - start_mass > mass * tol)
      {
        if (!dspec[local_max_index].isTargeted()) // targeted ones were already push_backed.
        {
          if (dspec[local_max_index].getTargetDummyType() == target_dummy_type_)
            filtered_pg_vec.push_back(dspec[local_max_index]);
        }
        start_mass = mass;
        local_max_SNR = 0;
      }

      if (local_max_SNR < dspec[i].getSNR())
      {
        local_max_SNR = dspec[i].getSNR();
        local_max_index = i;
      }
      if (dspec[i].isTargeted())
      {
        if (dspec[i].getTargetDummyType() == target_dummy_type_)
          filtered_pg_vec.push_back(dspec[i]);
      }
    }

    if (local_max_SNR > 0)
    {
      if (!dspec[local_max_index].isTargeted()) // targeted ones were already push_backed.
      {
        if (dspec[local_max_index].getTargetDummyType() == target_dummy_type_)
          filtered_pg_vec.push_back(dspec[local_max_index]);
      }
    }

    dspec.setPeakGroups(filtered_pg_vec);
    std::vector<PeakGroup>().swap(filtered_pg_vec);
  }

  void FLASHDeconvAlgorithm::setTargetMasses(const std::vector<double>& masses, bool excluded)
  {
    if (excluded)
    {
      excluded_masses_.clear();
      excluded_masses_.reserve(masses.size() * 30);
    }
    else
    {
      target_mono_masses_.clear();
      target_mono_masses_.reserve(masses.size() * 3);
    }
    for (auto& m : masses)
    {
      int start = 0;
      int end = 0;
      if (excluded)
      {
        end = (int)(avg_.getApexIndex(m) + avg_.getRightCountFromApex(m));
      }
      for (int j = start; j <= end + 1; j++)
      {
        if (excluded)
          excluded_masses_.push_back(m + iso_da_distance_ * j);
        else
          target_mono_masses_.push_back(m + iso_da_distance_ * j);
      }
    }
  }

  bool FLASHDeconvAlgorithm::registerPrecursor_(const std::vector<DeconvolvedSpectrum>& survey_scans, const std::map<int, std::vector<std::vector<float>>>& precursor_map_for_real_time_acquisition)
  {
    if (ms_level_ == 1 || (survey_scans.empty() && precursor_map_for_real_time_acquisition.empty()))
    {
      return false;
    }
    deconvolved_spectrum_.setPrecursorIntensity(.0);
    double start_mz = 0;
    double end_mz = 0;
    for (auto& precursor : deconvolved_spectrum_.getOriginalSpectrum().getPrecursors())
    {
      for (auto& activation_method : precursor.getActivationMethods())
      {
        deconvolved_spectrum_.setActivationMethod(activation_method);
        if (deconvolved_spectrum_.getActivationMethod() == Precursor::HCID)
        {
          deconvolved_spectrum_.setActivationMethod(Precursor::HCD);
        }
        break;
      }
      deconvolved_spectrum_.setPrecursor(precursor);

      double loffset = precursor.getIsolationWindowLowerOffset();
      double uoffest = precursor.getIsolationWindowUpperOffset();
      loffset = loffset <= 0 ? isolation_window_size_ / 2.0 : loffset;
      uoffest = uoffest <= 0 ? isolation_window_size_ / 2.0 : uoffest;

      start_mz = loffset > 100.0 ? loffset : -loffset + precursor.getMZ();
      end_mz = uoffest > 100.0 ? uoffest : uoffest + precursor.getMZ();
    }
    if (!precursor_map_for_real_time_acquisition.empty() && deconvolved_spectrum_.getPrecursorPeakGroup().empty())
    {
      auto map = precursor_map_for_real_time_acquisition.lower_bound(deconvolved_spectrum_.getScanNumber());
      while (map != precursor_map_for_real_time_acquisition.begin())
      {
        --map;
        if (map->first < deconvolved_spectrum_.getScanNumber() - 50)
        {
          return false;
        }

        if (map != precursor_map_for_real_time_acquisition.end())
        {
          for (auto& smap : map->second)
          {
            if (abs(start_mz - smap[3]) < .001 && abs(end_mz - smap[4]) < .001)
            {
              LogMzPeak precursor_log_mz_peak(deconvolved_spectrum_.getPrecursor(), is_positive_);
              precursor_log_mz_peak.abs_charge = (int)smap[1];
              precursor_log_mz_peak.isotopeIndex = 0;
              precursor_log_mz_peak.mass = smap[0];
              precursor_log_mz_peak.intensity = smap[6];

              Precursor precursor(deconvolved_spectrum_.getPrecursor());
              precursor.setCharge(precursor_log_mz_peak.abs_charge);
              precursor.setIntensity(smap[5]);
              deconvolved_spectrum_.setPrecursor(precursor);

              PeakGroup precursor_pg(deconvolved_spectrum_.getPrecursorPeakGroup());
              precursor_pg.push_back(precursor_log_mz_peak);
              precursor_pg.setAbsChargeRange((int)smap[7], (int)smap[8]);
              precursor_pg.setChargeIsotopeCosine(precursor_log_mz_peak.abs_charge, smap[9]);
              precursor_pg.setChargeSNR(precursor_log_mz_peak.abs_charge, smap[10]); // cnsr
              precursor_pg.setIsotopeCosine(smap[11]);
              precursor_pg.setSNR(smap[12]);
              precursor_pg.setChargeScore(smap[13]);
              precursor_pg.setAvgPPMError(smap[14]);
              precursor_pg.Qscore(smap[2]);
              precursor_pg.setRepAbsCharge((int)smap[1]);
              precursor_pg.updateMonoMassAndIsotopeIntensities();
              precursor_pg.setScanNumber(map->first);
              deconvolved_spectrum_.setPrecursorPeakGroup(precursor_pg);

              for (int i = (int)survey_scans.size() - 1; i >= 0; i--)
              {
                auto precursor_spectrum = survey_scans[i];
                if (precursor_spectrum.getScanNumber() != deconvolved_spectrum_.getPrecursorScanNumber())
                {
                  continue;
                }

                auto o_spec = precursor_spectrum.getOriginalSpectrum();
                auto spec_iterator = o_spec.PosBegin(start_mz);
                while (spec_iterator != o_spec.end() && spec_iterator->getMZ() < end_mz)
                {
                  spec_iterator++;
                }

                double max_score = -1.0;
                for (auto& pg : precursor_spectrum)
                {
                  if (pg[0].mz > end_mz || pg[pg.size() - 1].mz < start_mz)
                  {
                    continue;
                  }
                  if (abs(pg.getMonoMass() - smap[0]) > 5.0)
                  {
                    continue;
                  }
                  double max_intensity = .0;
                  const LogMzPeak* tmp_precursor = nullptr;
                  for (auto& tmp_peak : pg)
                  {
                    if (tmp_peak.mz < start_mz)
                    {
                      continue;
                    }
                    if (tmp_peak.mz > end_mz)
                    {
                      break;
                    }

                    if (tmp_peak.intensity < max_intensity)
                    {
                      continue;
                    }
                    max_intensity = tmp_peak.intensity;
                    tmp_precursor = &tmp_peak;
                  }

                  if (tmp_precursor == nullptr || tmp_precursor->abs_charge < min_abs_charge_ || tmp_precursor->abs_charge > max_abs_charge_)
                  {
                    continue;
                  }

                  auto score = pg.getChargeSNR(tmp_precursor->abs_charge);
                  // most intense one should determine the mass
                  if (score < max_score)
                  {
                    continue;
                  }

                  max_score = score;
                  deconvolved_spectrum_.setPrecursorPeakGroup(pg);
                }
              }
              return true;
            }
          }
        }
      }
      return false;
    }

    double max_score = -1.0;

    for (int i = (int)survey_scans.size() - 1; i >= 0; i--)
    {
      auto precursor_spectrum = survey_scans[i];
      if (deconvolved_spectrum_.getPrecursorScanNumber() == 0)
      {
        deconvolved_spectrum_.setPrecursorScanNumber(precursor_spectrum.getScanNumber());
      }
      if (precursor_spectrum.empty())
      {
        continue;
      }
      auto o_spec = precursor_spectrum.getOriginalSpectrum();

      for (auto& pg : precursor_spectrum)
      {
        if (pg[0].mz > end_mz || pg[pg.size() - 1].mz < start_mz)
        {
          continue;
        }

        double max_intensity = .0;
        const LogMzPeak* tmp_precursor = nullptr;

        int c = int(.5 + pg.getMonoMass() / start_mz);
        for (auto& tmp_peak : pg)
        {
          if (tmp_peak.abs_charge != c)
          {
            continue;
          }

          if (tmp_peak.mz < start_mz || tmp_peak.mz > end_mz)
          {
            continue;
          }

          if (tmp_peak.intensity < max_intensity)
          {
            continue;
          }
          max_intensity = tmp_peak.intensity;
          tmp_precursor = &tmp_peak;
        }

        if (tmp_precursor == nullptr)
        {
          continue;
        }

        auto score = pg.getChargeSNR(tmp_precursor->abs_charge); // most intense one should determine the mass

        if (score < max_score)
        {
          continue;
        }

        Precursor precursor_pg(deconvolved_spectrum_.getPrecursor());
        precursor_pg.setCharge(tmp_precursor->is_positive ? tmp_precursor->abs_charge : -tmp_precursor->abs_charge);

        deconvolved_spectrum_.setPrecursor(precursor_pg);
        max_score = score;

        deconvolved_spectrum_.setPrecursorPeakGroup(pg);
        deconvolved_spectrum_.setPrecursorScanNumber(precursor_spectrum.getScanNumber());
      }
      if (!deconvolved_spectrum_.getPrecursorPeakGroup().empty())
      {
        break;
      }
    }
    return deconvolved_spectrum_.getPrecursorPeakGroup().empty();
  }

  void FLASHDeconvAlgorithm::setAveragine(const FLASHDeconvAlgorithm::PrecalculatedAveragine& avg)
  {
    avg_ = avg;
  }
  double FLASHDeconvAlgorithm::getMassFromMassBin_(Size mass_bin, double bin_mul_factor) const
  {
    return exp(getBinValue_(mass_bin, mass_bin_min_value_, bin_mul_factor));
  }

  double FLASHDeconvAlgorithm::getMzFromMzBin_(Size mass_bin, double bin_mul_factor) const
  {
    return exp(getBinValue_(mass_bin, mz_bin_min_value_, bin_mul_factor));
  }

} // namespace OpenMS