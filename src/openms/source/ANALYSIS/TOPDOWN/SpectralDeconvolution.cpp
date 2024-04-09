// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim, Jaekwan Kim $
// $Authors: Kyowon Jeong, Jihyung Kim, Jaekwan Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  /// harmonic charge factors that will be considered for harmonic mass reduction.
  inline const std::vector<int> harmonic_charges_ {2, 3, 5, 7, 11};
  /// this is additional mass tolerance in Da to get more high signal-to-ratio peaks in this candidate peakgroup finding
  inline const double max_mass_dalton_tolerance = .16;
  /// high and low charges are differently deconvolved. This value determines the (inclusive) threshold for low charge.
  inline const int low_charge_ = 10; // 10 inclusive
  inline const double tol_div_factor = 2.5;                           // use narrow tolerance for deconvolution and at the end use the input tolerance to filter out overlapping masses.

  SpectralDeconvolution::SpectralDeconvolution(): DefaultParamHandler("SpectralDeconvolution")
  {
    defaults_.setValue("tol", DoubleList {10.0, 10.0},
                       "ppm tolerance for MS1, 2, ... (e.g., -tol 10.0 5.0 to specify 10.0 and 5.0 ppm for MS1 and MS2, respectively)");

    defaults_.setValue("min_mass", 50.0, "Minimum mass (Da)");
    defaults_.setValue("max_mass", 100000.0, "Maximum mass (Da)");

    defaults_.setValue("min_charge", 1, "Minimum charge state for MS1 spectra (can be negative for negative mode)");
    defaults_.setValue("max_charge", 100, "Maximum charge state for MS1 spectra (can be negative for negative mode)");

    defaults_.setValue("precursor_charge", 0,
                       "Charge state of the target precursor. All precursor charge is fixed to this value. "
                       "This parameter is useful for targeted studies where MS2 spectra are generated from a fixed precursor (e.g., Native-MS). ");
    defaults_.setMinInt("precursor_charge", 0);
    // defaults_.addTag("precursor_charge", "advanced");

    defaults_.setValue(
      "precursor_mz", 0.0,
      "Target precursor m/z value. This option must be used with -target_precursor_charge option. Otherwise, it will be ignored. "
      "If -precursor_charge option is used but this option is not used, the precursor m/z value written in MS2 spectra will be used by default. ");
    defaults_.setMinFloat("precursor_mz", 0.0);
    defaults_.addTag("precursor_mz", "advanced");

    defaults_.setValue("min_cos", DoubleList {.85, .85},
                       "Cosine similarity thresholds between avg. and observed isotope pattern for MS1, 2, ... (e.g., -min_cos 0.3 0.6 to specify 0.3 "
                       "and 0.6 for MS1 and MS2, respectively)");
    defaults_.addTag("min_cos", "advanced");
    defaults_.setValue("min_snr", DoubleList {1.0, 1.0},
                       "Minimum charge SNR (the SNR of the isotope pattern of a specific charge) thresholds for MS1, 2, ... (e.g., -min_snr 1.0 0.6 to "
                       "specify 1.0 and 0.6 for MS1 and MS2, respectively)");
    defaults_.addTag("min_snr", "advanced");
    defaults_.setValue("max_qvalue", DoubleList {1.0, 1.0},
                       "Qvalue thresholds for MS1, 2, ... Effective only when FDR estimation is active. (e.g., -max_qvalue 0.1 0.2 to specify 0.1 and "
                       "0.2 for MS1 and MS2, respectively)");
    defaults_.addTag("max_qvalue", "advanced");
    defaults_.setValue("allowed_isotope_error", 1,
                       "Allowed isotope index error for decoy and qvalue report. If it is set to 2, for example, +-2 isotope errors are "
                       "not counted as false. Beta version.");
    defaults_.addTag("allowed_isotope_error", "advanced");

    defaultsToParam_();
  }

  // Calculate the nominal (integer) mass from double mass. Multiplying 0.999497 to the original mass and then rounding reduce the error between the
  // original and nominal masses.
  int SpectralDeconvolution::getNominalMass(const double mass)
  {
    return (int)round(mass * 0.999497);
  }

  void SpectralDeconvolution::setTargetPrecursorCharge_()
  {
    const auto& spec = deconvolved_spectrum_.getOriginalSpectrum();
    if (spec.getPrecursors().empty() && target_precursor_mz_ == 0)
    {
      OPENMS_LOG_INFO << "Attempted to set target precursor charge but failed - no precursor is found in MS2 spectra. Specify target precursor m/z "
                         "with -target_precursor_mz option"
                      << std::endl;
      return;
    }

    auto precursor = spec.getPrecursors()[0];
    double target_precursor_mass
      = (precursor.getMZ() - FLASHDeconvHelperStructs::getChargeMass(target_precursor_charge_ > 0)) * std::abs(target_precursor_charge_);
    precursor.setCharge(target_precursor_charge_);
    PeakGroup precursorPeakGroup(1, std::abs(target_precursor_charge_), target_precursor_charge_ > 0);
    precursorPeakGroup.push_back(FLASHDeconvHelperStructs::LogMzPeak());
    precursorPeakGroup.setMonoisotopicMass(target_precursor_mass);
    precursorPeakGroup.setSNR(1.0);

    precursorPeakGroup.setRepAbsCharge(std::abs(target_precursor_charge_));
    precursorPeakGroup.setChargeSNR(std::abs(target_precursor_charge_), 1.0);
    precursorPeakGroup.setQscore(1.0);

    deconvolved_spectrum_.setPrecursor(precursor);
    deconvolved_spectrum_.setPrecursorPeakGroup(precursorPeakGroup);
  }

  // The main function called from outside. precursor_map_for_FLASHIda is used to read FLASHIda information
  void SpectralDeconvolution::performSpectrumDeconvolution(const MSSpectrum& spec, const int scan_number, const PeakGroup& precursor_peak_group)
  {
    // First prepare for decoy runs.
    iso_da_distance_ = target_decoy_type_ == PeakGroup::noise_decoy
                         ? Constants::ISOTOPE_MASSDIFF_55K_U * sqrt(13.0) / 5
                         : Constants::ISOTOPE_MASSDIFF_55K_U; // sqrt(13.0)/5.0 Da is used instead of C13 - C12 to make sure masses detected with this
                                                              // nonsensical mass difference are not true.
    previously_deconved_peak_masses_for_decoy_.clear();
    previously_deconved_mass_bins_for_decoy_.reset();

    if (target_decoy_type_ == PeakGroup::charge_decoy) // charge decoy
    {
      for (const auto& pg : *target_dspec_for_decoy_calculation_)
      {
        for (const auto& p : pg)
        {
          for (int i = -allowed_iso_error_ - 6; i <= allowed_iso_error_ + 6; i++)
            previously_deconved_peak_masses_for_decoy_.push_back(i * iso_da_distance_ + p.getUnchargedMass());
        }
      }
    }

    ms_level_ = spec.getMSLevel();
    deconvolved_spectrum_ = DeconvolvedSpectrum(scan_number);
    deconvolved_spectrum_.setOriginalSpectrum(spec);

    // here register targeted peak mzs etc.
    // for MSn (n>1) register precursor peak and peak group.
    if (ms_level_ > 1)
    {
      for (const auto& precursor : deconvolved_spectrum_.getOriginalSpectrum().getPrecursors())
      {
        for (const auto& activation_method : precursor.getActivationMethods())
        {
          deconvolved_spectrum_.setActivationMethod(activation_method);
          if (deconvolved_spectrum_.getActivationMethod() == Precursor::HCID) { deconvolved_spectrum_.setActivationMethod(Precursor::HCD); }
          break;
        }
        deconvolved_spectrum_.setPrecursor(precursor);
      }

      if (target_precursor_charge_ != 0 || target_precursor_mz_ > 0) { setTargetPrecursorCharge_(); }

      if (deconvolved_spectrum_.getPrecursorPeakGroup().empty())
      {
        if (! precursor_peak_group.empty())
        {
          deconvolved_spectrum_.setPrecursorPeakGroup(precursor_peak_group);
          deconvolved_spectrum_.setPrecursorScanNumber(precursor_peak_group.getScanNumber());
          Precursor precursor(deconvolved_spectrum_.getPrecursor());
          int abs_charge = (int)round(precursor_peak_group.getMonoMass() / precursor.getMZ());
          precursor.setCharge(precursor_peak_group.isPositive() ? abs_charge : -abs_charge);
          deconvolved_spectrum_.setPrecursor(precursor);
        }
      }
    }

    // based on MS level, adjust charge and mass ranges. Precursor charge and mass determine those.
    current_max_charge_ = deconvolved_spectrum_.getCurrentMaxAbsCharge(max_abs_charge_); //
    current_max_mass_ = deconvolved_spectrum_.getCurrentMaxMass(max_mass_);
    current_min_mass_ = deconvolved_spectrum_.getCurrentMinMass(min_mass_);

    // set universal pattern filter and harmonic pattern filters
    setFilters_();
    // LogMzPeaks are generated from raw peaks

    updateLogMzPeaks_();
    if (log_mz_peaks_.empty()) { return; }

    // This is the main FLASHDeconv function in which deconvolution is performed.
    generatePeakGroupsFromSpectrum_();
  }

  void SpectralDeconvolution::updateMembers_()
  {
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

    current_max_mass_ = max_mass_ = param_.getValue("max_mass");
    current_min_mass_ = min_mass_ = param_.getValue("min_mass");

    bin_mul_factors_.clear();
    tolerance_ = param_.getValue("tol");

    for (double& j : tolerance_)
    {
      j *= 1e-6;
      // j /= tol_div_factor; // finder bins are far better.
      bin_mul_factors_.push_back(1.0 / j * tol_div_factor);
    }

    min_isotope_cosine_ = param_.getValue("min_cos");
    min_snr_ = param_.getValue("min_snr");
    max_qvalue_ = param_.getValue("max_qvalue");

    allowed_iso_error_ = param_.getValue("allowed_isotope_error");

    target_precursor_mz_ = param_.getValue("precursor_mz");
    target_precursor_charge_ = param_.getValue("precursor_charge");
  }

  const FLASHDeconvHelperStructs::PrecalculatedAveragine& SpectralDeconvolution::getAveragine()
  {
    return avg_;
  }

  void SpectralDeconvolution::calculateAveragine(const bool use_RNA_averagine)
  {
    CoarseIsotopePatternGenerator generator(300);

    auto iso = use_RNA_averagine ? generator.estimateFromRNAWeight(current_max_mass_) : generator.estimateFromPeptideWeight(current_max_mass_);
    iso.trimRight(0.0001 * iso.getMostAbundant().getIntensity());
    auto max_isotope = std::max(200, (int)iso.size());

    generator.setMaxIsotope(max_isotope);
    avg_ = FLASHDeconvHelperStructs::PrecalculatedAveragine(50, current_max_mass_, 25, generator, use_RNA_averagine);
    avg_.setMaxIsotopeIndex((int)(max_isotope - 1));
  }

  // generate filters
  void SpectralDeconvolution::setFilters_()
  {
    filter_.clear();
    int charge_range = current_max_charge_;
    for (int i = 0; i < charge_range; i++)
    {
      filter_.push_back(-log(i + 1)); //+
    }

    harmonic_filter_matrix_ = Matrix<double>(harmonic_charges_.size(), charge_range, .0);

    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      int hc = harmonic_charges_[k];
      int n = hc / 2;
      // int m = (1 + hc) / 2;
      for (int i = 0; i < charge_range; i++)
      {
        double a = i > 0 ? exp(-filter_[i - 1]) : 0;
        double b = exp(-filter_[i]);
        harmonic_filter_matrix_.setValue(k, i, -log(b - (b - a) * n / hc));
        // harmonic_filter_matrix_.setValue(k * 2 + 1, i, -log(b - (b - a) * m / hc));
      }
    }
  }

  // Generate uncharged log mz transformated peaks
  void SpectralDeconvolution::updateLogMzPeaks_()
  {
    log_mz_peaks_.clear();
    log_mz_peaks_.reserve(deconvolved_spectrum_.getOriginalSpectrum().size());

    for (const auto& peak : deconvolved_spectrum_.getOriginalSpectrum())
    {
      if (peak.getIntensity() <= 0) //
      {
        continue;
      }

      LogMzPeak log_mz_peak(peak, is_positive_);
      log_mz_peaks_.push_back(log_mz_peak);
    }
  }

  // from bin to raw value
  double SpectralDeconvolution::getBinValue_(const Size bin, const double min_value, const double bin_mul_factor)
  {
    return min_value + (double)bin / bin_mul_factor;
  }

  // from value to bin number
  Size SpectralDeconvolution::getBinNumber_(const double value, const double min_value, const double bin_mul_factor)
  {
    if (value < min_value) { return 0; }
    return (Size)round((value - min_value) * bin_mul_factor);
  }

  // From log mz to mz bins.
  void SpectralDeconvolution::updateMzBins_(const Size bin_number, std::vector<float>& mz_bin_intensities)
  {
    mz_bins_ = boost::dynamic_bitset<>(bin_number);
    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];

    for (const auto& p : log_mz_peaks_)
    {
      Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_mul_factor);
      if (bi >= bin_number) { break; }
      mz_bins_.set(bi);

      mz_bin_intensities[bi] += p.intensity;
    }
  }

  // Find candidate mass bins from the current spectrum. The runtime of FLASHDeconv is determined by this function.
  void SpectralDeconvolution::updateCandidateMassBins_(std::vector<float>& mass_intensities, const std::vector<float>& mz_intensities)
  { //
    Size mz_bin_index = mz_bins_.find_first();
    auto mz_bin_index_reverse = std::vector<Size>();
    mz_bin_index_reverse.reserve(mz_bins_.count());
    // invert mz bins so charges are counted from small to large given a mass

    while (mz_bin_index != mz_bins_.npos)
    {
      mz_bin_index_reverse.push_back(mz_bin_index);
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }
    size_t h_charge_size = harmonic_charges_.size();
    long bin_end = (long)mass_bins_.size();

    auto support_peak_count = std::vector<unsigned short>(mass_bins_.size(), 0); // per mass bin how many peaks are present

    // to calculate continuous charges, the previous charge value per mass should be stored
    auto prev_charges = std::vector<unsigned short>(mass_bins_.size(), current_max_charge_ + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prev_intensities = std::vector<float>(mass_bins_.size(), .0f);

    mass_intensities = std::vector<float>(mass_bins_.size(), .0f);

    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];
    std::vector<float> sub_max_h_intensity(h_charge_size, .0f);

    for (auto iter = mz_bin_index_reverse.rbegin(); iter < mz_bin_index_reverse.rend(); iter++)
    {
      mz_bin_index = *iter;
      const float intensity = mz_intensities[mz_bin_index];
      const double log_mz = getBinValue_(mz_bin_index, mz_bin_min_value_, bin_mul_factor); // uncharged log mz;
      const double mz = exp(log_mz);                                                       // uncharged mz
      const double iso_div_by_mz = iso_da_distance_ / mz;
      // scan through charges
      for (int j = 0; j < current_max_charge_; j++) // take all charge one ?
      {
        // mass is given by shifting by bin_offsets_[j]
        const long mass_bin_index = (long)mz_bin_index + bin_offsets_[j];

        if (mass_bin_index < 0) { continue; }
        if (mass_bin_index >= bin_end) { break; }

        if (! previously_deconved_mass_bins_for_decoy_.empty() && previously_deconved_mass_bins_for_decoy_[mass_bin_index]) { continue; }

        auto& spc = support_peak_count[mass_bin_index];
        const int abs_charge = (j + 1);
        float& prev_intensity = prev_intensities[mass_bin_index];
        auto& prev_charge = prev_charges[mass_bin_index];
        const bool charge_not_continous = prev_charge - j != -1 && (prev_charge <= current_max_charge_);
        bool pass_first_check = false;

        // intensity ratio between consecutive charges should not exceed the factor.
        const float highest_factor = 10.0f;
        const float factor = abs_charge <= low_charge_ ? highest_factor : (highest_factor / 2 + highest_factor / 2 * low_charge_ / (float)abs_charge);
        // intensity ratio between consecutive charges for possible harmonic should be within this factor

        const float hfactor = factor / 2.0f;
        // intensity of previous charge
        // intensity ratio between current and previous charges
        float intensity_ratio = prev_intensity <= 0 ? (factor + 1) : (intensity / prev_intensity);
        intensity_ratio = intensity_ratio < 1 ? 1.0f / intensity_ratio : intensity_ratio;
        float support_peak_intensity = 0;
        // check if peaks of continuous charges are present
        std::fill(sub_max_h_intensity.begin(), sub_max_h_intensity.end(), .0f);

        // if charge not continuous or intensity ratio is too high reset support_peak_count
        if (charge_not_continous || intensity_ratio > factor) { spc = 0; }
        else
        {
          pass_first_check = true;
          if (spc == 0 && abs_charge > low_charge_) { support_peak_intensity = prev_intensity; }
        }

        // for low charges, check isotope peak presence.
        if (! pass_first_check && abs_charge <= low_charge_)
        { // for low charges
          for (int d = 1; d >= -1; d -= 2)
          {
            bool iso_exist = false;
            int next_iso_bin = 0;
            const int nib = (int)getBinNumber_(log_mz + d * iso_div_by_mz / abs_charge, mz_bin_min_value_, bin_mul_factor);
            const int nibr
              = abs_charge > 1 ? (int)getBinNumber_(log_mz + (d * iso_div_by_mz / (abs_charge - 1)), mz_bin_min_value_, bin_mul_factor) : 0;
            const int nibl = (int)getBinNumber_(log_mz + (d * iso_div_by_mz / (abs_charge + 1)), mz_bin_min_value_, bin_mul_factor);
            if (abs(nib - nibr) < tol_div_factor
                || abs(nib - nibl)
                     < tol_div_factor) // if different charges are not distinguishable, we ignore. Not informative and the source of the errors.
              break;

            for (int t = -1; t < 2; t++)
            {
              int nibt = nib + t;
              if (std::abs(nibt - (int)mz_bin_index) >= tol_div_factor && nibt > 0 && nibt < (int)mz_bins_.size() && mz_bins_[nibt])
              {
                iso_exist = true;
                pass_first_check = true;

                if (next_iso_bin == 0 || mz_intensities[next_iso_bin] < mz_intensities[nibt]) { next_iso_bin = nibt; }
              }
            }

            // harmonic check
            if (iso_exist)
            {
              const double h_threshold = intensity + mz_intensities[next_iso_bin]; // std::min(intensity, mz_intensities[next_iso_bin]); //

              for (size_t k = 0; k < h_charge_size; k++)
              {
                const int hc = harmonic_charges_[k];
                int harmonic_cntr = 0;
                if (ms_level_ > 1 && hc * abs_charge > current_max_charge_) { break; }

                const int hdiff = (int)round((double)(next_iso_bin - mz_bin_index)) / hc * (hc / 2);
                const int next_harmonic_iso_bin = (int)mz_bin_index + hdiff; //(int)getBinNumber_(log_mz + hdiff, mz_bin_min_value_, bin_mul_factor);
                // check if there are harmonic peaks between the current peak and the next isotope peak.

                // no perfect filtration. Just obvious ones are filtered out by checking if a peak is in the harmonic position and the intensity ratio
                // is within two folds from the current peak (specified by mz_bin_index)
                if (std::abs(next_harmonic_iso_bin - (int)mz_bin_index) >= tol_div_factor && next_harmonic_iso_bin >= 0
                    && next_harmonic_iso_bin < (int)mz_bins_.size() && mz_bins_[next_harmonic_iso_bin]
                    && mz_intensities[next_harmonic_iso_bin] > h_threshold / 2 && mz_intensities[next_harmonic_iso_bin] < h_threshold * 2)
                {
                  harmonic_cntr++;
                  sub_max_h_intensity[k] += mz_intensities[next_harmonic_iso_bin];
                  // sub_max_h_intensity[k] = std::max(sub_max_h_intensity[k] , mz_intensities[next_harmonic_iso_bin]);
                }

                if (harmonic_cntr > 0) { pass_first_check = false; }
              }
            }
            if (pass_first_check)
            {
              support_peak_intensity += mz_intensities[next_iso_bin];
              // support_peak_intensity = std::max(support_peak_intensity, mz_intensities[next_iso_bin]);
            }
          }
          pass_first_check &= *std::max_element(sub_max_h_intensity.begin(), sub_max_h_intensity.end()) <= 0;
        }

        if (pass_first_check)
        {
          if (prev_charge - j == -1) // check harmonic artifacts for high charge ranges
          {
            float max_intensity = intensity;
            float min_intensity = prev_intensity;
            if (prev_intensity <= .0)
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

            const float high_threshold = max_intensity * hfactor;
            const float low_threshold = min_intensity / hfactor;

            bool is_harmonic = false;

            // check if harmonic peaks are present with different harmonic multiple factors (2, 3, 5, 7, 11  defined in harmonic_charges_).
            int min_dis = tol_div_factor + 1;
            for (size_t k = 0; k < h_charge_size; k++)
            {
              if (ms_level_ > 1 && harmonic_charges_[k] * abs_charge > current_max_charge_) break;
              float harmonic_intensity = 0;
              for (int t = -tol_div_factor; t <= tol_div_factor; t++)
              {
                long hmz_bin_index = mass_bin_index - harmonic_bin_offset_matrix_.getValue(k, j) + t;
                if (hmz_bin_index > 0 && hmz_bin_index != (long)mz_bin_index && hmz_bin_index < (int)mz_bins_.size() && mz_bins_[hmz_bin_index])
                {
                  float h_intensity = mz_intensities[hmz_bin_index];
                  if (h_intensity > low_threshold && h_intensity < high_threshold)
                  {
                    is_harmonic = true;
                    if (abs(t) < min_dis)
                    {
                      harmonic_intensity = h_intensity;
                      min_dis = abs(t);
                    }
                  }
                }
              }
              sub_max_h_intensity[k] += harmonic_intensity;
            }

            if (! is_harmonic) // if it is not harmonic
            {
              mass_intensities[mass_bin_index] += intensity + support_peak_intensity;

              if (! mass_bins_[mass_bin_index])
              {
                spc++;
                if (spc >= min_support_peak_count_ || spc >= abs_charge / 2) {
                   mass_bins_[mass_bin_index] = true;
                }
              }
            }
            else // if harmonic
            {
              mass_intensities[mass_bin_index]
                -= *std::max_element(sub_max_h_intensity.begin(), sub_max_h_intensity.end()); // std::min(max_h_intensity, intensity);
              if (spc > 0) { spc--; }
            }
          }
          else if (abs_charge <= low_charge_) // for low charge, include the mass if isotope is present
          {
            mass_intensities[mass_bin_index] += intensity + support_peak_intensity;
            if (! mass_bins_[mass_bin_index])
            {
              spc++;
              mass_bins_[mass_bin_index] = true;
            }
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
  Matrix<int> SpectralDeconvolution::filterMassBins_(const std::vector<float>& mass_intensities)
  {
    Matrix<int> abs_charge_ranges(2, mass_bins_.size(), INT_MAX);
    for (Size i = 0; i < mass_bins_.size(); i++)
    {
      abs_charge_ranges.setValue(1, (int)i, INT_MIN);
    }
    Size mz_bin_index = mz_bins_.find_first();
    long bin_size = (long)mass_bins_.size();

    auto to_skip = mass_bins_.flip();
    mass_bins_.reset();

    const int select_top_N = 2; // select top N charges per peak. We allow up to 2 just to consider frequent coelution.
    std::vector<long> max_indices(select_top_N, -1);
    std::vector<int> max_intensity_abs_charge_ranges(select_top_N, -1);

    while (mz_bin_index != mz_bins_.npos)
    {
      std::fill(max_indices.begin(), max_indices.end(), -1);
      std::fill(max_intensity_abs_charge_ranges.begin(), max_intensity_abs_charge_ranges.end(), -1);

      float max_intensity = 0;

      for (int j = 0; j < current_max_charge_; j++)
      {
        long mass_bin_index = (long)mz_bin_index + bin_offsets_[j];

        if (mass_bin_index < 0) { continue; }
        if (mass_bin_index >= bin_size) { break; }

        if (! target_mono_masses_.empty() && target_mass_bins_[mass_bin_index])
        {
          float t = mass_intensities[mass_bin_index];

          if (t == 0)
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
          if (to_skip[mass_bin_index]) { continue; }

          float t = mass_intensities[mass_bin_index];
          if (t == 0) // no signal
          {
            continue;
          }

          if (max_intensity == 0 || max_intensity < t)
          {
            // store best values after shift by 1.
            bool big_diff = t > max_intensity * 2.0;
            for (int i = select_top_N - 1; i > 0; i--)
            {
              max_indices[i] = big_diff ? -1 : max_indices[i - 1]; // if big change, only take the best charge
              max_intensity_abs_charge_ranges[i] = max_intensity_abs_charge_ranges[i - 1];
            }
            max_indices[0] = mass_bin_index;
            max_intensity_abs_charge_ranges[0] = j;
            max_intensity = t;
          }
        }
      }
      for (int i = 0; i < select_top_N; i++)
      {
        long max_index = max_indices[i];
        int max_intensity_abs_charge_range = max_intensity_abs_charge_ranges[i];
        if (max_index >= 0 && max_index < bin_size)
        {
          abs_charge_ranges.setValue(0, max_index, std::min(abs_charge_ranges.getValue(0, max_index), max_intensity_abs_charge_range));
          abs_charge_ranges.setValue(1, max_index, std::max(abs_charge_ranges.getValue(1, max_index), max_intensity_abs_charge_range));
          mass_bins_[max_index] = true;
        }
      }
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }

    return abs_charge_ranges;
  }

  // update mass bins which will be used to select peaks in the input spectrum...
  Matrix<int> SpectralDeconvolution::updateMassBins_(const std::vector<float>& mz_intensities)
  {
    std::vector<float> mass_intensities;
    updateCandidateMassBins_(mass_intensities, mz_intensities);

    auto per_mass_abs_charge_ranges = filterMassBins_(mass_intensities);

    return per_mass_abs_charge_ranges;
  }

  // With mass_bins_ from updateMassBins_ function, select peaks from the same mass in the original input spectrum
  void SpectralDeconvolution::getCandidatePeakGroups_(const Matrix<int>& per_mass_abs_charge_ranges)
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
    // per peak, store the m/z bin number for fast processing
    for (int i = 0; i < log_mz_peak_size; i++)
    {
      peak_bin_numbers[i] = getBinNumber_(log_mz_peaks_[i].logMz, mz_bin_min_value_, bin_mul_factor);
    }

    // main iteration. per_mass_abs_charge_ranges gives the range of charges for each mass bin
    while (mass_bin_index != mass_bins_.npos)
    {
      double log_m = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_mul_factor);
      double mass = exp(log_m);

      PeakGroup pg(1, per_mass_abs_charge_ranges.getValue(1, mass_bin_index) + 1, // make an empty peakGroup (mass)
                   is_positive_);

      pg.reserve(charge_range * 12);
      pg.setTargetDecoyType(target_decoy_type_);
      pg.setIsotopeDaDistance(iso_da_distance_);
      // the range of isotope span. For a given peak the peaks within the span are searched.
      Size right_index = avg_.getRightCountFromApex(mass);
      Size left_index = avg_.getLeftCountFromApex(mass);

      // scan through charge - from mass to m/z
      for (size_t j = per_mass_abs_charge_ranges.getValue(0, mass_bin_index); j <= (size_t)per_mass_abs_charge_ranges.getValue(1, mass_bin_index); j++)
      {
        int max_peak_index = -1;
        size_t abs_charge = j + 1;
        int bin_offset = bin_offsets_[j];

        if (mass_bin_index < (size_t)bin_offset) { continue; }

        Size b_index = mass_bin_index - bin_offset; // m/z bin
        int& cpi = current_peak_index[j];           // in this charge which peak is to be considered?
        double max_intensity = -1;

        while (cpi < log_mz_peak_size - 1)      // scan through peaks from cpi
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
          else if (peak_bin_numbers[cpi] > b_index) { break; }
          cpi++;
        }

        if (max_peak_index < 0) { continue; }

        // Search for local max.
        if (max_peak_index > 0 && max_peak_index <= log_mz_peak_size && peak_bin_numbers[max_peak_index - 1] == b_index - 1
            && log_mz_peaks_[max_peak_index - 1].intensity > max_intensity)
        {
          continue;
        }
        if (max_peak_index >= 0 && max_peak_index < log_mz_peak_size - 1 && peak_bin_numbers[max_peak_index + 1] == b_index + 1
            && log_mz_peaks_[max_peak_index + 1].intensity > max_intensity)
        {
          continue;
        }

        // now we have a matching peak for this mass of charge  abs_charge. From here, isotope peaks are collected
        const double mz = log_mz_peaks_[max_peak_index].mz;                                   // charged mz
        const double iso_delta = iso_da_distance_ / (double)abs_charge;
        double mz_delta = std::min(max_mass_dalton_tolerance / (double)abs_charge, tol * mz); //

        double max_mz = mz;

        for (int peak_index = max_peak_index; peak_index < log_mz_peak_size; peak_index++)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          double mz_diff = observed_mz - mz;
          int tmp_i = (int)round(mz_diff / iso_delta);

          if (observed_mz - max_mz > (double)right_index * iso_delta + mz_delta) { break; }

          if (abs(mz_diff - tmp_i * iso_delta) < mz_delta) // if peak is signal
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size && ! (bin < previously_deconved_mass_bins_for_decoy_.size() && previously_deconved_mass_bins_for_decoy_[bin]))
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = (int)abs_charge;
              p.isotopeIndex = tmp_i;
              pg.push_back(p);
            }
          }
        }

        for (int peak_index = max_peak_index - 1; peak_index >= 0; peak_index--)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          double mz_diff = observed_mz - mz;
          int tmp_i = (int)round(mz_diff / iso_delta);

          if (max_mz - observed_mz > (float)left_index * iso_delta + mz_delta) { break; }
          if (abs(mz_diff - tmp_i * iso_delta) < mz_delta)
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size && ! (bin < previously_deconved_mass_bins_for_decoy_.size() && previously_deconved_mass_bins_for_decoy_[bin]))
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = (int)abs_charge;
              p.isotopeIndex = tmp_i;
              pg.push_back(p);
            }
          }
        }
      }

      if (! pg.empty()) // total_signal_intensity > 0)// *std::max_element(total_harmonic_intensity.begin(), total_harmonic_intensity.end())) //
      {
        double max_intensity = -1.0;
        double t_mass = .0;
        auto new_peaks = std::vector<LogMzPeak>();
        new_peaks.reserve(pg.size());
        for (const auto& p : pg)
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

          if (abs(t_mass - p.getUnchargedMass() + iso_da_distance_ * p.isotopeIndex) > iso_tolerance) { continue; }
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

          if (pg.getMonoMass() < current_min_mass_ || pg.getMonoMass() > current_max_mass_) { continue; }
          pg.setScanNumber(deconvolved_spectrum_.getScanNumber());
          pg.sort();

          deconvolved_spectrum_.push_back(pg); //
        }
      }
      mass_bin_index = mass_bins_.find_next(mass_bin_index);
    }
  }

  DeconvolvedSpectrum& SpectralDeconvolution::getDeconvolvedSpectrum()
  {
    return deconvolved_spectrum_;
  }

  void SpectralDeconvolution::setTargetDecoyType(PeakGroup::TargetDecoyType target_decoy_type,
                                                 const DeconvolvedSpectrum& target_dspec_for_decoy_calcualtion)
  {
    target_decoy_type_ = target_decoy_type;
    target_dspec_for_decoy_calculation_ = &target_dspec_for_decoy_calcualtion;
  }

  // spectral deconvolution main function
  void SpectralDeconvolution::generatePeakGroupsFromSpectrum_()
  {
    deconvolved_spectrum_.clear();
    int current_charge_range = current_max_charge_;
    int tmp_peak_cntr = current_charge_range - min_support_peak_count_;

    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;
    double mass_bin_max_value = std::min(log_mz_peaks_.back().logMz - filter_[tmp_peak_cntr],
                                         log(current_max_mass_ + (double)avg_.getRightCountFromApex(current_max_mass_) + 1.0));

    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];

    // mass_bin_min_value_ = log(std::max(1.0, current_min_mass_ - avg_.getAverageMassDelta(current_min_mass_)));
    mass_bin_min_value_ = log(std::max(1.0, 50 - avg_.getAverageMassDelta(50)));
    mz_bin_min_value_ = log_mz_peaks_[0].logMz;

    double mz_bin_max_value = log_mz_peaks_.back().logMz;
    Size mass_bin_number = getBinNumber_(mass_bin_max_value, mass_bin_min_value_, bin_mul_factor) + 1;
    bin_offsets_.clear();

    for (int i = 0; i < current_charge_range; i++)
    {
      bin_offsets_.push_back((int)round((mz_bin_min_value_ - filter_[i] - mass_bin_min_value_) * bin_mul_factor));
    }

    harmonic_bin_offset_matrix_ = Matrix<int>(harmonic_filter_matrix_.rows(), harmonic_filter_matrix_.cols(), 0);

    for (Size k = 0; k < harmonic_filter_matrix_.rows(); k++)
    {
      for (Size i = 0; i < harmonic_filter_matrix_.cols(); i++)
      {
        harmonic_bin_offset_matrix_.setValue(
          k, i, (int)round((mz_bin_min_value_ - harmonic_filter_matrix_.getValue(k, i) - mass_bin_min_value_) * bin_mul_factor));
      }
    }

    Size mz_bin_number = getBinNumber_(mz_bin_max_value, mz_bin_min_value_, bin_mul_factor) + 1;
    auto mz_bin_intensities = std::vector<float>(mz_bin_number, .0f);
    updateMzBins_(mz_bin_number, mz_bin_intensities);
    mass_bins_ = boost::dynamic_bitset<>(mass_bin_number);

    if (! previously_deconved_peak_masses_for_decoy_.empty())
    {
      std::sort(previously_deconved_peak_masses_for_decoy_.begin(), previously_deconved_peak_masses_for_decoy_.end());
      previously_deconved_mass_bins_for_decoy_ = boost::dynamic_bitset<>(mass_bins_.size());
      // always positive
      int bin_offset = (int)round(tol_div_factor) + 1;
      for (double m : previously_deconved_peak_masses_for_decoy_)
      {
        if (m <= 0) { continue; }
        Size j = getBinNumber_(log(m), mass_bin_min_value_, bin_mul_factors_[ms_level_ - 1]);
        if ((int)j >= bin_offset && j < previously_deconved_mass_bins_for_decoy_.size() - bin_offset)
        {
          for (int k = -bin_offset; k <= bin_offset; k++)
            previously_deconved_mass_bins_for_decoy_[j + k] = true;
        }
      }
    }

    if (! target_mono_masses_.empty())
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
          if (j < 1) { continue; }

          if (j >= target_mass_bins_.size() - 2) { break; }

          target_mass_bins_[j - 1] = true;
          target_mass_bins_[j] = true;
          target_mass_bins_[j + 1] = true;
        }
      }
    }

    if (target_decoy_type_ != PeakGroup::isotope_decoy)
    {
      auto per_mass_abs_charge_ranges = updateMassBins_(mz_bin_intensities);
      getCandidatePeakGroups_(per_mass_abs_charge_ranges);
    }
    else { deconvolved_spectrum_ = *target_dspec_for_decoy_calculation_; }
    scoreAndFilterPeakGroups_();
  }

  void SpectralDeconvolution::scoreAndFilterPeakGroups_()
  {
    double tol = tolerance_[ms_level_ - 1];
    auto selected = boost::dynamic_bitset<>(deconvolved_spectrum_.size());

  #pragma omp parallel for default(none) shared(tol, selected)
    for (int i = 0; i < (int)deconvolved_spectrum_.size(); i++)
    {
      int offset = 0;
      auto& peak_group = deconvolved_spectrum_[i];

      auto prev_peak_group(peak_group);
      peak_group.setTargetDecoyType(target_decoy_type_);
      bool is_isotope_decoy = target_decoy_type_ == PeakGroup::TargetDecoyType::isotope_decoy;
      int window_width = is_isotope_decoy ? 0 : -1;

      float cos = getIsotopeCosineAndIsoOffset(peak_group.getMonoMass(), peak_group.getIsotopeIntensities(), offset, avg_,
                                               -peak_group.getMinNegativeIsotopeIndex(), window_width, allowed_iso_error_, target_decoy_type_);
      peak_group.setIsotopeCosine(cos);
      // first filtration to remove false positives before further processing.

      if (cos < std::min(min_isotope_cosine_[ms_level_ - 1], 0.5)) { continue; }

      int num_iteration = 10;
      bool mass_determined = false;
      for (int k = 0; k < num_iteration; k++)
      {
        auto noisy_peaks = peak_group.recruitAllPeaksInSpectrum(deconvolved_spectrum_.getOriginalSpectrum(), tol, avg_,
                                                                peak_group.getMonoMass() + offset * iso_da_distance_);
        // min cosine is checked in here. mono mass is also updated one last time. SNR, per charge SNR, and avg errors are updated here.
        const auto& [z1, z2] = peak_group.getAbsChargeRange();

        offset = peak_group.updateQscore(noisy_peaks, deconvolved_spectrum_.getOriginalSpectrum(), avg_, min_isotope_cosine_[ms_level_ - 1], tol,
                                         (z1 + z2) < 2 * low_charge_, allowed_iso_error_, false);

        if (offset == 0 || is_isotope_decoy)
        {
          peak_group.updateQscore(noisy_peaks, deconvolved_spectrum_.getOriginalSpectrum(), avg_, min_isotope_cosine_[ms_level_ - 1], tol,
                                  (z1 + z2) < 2 * low_charge_, allowed_iso_error_, true);
          mass_determined = true;
          break;
        }
      }
      if (! mass_determined || peak_group.empty() || peak_group.getQscore() <= 0 || peak_group.getMonoMass() < current_min_mass_
          || peak_group.getMonoMass() > current_max_mass_)
        continue;

      auto [z1, z2] = peak_group.getAbsChargeRange();

      if (z1 > low_charge_ && (z2 - z1) < min_support_peak_count_) { continue; }

      if (is_isotope_decoy)
      {
        if (prev_peak_group.getIsotopeCosine() * (1 - .0026) > peak_group.getIsotopeCosine()) // a magic number to generate isotope decoy masses.
          continue;
        peak_group.setQscore(prev_peak_group.getQscore());
      }

      if (! target_mono_masses_.empty())
      {
        double delta = peak_group.getMonoMass() * tolerance_[ms_level_ - 1] * 2;
        auto upper = std::upper_bound(target_mono_masses_.begin(), target_mono_masses_.end(), peak_group.getMonoMass() + delta);

        while (! peak_group.isTargeted())
        {
          if (upper != target_mono_masses_.end())
          {
            if (std::abs(*upper - peak_group.getMonoMass()) < delta) { peak_group.setTargeted(); }
            if (peak_group.getMonoMass() - *upper > delta) { break; }
          }
          if (upper == target_mono_masses_.begin()) { break; }
          --upper;
        }
      }

      double snr_threshold = min_snr_[ms_level_ - 1];
      if (! peak_group.isTargeted()
          && (peak_group.getChargeSNR(peak_group.getRepAbsCharge()) < snr_threshold)) // snr check prevents harmonics or noise.
      {
        continue;
      }

  #pragma omp critical
      selected[i] = true;
    }

    Size selected_count = selected.count();
    if (selected_count == 0)
    {
      deconvolved_spectrum_.clear();
      return;
    }

    std::vector<PeakGroup> filtered_peak_groups;
    filtered_peak_groups.reserve(selected_count
                                 + (target_decoy_type_ != PeakGroup::TargetDecoyType::charge_decoy ? 0 : target_dspec_for_decoy_calculation_->size()));
    // filtered_peak_groups.reserve(selected_count);
    Size index = selected.find_first();
    while (index != selected.npos)
    {
      filtered_peak_groups.push_back(deconvolved_spectrum_[index]);
      index = selected.find_next(index);
    }

    if (target_decoy_type_ == PeakGroup::TargetDecoyType::charge_decoy)
    {
      filtered_peak_groups.insert(filtered_peak_groups.end(), target_dspec_for_decoy_calculation_->begin(), target_dspec_for_decoy_calculation_->end());
    }
    deconvolved_spectrum_.setPeakGroups(filtered_peak_groups);
    deconvolved_spectrum_.sort();

    removeOverlappingPeakGroups(deconvolved_spectrum_, 0, PeakGroup::TargetDecoyType::non_specific);
    removeChargeErrorPeakGroups_(deconvolved_spectrum_, target_decoy_type_);
    removeOverlappingPeakGroups(deconvolved_spectrum_, tol, target_decoy_type_);

    removeExcludedMasses_(deconvolved_spectrum_, previously_deconved_peak_masses_for_decoy_);
    removeExcludedMasses_(deconvolved_spectrum_, excluded_masses_);

    // final harmonics removal
    selected = boost::dynamic_bitset<>(deconvolved_spectrum_.size());
    filtered_peak_groups = std::vector<PeakGroup>();

  #pragma omp parallel for default(none) shared(filtered_peak_groups, tol, selected, harmonic_charges_)
    for (int i = 0; i < (int)deconvolved_spectrum_.size(); i++)
    {
      auto peak_group = deconvolved_spectrum_[i];
      bool pass = true;
      const auto& [z1, z2] = peak_group.getAbsChargeRange();

      for (auto hz : harmonic_charges_)
      {
        PeakGroup pg(std::min(current_max_charge_, z1 * hz), std::min(current_max_charge_, z2 * hz), is_positive_);
        pg.setTargetDecoyType(peak_group.getTargetDecoyType());
        pg.setMonoisotopicMass(peak_group.getMonoMass() * hz);
        auto nps = pg.recruitAllPeaksInSpectrum(deconvolved_spectrum_.getOriginalSpectrum(), tol, avg_, pg.getMonoMass());
        pg.updateQscore(nps, deconvolved_spectrum_.getOriginalSpectrum(), avg_, min_isotope_cosine_[ms_level_ - 1], tol,
                        (z1 + z2) * hz < 2 * low_charge_, allowed_iso_error_, true);

        if (pg.getQscore() > 0 && pg.getSNR() > peak_group.getSNR())
        {
          pass = false;
          break;
        }
      }
      if (! pass) continue;

  #pragma omp critical
      selected[i] = true;
    }

    filtered_peak_groups.reserve(selected.count());
    index = selected.find_first();

    while (index != selected.npos)
    {
      filtered_peak_groups.push_back(deconvolved_spectrum_[index]);
      index = selected.find_next(index);
    }

    deconvolved_spectrum_.setPeakGroups(filtered_peak_groups);
  }

  float SpectralDeconvolution::getIsotopeCosineAndIsoOffset(double mono_mass,
                                                            const std::vector<float>& per_isotope_intensities,
                                                            int& offset,
                                                            const PrecalculatedAveragine& avg,
                                                            int iso_int_shift,
                                                            int window_width,
                                                            int allowed_isotope_error,
                                                            PeakGroup::TargetDecoyType target_decoy_type)
  {
    offset = 0;
    if ((int)per_isotope_intensities.size() < min_iso_size + iso_int_shift) { return .0; }
    auto iso = avg.get(mono_mass);

    int right = (int)avg.getApexIndex(mono_mass) / 4 + 1;
    int left = right;
    if (target_decoy_type == PeakGroup::TargetDecoyType::isotope_decoy && allowed_isotope_error >= 0 && right < allowed_isotope_error) return 0;

    right += iso_int_shift;
    left -= iso_int_shift;
    float max_cos = -1000;
    int max_isotope_index = (int)per_isotope_intensities.size(); // exclusive
    int min_isotope_index = -1;                                  // inclusive

    for (int i = 0; i < max_isotope_index; i++)
    {
      if (per_isotope_intensities[i] <= 0) { continue; }

      if (min_isotope_index < 0) { min_isotope_index = i; }
    }
    if (max_isotope_index - min_isotope_index < min_iso_size) { return .0; }

    std::vector<std::pair<int, float>> offset_cos;
    offset_cos.reserve(right + left + 1);

    for (int tmp_offset = -left; tmp_offset <= right; tmp_offset++)
    {
      if (target_decoy_type != PeakGroup::TargetDecoyType::isotope_decoy && window_width >= 0 && abs(tmp_offset - iso_int_shift) > window_width)
        continue;
      float tmp_cos = getCosine(per_isotope_intensities, min_isotope_index, max_isotope_index, iso, tmp_offset, min_iso_size,
                                target_decoy_type == PeakGroup::TargetDecoyType::noise_decoy);
      offset_cos.emplace_back(tmp_offset, tmp_cos);
    }

    std::sort(offset_cos.begin(), offset_cos.end(),
              [](const std::pair<int, float>& p1, const std::pair<int, float>& p2) { return p1.second > p2.second; });

    for (const auto& [o, c] : offset_cos)
    {
      if (o > right || o < -left) continue;
      if (window_width >= 0 && abs(o - iso_int_shift) > window_width) //
        continue;
      offset = o;
      max_cos = c;
      break;
    }

    if (target_decoy_type == PeakGroup::TargetDecoyType::isotope_decoy && allowed_isotope_error >= 0)
    {
      int original_offset = offset;
      max_cos = -1000;

      for (const auto& [o, c] : offset_cos)
      {
        if (abs(original_offset - o) <= allowed_isotope_error) //
          continue;

        if (max_cos < 0)
        {
          offset = o;
          max_cos = c;
          break;
        }
      }
    }

    max_cos = std::max(max_cos, .0f);
    offset -= iso_int_shift;

    return max_cos;
  }

  float SpectralDeconvolution::getCosine(const std::vector<float>& a,
                                         int a_start,
                                         int a_end,
                                         const IsotopeDistribution& b,
                                         int offset,
                                         int min_iso_len,
                                         bool decoy)
  {
    float n = .0, a_norm = .0, b_norm = decoy ? .0 : 1.0f;
    a_start = std::max(0, a_start);
    a_end = std::min((int)a.size(), a_end);

    if (a_end - a_start < min_iso_len) { return 0; }

    float max_intensity = 0;

    if (decoy)
    {
      for (int i = 0; i < (int)b.size(); i++)
      {
        int ni = (i % 2 == 0 ? i + 1 : i - 1) % (int)b.size();
        b_norm += b[ni].getIntensity() * b[ni].getIntensity();
      }
    }

    for (int j = a_start; j < a_end; j++)
    {
      int i = j - offset;
      a_norm += a[j] * a[j];

      if (max_intensity < a[j])
      {
        max_intensity = a[j];
      }
      //
      if (i < 0 && a[j] > 0)
      {
        // n -= a[j] * b[0].getIntensity();
      }
      else if (i >= (int)b.size() || i < 0 || b[i].getIntensity() <= 0) { continue; }
      else
      {
        if (decoy)
        {
          int ni = (i % 2 == 0 ? i + 1 : i - 1) % (int)b.size();
          n += a[j] * b[ni].getIntensity(); //
        }
        else
          n += a[j] * b[i].getIntensity();
      }
    }

    if (a_norm <= 0) { return 0; }

    return n / sqrt(a_norm * b_norm);
  }


  void SpectralDeconvolution::removeChargeErrorPeakGroups_(DeconvolvedSpectrum& dspec, const PeakGroup::TargetDecoyType& target_decoy_type) const
  {
    std::map<double, std::set<int>> peak_to_pgs;
    std::map<double, float> mz_to_intensities;
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(dspec.size());
    std::vector<float> overlap_intensity(dspec.size(), .0f);

    for (Size i = 0; i < dspec.size(); i++)
    {
      const auto& pg = dspec[i];
      for (auto& p : pg)
      {
        peak_to_pgs[p.mz].insert((int)i);
        mz_to_intensities[p.mz] = p.intensity;
      }
    }

    for (const auto& e : peak_to_pgs)
    {
      const auto& pg_is = e.second;
      double pmz = e.first;
      float pint = mz_to_intensities[pmz];

      if (pg_is.size() == 1) { continue; }

      for (auto i : pg_is)
      {
        bool is_overlap = false;
        double mass1 = dspec[i].getMonoMass();
        int repz1 = (int)round(mass1 / (pmz - FLASHDeconvHelperStructs::getChargeMass(is_positive_)));
        const double high_snr = 10;
        if (dspec[i].getSNR() > high_snr && dspec[i].getChargeSNR(repz1) > high_snr) // if signal is almost perfect, skip the filtration
          continue;

        for (auto j : pg_is)
        {
          if (i == j) { continue; }
          double mass2 = dspec[j].getMonoMass();
          int repz2 = (int)round(mass2 / (pmz - FLASHDeconvHelperStructs::getChargeMass(is_positive_)));
          if (repz1 == repz2) { continue; }

          if (dspec[i].getChargeSNR(repz1) > dspec[j].getChargeSNR(repz2)) { continue; }

          is_overlap = true;

          break;
        }
        if (is_overlap) overlap_intensity[i] += pint;
      }
    }

    for (Size i = 0; i < dspec.size(); i++)
    {
      if (target_decoy_type != PeakGroup::TargetDecoyType::non_specific && dspec[i].getTargetDecoyType() != target_decoy_type) { continue; }
      if ((ms_level_ == 1 && dspec[i].getRepAbsCharge() < min_abs_charge_) || dspec[i].getRepAbsCharge() > max_abs_charge_) { continue; }

      double mul_factor = .5;
      if (! dspec[i].isTargeted() &&                                    // z1 != z2 &&
          overlap_intensity[i] >= dspec[i].getIntensity() * mul_factor) // If the overlapped intensity takes more than mul_factor * total intensity then
                                                                        // it is a peakgroup with a charge error. the smaller, the harsher
      {
        continue;
      }
      filtered_pg_vec.push_back(dspec[i]);
    }
    dspec.setPeakGroups(filtered_pg_vec);
  }

  void SpectralDeconvolution::removeOverlappingPeakGroups(DeconvolvedSpectrum& dspec, double tol, PeakGroup::TargetDecoyType target_decoy_type)
  {
    if (dspec.empty()) { return; }

    std::vector<PeakGroup> filtered_pg_vec; //
    filtered_pg_vec.reserve(dspec.size());

    double start_mass = dspec[0].getMonoMass();
    float local_max_SNR = -1.0;
    Size local_max_index = 0;
    Size last_local_max_index = dspec.size();

    for (Size i = 0; i < dspec.size(); i++)
    {
      double mass = dspec[i].getMonoMass();
      if (mass - start_mass > mass * tol * 1.5)
      {
        if (! dspec[local_max_index].isTargeted()) // targeted ones were already push_backed.
        {
          if ((target_decoy_type == PeakGroup::TargetDecoyType::non_specific || dspec[local_max_index].getTargetDecoyType() == target_decoy_type)
              && last_local_max_index != local_max_index)
            filtered_pg_vec.push_back(dspec[local_max_index]);
        }
        last_local_max_index = local_max_index;
        start_mass = mass;
        local_max_SNR = -1.0;
      }

      float snr = dspec[i].getSNR();

      if (local_max_SNR < snr)
      {
        local_max_SNR = snr;
        local_max_index = i;
      }
      if (dspec[i].isTargeted())
      {
        if (target_decoy_type == PeakGroup::TargetDecoyType::non_specific || dspec[i].getTargetDecoyType() == target_decoy_type)
          filtered_pg_vec.push_back(dspec[i]);
      }
    }

    if (local_max_SNR >= 0)
    {
      if (! dspec[local_max_index].isTargeted()) // targeted ones were already push_backed.
      {
        if ((target_decoy_type == PeakGroup::TargetDecoyType::non_specific || dspec[local_max_index].getTargetDecoyType() == target_decoy_type)
            && last_local_max_index != local_max_index)
          filtered_pg_vec.push_back(dspec[local_max_index]);
      }
    }

    dspec.setPeakGroups(filtered_pg_vec);
    std::vector<PeakGroup>().swap(filtered_pg_vec);
  }

  void SpectralDeconvolution::removeExcludedMasses_(DeconvolvedSpectrum& dspec, std::vector<double> excluded_masses) const
  {
    if (excluded_masses.empty()) return;
    std::vector<PeakGroup> filtered_pg_vec; //
    filtered_pg_vec.reserve(dspec.size());

    for (const auto& peak_group : dspec)
    {
      if (peak_group.getTargetDecoyType() != target_decoy_type_) continue;
      double delta = peak_group.getMonoMass() * tolerance_[ms_level_ - 1] * 2;
      auto upper = std::upper_bound(excluded_masses.begin(), excluded_masses.end(), peak_group.getMonoMass() + delta);
      bool exclude = false;
      while (! exclude)
      {
        if (upper != excluded_masses.end())
        {
          if (std::abs(*upper - peak_group.getMonoMass()) < delta) { exclude = true; }
          if (peak_group.getMonoMass() - *upper > delta) { break; }
        }
        if (upper == excluded_masses.begin()) { break; }
        --upper;
      }
      if (exclude) { continue; }
      filtered_pg_vec.push_back(peak_group);
    }
    dspec.setPeakGroups(filtered_pg_vec);
    std::vector<PeakGroup>().swap(filtered_pg_vec);
  }

  void SpectralDeconvolution::setTargetMasses(const std::vector<double>& masses, bool excluded)
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
    for (const auto& m : masses)
    {
      int start = 0;
      int end = 0;
      if (excluded) { end = (int)(avg_.getApexIndex(m) + avg_.getRightCountFromApex(m)); }
      for (int j = start; j <= end + 1; j++)
      {
        if (excluded) excluded_masses_.push_back(m + iso_da_distance_ * j);
        else
          target_mono_masses_.push_back(m + iso_da_distance_ * j);
      }
    }
  }

  void SpectralDeconvolution::setAveragine(const SpectralDeconvolution::PrecalculatedAveragine& avg)
  {
    avg_ = avg;
  }
  double SpectralDeconvolution::getMassFromMassBin_(Size mass_bin, double bin_mul_factor) const
  {
    return exp(getBinValue_(mass_bin, mass_bin_min_value_, bin_mul_factor));
  }

  double SpectralDeconvolution::getMzFromMzBin_(Size mass_bin, double bin_mul_factor) const
  {
    return exp(getBinValue_(mass_bin, mz_bin_min_value_, bin_mul_factor));
  }
} // namespace OpenMS