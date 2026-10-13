// Copyright (c) 2002-2025, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim, Jaekwan Kim $
// $Authors: Kyowon Jeong, Jihyung Kim, Jaekwan Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <cmath>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  namespace
  {
    /// round(a / b) given inv_b = 1 / b, without the division; false when the product is within 1e-7 of a rounding
    /// boundary (or large or non-finite), in which case the caller divides. Same reasoning as in PeakGroup.cpp.
    inline bool roundedQuotient(double a, double inv_b, int& rounded)
    {
      const double q = a * inv_b;
      if (! (std::fabs(q) < 1e8)) { return false; }
      const long long t = (long long)q;
      const long long f = t - (long long)(q < (double)t);
      const double frac = q - (double)f;
      if (! (frac > 1e-7 && frac < 1 - 1e-7 && std::fabs(frac - 0.5) > 1e-7)) { return false; }
      rounded = (int)(f + (long long)(frac > 0.5));
      return true;
    }
  } // namespace

  /// harmonic charge factors that will be considered for harmonic mass reduction.
  inline const std::vector<int> harmonic_charges_ {2, 3, 5, 7, 11};
  /// high and low charges are differently deconvolved. This value determines the (inclusive) threshold for low charge.
  inline const int low_charge_ = 10; // 10 inclusive
  inline const double tol_div_factor = 2.5; // use narrow tolerance for deconvolution and at the end use the input tolerance to filter out overlapping masses.

  /// whether one of the @p n values at @p p equals @p value; four at a time, as the 16-bit lanes of 64-bit words (w has a zero lane iff
  /// (w - 0x0001..) & ~w & 0x8000.. is not 0)
  static inline bool anyEquals_(const unsigned short* p, const long n, const int value)
  {
    if (value < 0 || value > (int)std::numeric_limits<unsigned short>::max()) { return false; }
    const std::uint64_t lanes = 0x0001000100010001ULL;
    const std::uint64_t pattern = lanes * (std::uint64_t)value;
    long i = 0;
    for (; i + 4 <= n; i += 4)
    {
      std::uint64_t w;
      std::memcpy(&w, p + i, sizeof(w));
      w ^= pattern;
      if (((w - lanes) & ~w & (lanes << 15)) != 0) { return true; }
    }
    for (; i < n; i++)
    {
      if (p[i] == value) { return true; }
    }
    return false;
  }

  /// a hint to the processor to fetch the cache line of @p p for writing (no effect on any result)
  static inline void prefetchForWrite_(const void* p)
  {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(p, 1);
#else
    (void)p;
#endif
  }

  SpectralDeconvolution::SpectralDeconvolution(): DefaultParamHandler("SpectralDeconvolution")
  {
    defaults_.setValue(
      "tol", DoubleList {10.0, 10.0},
      "PPM Tolerance for MS1, 2, ...: Specify the tolerance values in parts per million (ppm) using the this option. For example, -tol 10.0 5.0 sets "
      "the tolerance to 10.0 ppm for MS1 and 5.0 ppm for MS2. "
      "If a negative value (e.g., -1) is given for an MS level, its tolerance is determined from the data: spectra of the level are deconvolved at 30 ppm "
      "(every 64th, then more until 20,000 peak errors or every 4th spectrum), and the tolerance is set so that +-0.8 x tol (the window in which a peak is "
      "taken) holds 98% of the peak errors of the masses of Qscore >= 0.9, after the random peaks are subtracted (no fit). If the errors reach beyond "
      "a third of the window, at 100 and then 300 ppm. With fewer than 1,000 errors or no error distribution standing out of the random peaks, "
      "10 ppm. The determined tolerance is logged.");

    defaults_.setValue("min_mass", 50.0, "Minimum mass (Da)");
    defaults_.setValue("max_mass", 100000.0, "Maximum mass (Da)");

    defaults_.setValue("min_charge", 1,
                       "Minimum charge state for MS1 spectra (can be negative for negative mode). For MSn (n > 1), minimum charge is fixed to 1.");
    defaults_.setValue(
      "max_charge", 100,
      "Maximum charge state for spectra (can be negative for negative mode). Apart from min_charge, this option applies to all MS levelss.");

    defaults_.setValue("precursor_charge", 0,
                       "Charge state of the target precursor. All precursor charge for MSn (n > 1) is fixed to this value. When precursor m/z is "
                       "provided within the input mzML file or is specified using precursor_mz option, "
                       "precursor mass is calculated accordingly.");
    defaults_.setMinInt("precursor_charge", 0);

    defaults_.setValue(
      "precursor_mz", 0.0,
      "Target precursor m/z value. This option must be used in combination with the -precursor_charge option; otherwise, it will be ignored. "
      "If the -precursor_charge option is specified but this option is not, the precursor m/z value from the input mzML file will be used as the "
      "default.");
    defaults_.setMinFloat("precursor_mz", 0.0);
    defaults_.addTag("precursor_mz", "advanced");

    defaults_.setValue("min_cos", DoubleList {.85, .85},
                       "Cosine similarity thresholds between avg. and observed isotope pattern for MS1, 2, ...: e.g., -min_cos 0.3 0.6 to specify 0.3 "
                       "and 0.6 for MS1 and MS2, respectively.");
    defaults_.addTag("min_cos", "advanced");
    defaults_.setMinFloat("min_cos", 0.5);
    defaults_.setMaxFloat("min_cos", 1.0);

    defaults_.setValue("min_snr", DoubleList {.25, .25},
                       "Minimum charge SNR (the SNR of the isotope pattern of a specific charge) thresholds for MS1, 2, ...: e.g., -min_snr 1.0 0.6 to "
                       "specify 1.0 and 0.6 for MS1 and MS2, respectively.");
    defaults_.addTag("min_snr", "advanced");
    defaults_.setMinFloat("min_snr", 0.25);

    defaults_.setValue("allowed_isotope_error", 0,
                       "Tolerance for isotope index errors when calculating FDR. For instance, setting a value of 2 permits the inclusion of up to 2 "
                       "isotope errors as valid matches. Beta version.");
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

    // Get precursor from spectrum if available, to keep its isolation window and other metadata
    Precursor precursor;
    if (!spec.getPrecursors().empty())
    {
      precursor = spec.getPrecursors()[0];
    }
    // An explicitly specified target m/z takes precedence over the one from the spectrum
    if (target_precursor_mz_ > 0)
    {
      precursor.setMZ(target_precursor_mz_);
    }

    double target_precursor_mass
      = (precursor.getMZ() - FLASHHelperClasses::getChargeMass(target_precursor_charge_ > 0)) * std::abs(target_precursor_charge_);
    precursor.setCharge(target_precursor_charge_);
    PeakGroup precursorPeakGroup(1, std::abs(target_precursor_charge_), target_precursor_charge_ > 0);
    precursorPeakGroup.push_back(FLASHHelperClasses::LogMzPeak());
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
    // First prepare for decoy runs. if it is noisy decoy, change the distance between the isotopes
    iso_da_distance_ = target_decoy_type_ == PeakGroup::noise_decoy
                         ? Constants::ISOTOPE_MASSDIFF_55K_U * noise_iso_delta_
                         : Constants::ISOTOPE_MASSDIFF_55K_U;

    // excluded peak masses or mass bins or monoisotopic masses for signal decoy runs
    excluded_peak_masses_for_decoy_runs_.clear();
    excluded_mass_bins_for_decoy_runs_.reset();
    excluded_masses_for_decoy_runs_.clear();

    if (target_decoy_type_ == PeakGroup::signal_decoy) // charge decoy
    {
      for (const auto& pg : *target_dspec_for_decoy_calculation_) // pg are the target peak groups from normal spectrum deconvolution
      {
        for (int i = -0; i <= 0; i++)
        {
          excluded_masses_for_decoy_runs_.push_back(pg.getMonoMass() + i * iso_da_distance_);
          excluded_peak_masses_for_decoy_runs_.push_back(pg.getMonoMass() + activeAveragine_().getAverageMassDelta(pg.getMonoMass()) + i * iso_da_distance_);
        }
      }
      std::sort(excluded_masses_for_decoy_runs_.begin(), excluded_masses_for_decoy_runs_.end());
      excluded_masses_for_decoy_runs_.erase( unique(excluded_masses_for_decoy_runs_.begin(), excluded_masses_for_decoy_runs_.end() ),
                                            excluded_masses_for_decoy_runs_.end() );

      std::sort(excluded_peak_masses_for_decoy_runs_.begin(), excluded_peak_masses_for_decoy_runs_.end());
      excluded_peak_masses_for_decoy_runs_.erase( unique(excluded_peak_masses_for_decoy_runs_.begin(), excluded_peak_masses_for_decoy_runs_.end() ),
                                                 excluded_peak_masses_for_decoy_runs_.end() );
    }

    ms_level_ = spec.getMSLevel();
    deconvolved_spectrum_ = DeconvolvedSpectrum(scan_number);
    deconvolved_spectrum_.setOriginalSpectrum(spec);

    // for noise decoy, we exclude all raw peaks used for the normal run
    if (target_decoy_type_ == PeakGroup::noise_decoy)
    {
      std::set<double> signal_mzs;
      for (const auto& pg : *target_dspec_for_decoy_calculation_)
      {
        for (const auto& p : pg)
          signal_mzs.insert(p.mz);
      }

      auto nspec = spec;
      nspec.clear(false);

      for (const auto& p : spec)
      {
        if (signal_mzs.find(p.getMZ()) != signal_mzs.end()) continue;
        nspec.push_back(p);
      }
      deconvolved_spectrum_.setOriginalSpectrum(nspec);
    }

    // here register targeted peak mzs etc.
    // for MSn (n>1) register precursor peak and peak group.
    if (ms_level_ > 1)
    {
      for (const auto& precursor : deconvolved_spectrum_.getOriginalSpectrum().getPrecursors())
      {
        for (const auto& activation_method : precursor.getActivationMethods())
        {
          deconvolved_spectrum_.setActivationMethod(activation_method);
          if (deconvolved_spectrum_.getActivationMethod() == Precursor::ActivationMethod::HCID) { deconvolved_spectrum_.setActivationMethod(Precursor::ActivationMethod::HCD); }
          break;
        }
        deconvolved_spectrum_.setPrecursor(precursor);
      }

      // a target precursor m/z needs its charge (help of precursor_mz): without one it built a precursor of mass 0, which capped the MSn masses
      if (target_precursor_charge_ != 0) { setTargetPrecursorCharge_(); }

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
        else
        {

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
      bin_mul_factors_.push_back(1.0 / j * tol_div_factor);
    }

    min_isotope_cosine_ = param_.getValue("min_cos");
    min_snr_ = param_.getValue("min_snr");

    target_precursor_mz_ = param_.getValue("precursor_mz");
    target_precursor_charge_ = param_.getValue("precursor_charge");
  }

  const FLASHHelperClasses::PrecalculatedAveragine& SpectralDeconvolution::getAveragine()
  {
    return activeAveragine_();
  }

  void SpectralDeconvolution::calculateAveragine(const bool use_RNA_averagine)
  {
    CoarseIsotopePatternGenerator generator(300);

    auto iso = use_RNA_averagine ? generator.estimateFromRNAWeight(current_max_mass_) : generator.estimateFromPeptideWeight(current_max_mass_);
    iso.trimRight(0.0001 * iso.getMostAbundant().getIntensity());
    auto max_isotope = std::max(200, (int)iso.size());

    generator.setMaxIsotope(max_isotope);
    avg_ = FLASHHelperClasses::PrecalculatedAveragine(50, current_max_mass_, 25, generator, use_RNA_averagine,
                                                      target_decoy_type_ == PeakGroup::noise_decoy ? Constants::ISOTOPE_MASSDIFF_55K_U * noise_iso_delta_ : -1);
    avg_.setMaxIsotopeIndex((int)(max_isotope - 1));
    shared_avg_.reset();
  }

  // generate filters
  void SpectralDeconvolution::setFilters_()
  {
    universal_pattern_.clear();
    int charge_range = current_max_charge_;
    for (int i = 0; i < charge_range; i++)
    {
      universal_pattern_.push_back(-log(i + 1));
    }

    harmonic_pattern_matrix_ = Matrix<double>(harmonic_charges_.size(), charge_range, .0);

    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      int hc = harmonic_charges_[k];
      int n = hc / 2;
      // int m = (1 + hc) / 2;
      for (int i = 0; i < charge_range; i++)
      {
        double a = i > 0 ? exp(-universal_pattern_[i - 1]) : 0;
        double b = exp(-universal_pattern_[i]);
        harmonic_pattern_matrix_.setValue(k, i, -log(b - (b - a) * n / hc));
      }
    }
  }

  // Generate uncharged log mz transformated peaks
  void SpectralDeconvolution::updateLogMzPeaks_()
  {
    const auto& spec = deconvolved_spectrum_.getOriginalSpectrum();
    // the log m/z of every peak once per spectrum; the peak recruitment of the scoring reuses them
    original_log_mzs_.resize(spec.size());
    for (Size i = 0; i < spec.size(); i++)
    {
      original_log_mzs_[i] = FLASHHelperClasses::getLogMz(spec[i].getMZ(), is_positive_);
    }
    original_nearest_peaks_.build(spec);
    log_mz_peaks_.clear();
    log_mz_peaks_.reserve(spec.size());

    for (Size i = 0; i < spec.size(); i++)
    {
      const auto& peak = spec[i];
      if (peak.getIntensity() <= 0) //
      {
        continue;
      }
      log_mz_peaks_.emplace_back(peak, is_positive_, original_log_mzs_[i]);
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
  void SpectralDeconvolution::binLogMzPeaks_(const Size bin_number, std::vector<float>& binned_log_mz_peak_intensities)
  {
    // reuse the storage of the previous spectrum (all bits cleared)
    binned_log_mz_peaks_.resize(bin_number);
    binned_log_mz_peaks_.reset();
    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];

    for (const auto& p : log_mz_peaks_)
    {
      Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_mul_factor);
      if (bi >= bin_number) { break; }
      binned_log_mz_peaks_.set(bi);

      binned_log_mz_peak_intensities[bi] += p.intensity;
    }
  }

  // Find candidate mass bins from the current spectrum. The runtime of FLASHDeconv is determined by this function.
  void SpectralDeconvolution::updateCandidateMassBins_(const std::vector<float>& mz_intensities, std::vector<std::pair<Size, long>>& chain_redirects)
  { //
    Size mz_bin_index = binned_log_mz_peaks_.find_first();
    auto& mz_bin_index_reverse = mz_bin_index_reverse_pool_; // pooled: refilled here
    mz_bin_index_reverse.clear();
    mz_bin_index_reverse.reserve(binned_log_mz_peaks_.count());
    // invert mz bins so charges are counted from small to large given a mass

    while (mz_bin_index != binned_log_mz_peaks_.npos)
    {
      mz_bin_index_reverse.push_back(mz_bin_index);
      mz_bin_index = binned_log_mz_peaks_.find_next(mz_bin_index);
    }
    size_t h_charge_size = harmonic_charges_.size();
    long bin_end = (long)binned_log_masses_.size();

    // The state per mass bin is kept from spectrum to spectrum (mass_bin_scratch_), so that the three mass-bin sized arrays (~18 bytes per
    // bin, 46-68 MB per spectrum at 5 ppm) are neither allocated nor filled per spectrum:
    // - to calculate continuous charges, the previous charge value per mass should be stored. prev_charges holds it as a tag,
    //   charge_tag_base + charge index, with a charge_tag_base above all tags of the previous spectra: a lower value is "no charge"
    //   (no_charge resets a bin). The tags are reset (to 0) only when they would overflow, every few hundred spectra.
    // - the rest of the state per mass bin in one struct (one cache line per peak and charge, not one per array): not just charges but
    //   intensities are stored to see the intensity fold change, how many peaks are present to support that mass bin, and (see below) the
    //   unrounded position of the last peak and the offset to the chain's bin. prev_charges stays apart for the chain search. A ChainState is
    //   read only where prev_charges holds a tag of this spectrum, which is written together with all of its fields, so it is not reset.
    // - mass_intensities is all zero between spectra: updateMassBins_ zeroes the bins written here (written_mass_bins) after filterMassBins_.
    auto& scratch = mass_bin_scratch_;
    auto& prev_charges = scratch.prev_charges;
    auto& chain_states = scratch.chain_states;
    auto& mass_intensities = scratch.mass_intensities;
    if (prev_charges.size() < (Size)bin_end)
    {
      prev_charges.resize(bin_end, 0);
      chain_states.resize(bin_end);
    }
    int charge_tag_base = 0; // 0: as without tags (the charge index itself; a charge range too wide for tags)
    unsigned short no_charge = 0;
    if (current_max_charge_ + 2 > (int)std::numeric_limits<unsigned short>::max())
    {
      no_charge = (unsigned short)(current_max_charge_ + 2);
      std::fill(prev_charges.begin(), prev_charges.begin() + bin_end, no_charge);
      std::fill(chain_states.begin(), chain_states.begin() + bin_end, MassBinScratch_::ChainState {});
      scratch.charge_tag_base = 0;
    }
    else
    {
      if (scratch.charge_tag_base < 1 || scratch.charge_tag_base + current_max_charge_ > (int)std::numeric_limits<unsigned short>::max() + 1)
      {
        std::fill(prev_charges.begin(), prev_charges.end(), (unsigned short)0);
        scratch.charge_tag_base = 1;
      }
      charge_tag_base = scratch.charge_tag_base;
      scratch.charge_tag_base += std::max(current_max_charge_, 0);
    }

    if (! scratch.mass_intensities_clean) { std::fill(mass_intensities.begin(), mass_intensities.end(), .0f); }
    if (mass_intensities.size() < (Size)bin_end) { mass_intensities.resize(bin_end, .0f); }
    scratch.mass_intensities_clean = false;
    auto& written_mass_bins = scratch.written_mass_bins;
    written_mass_bins.clear();

    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];
    // per charge index: the intensity ratio between consecutive charges should not exceed the factor; for possible harmonics, half of it
    std::vector<float> factors(current_max_charge_), hfactors(current_max_charge_);
    for (int j = 0; j < current_max_charge_; j++)
    {
      const int abs_charge = (j + 1);
      const float highest_factor = 10.0f;
      factors[j] = abs_charge <= low_charge_ ? highest_factor : (highest_factor / 2 + highest_factor / 2 * low_charge_ / (float)abs_charge);
      hfactors[j] = factors[j] / 2.0f;
    }

    // A high charge (> low_charge_) is accepted only as part of a chain of consecutive charges of one mass, and that continuity is decided
    // on the unrounded masses, not by mass bin identity. A mass bin is the m/z bin of the peak plus the rounded charge offset, so two peaks
    // of one isotope at consecutive charges land in the same mass bin only by chance (two roundings and the m/z error of each peak): on
    // Eclip_1480, 70% of such pairs of real high-charge envelopes (median difference 1.5 ppm) fell in different 2-ppm bins, which broke the
    // chain of their mass. A high-charge peak continues the chain of the previous charge whose last peak is within the tolerance
    // (tol_div_factor bins) in unrounded mass, the closest one, and adds to the mass bin of that chain. A low-charge peak (and every peak in
    // MSn) is accepted by its own isotope peaks; there the chain keeps the mass bin identity: with the unrounded rule at every charge,
    // coincidental chains of unrelated peaks took peaks from their own masses (Eclip_1480 MSn: 5% fewer complementary fragment pairs).
    // Per mass bin (ChainState): the unrounded position (in bins) of its last peak relative to the bin and the offset from the bin to the
    // chain's bin.
    const double chain_tol = tol_div_factor;                // tolerance in mass bins (= the input ppm tolerance)
    const long chain_search = (long)(chain_tol + 2.0);      // + the two rounding errors (each <= 1/2 bin) of each of the two peaks
    std::vector<double> exact_mass_offsets(current_max_charge_); // unrounded binned_universal_pattern_
    for (int j = 0; j < current_max_charge_; j++)
    {
      exact_mass_offsets[j] = (mz_bin_min_value_ - universal_pattern_[j] - mass_bin_min_value_) * bin_mul_factor;
    }
    // unrounded position of each m/z bin (in bins): its most intense peak
    auto& exact_mz_positions = exact_mz_positions_pool_; // pooled: refilled here
    exact_mz_positions.assign(mz_bin_index_reverse.size(), .0);
    {
      Size k = 0;
      float max_intensity = -1;
      for (const auto& p : log_mz_peaks_)
      {
        Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_mul_factor);
        while (k < mz_bin_index_reverse.size() && mz_bin_index_reverse[k] < bi)
        {
          k++;
          max_intensity = -1;
        }
        if (k >= mz_bin_index_reverse.size()) { break; }
        if (mz_bin_index_reverse[k] != bi || p.intensity <= max_intensity) { continue; }
        max_intensity = p.intensity;
        exact_mz_positions[k] = (p.logMz - mz_bin_min_value_) * bin_mul_factor;
      }
    }
    chain_redirects.clear();

    // The state of the mass bins of one peak (one per charge) is spread over the whole mass bin range: a cache miss per peak and charge.
    // While a peak is processed, the state of the next peak's mass bin at the same charge is fetched (prev_charges over the chain search
    // window and the ChainState written), one peak ahead.
    auto prefetch_mass_bin = [&](const long next_mz_bin_index, const int j) {
      const long b = next_mz_bin_index + binned_universal_pattern_[j];
      if (b < 0 || b >= bin_end) { return; }
      const long search = (j > 0 && j + 1 > low_charge_ && ms_level_ == 1) ? chain_search : 0;
      prefetchForWrite_(prev_charges.data() + std::max(0L, b - search));
      prefetchForWrite_(prev_charges.data() + std::min(bin_end - 1, b + search));
      prefetchForWrite_(chain_states.data() + b);
      prefetchForWrite_(reinterpret_cast<const char*>(chain_states.data() + b) + sizeof(MassBinScratch_::ChainState) - 1);
    };

    // traverse from right to left in the log mz bin space
    for (auto iter = mz_bin_index_reverse.rbegin(); iter < mz_bin_index_reverse.rend(); iter++)
    {
      mz_bin_index = *iter;
      const long next_mz_bin_index = iter + 1 < mz_bin_index_reverse.rend() ? (long)*(iter + 1) : -1; // the next peak (-1: none)
      const double exact_mz_position = exact_mz_positions[mz_bin_index_reverse.rend() - iter - 1];
      const float intensity = mz_intensities[mz_bin_index];
      const double log_mz = getBinValue_(mz_bin_index, mz_bin_min_value_, bin_mul_factor); // uncharged log mz;
      const double mz = exp(log_mz);                                                       // uncharged mz
      const double iso_div_by_mz = iso_da_distance_ / mz;
      // iso_div_by_mz / (abs_charge - 1), / abs_charge and / (abs_charge + 1) for the charge iso_div_charge (the isotope check below); for the
      // next charge only the last one is new, and for d = -1 they are negated (the same values as of d * iso_div_by_mz)
      int iso_div_charge = -1;
      double iso_div_lower = 0, iso_div_this = 0, iso_div_upper = 0;
      // scan through charges
      for (int j = 0; j < current_max_charge_; j++) // loop over all charges
      {
        if (next_mz_bin_index >= 0) { prefetch_mass_bin(next_mz_bin_index, j); }
        // mass is given by shifting by binned_universal_pattern_[j]
        const long mass_bin_index = (long)mz_bin_index + binned_universal_pattern_[j];

        if (mass_bin_index < 0) { continue; }
        if (mass_bin_index >= bin_end)
        {
          // the next peak (lower m/z) may have mass bins at higher charges
          for (int k = j + 1; next_mz_bin_index >= 0 && k < current_max_charge_ && next_mz_bin_index + binned_universal_pattern_[k] < bin_end; k++)
          {
            prefetch_mass_bin(next_mz_bin_index, k);
          }
          break;
        }

        if (! excluded_mass_bins_for_decoy_runs_.empty() && excluded_mass_bins_for_decoy_runs_[mass_bin_index]) { continue; }

        const int abs_charge = (j + 1);
        bool pass_first_check = false;

        // intensity ratio between consecutive charges should not exceed the factor.
        const float factor = factors[j];
        // intensity ratio between consecutive charges for possible harmonic should be within this factor
        const float hfactor = hfactors[j];
        auto get_intensity_ratio = [&](const float prev_int) {
          float r = prev_int <= 0 ? (factor + 1) : (intensity / prev_int);
          return r < 1 ? 1.0f / r : r;
        };

        // the chain of the previous charge this peak continues: the mass bin (for a high charge within chain_search bins, else its own) whose
        // last peak is at the previous charge and within chain_tol in unrounded mass; of those, one within the intensity ratio factor first,
        // then the closest
        const float residual = (float)(exact_mz_position + exact_mass_offsets[j] - (double)mass_bin_index);
        long chain_bin = -1;
        double chain_dist = 0;
        bool chain_ratio_ok = false;
        if (j > 0)
        {
          const long search = (abs_charge > low_charge_ && ms_level_ == 1) ? chain_search : 0; // MSn keeps the mass bin identity
          const int prev_charge_tag = charge_tag_base + j - 1;
          const long first_s = std::max(0L, mass_bin_index - search), last_s = std::min(bin_end - 1, mass_bin_index + search);
          // most windows hold no peak of the previous charge (Eclip_1480 MS1: 0.18 per window of 9 bins): all of them are compared at once first
          const bool any_chain = anyEquals_(prev_charges.data() + first_s, last_s - first_s + 1, prev_charge_tag);
          for (long s = first_s; any_chain && s <= last_s; s++)
          {
            if (prev_charges[s] != prev_charge_tag) { continue; }
            const double d = std::abs((double)(mass_bin_index - s) + residual - chain_states[s].prev_residual);
            if (d > chain_tol) { continue; }
            const bool ratio_ok = get_intensity_ratio(chain_states[s].prev_intensity) <= factor;
            if (chain_bin < 0 || (ratio_ok && ! chain_ratio_ok) || (ratio_ok == chain_ratio_ok && d < chain_dist))
            {
              chain_bin = s;
              chain_dist = d;
              chain_ratio_ok = ratio_ok;
            }
          }
        }
        const bool charge_continous = chain_bin >= 0;
        // the mass bin this peak adds to: the bin of the chain it continues, else its own
        const long target_bin = charge_continous ? chain_bin + chain_states[chain_bin].chain_offset : mass_bin_index;
        unsigned short spc = charge_continous ? chain_states[chain_bin].support_peak_count : 0;
        // read only for a continued chain (below); not loaded otherwise (a cache miss per peak and charge)
        const float prev_intensity = charge_continous ? chain_states[chain_bin].prev_intensity : .0f;

        float support_peak_intensity = 0;

        // if charge not continuous or intensity ratio (between current and previous charges) is too high reset support_peak_count
        if (! charge_continous || get_intensity_ratio(prev_intensity) > factor) { spc = 0; }
        else
        {
          pass_first_check = true;
          if (spc == 0 && abs_charge > low_charge_) { support_peak_intensity = prev_intensity; }
        }

        // for low charges, check isotope peak presence. In MSn spectra for every charge: a fragment ion mostly appears at a single
        // charge state, so the consecutive-charge check above rarely passes for it (Eclip_1480: of 588 fragments of the TopPIC-only
        // PrSMs seen at one charge above 10, 5 were reported without this).
        if (! pass_first_check && (abs_charge <= low_charge_ || ms_level_ > 1))
        {
          // find the next isotope peak(s) and try to avoid harmonic masses
          const int last_d = activeAveragine_().getApexIndex(mz * abs_charge) > 0 ? -1 : 1; // once, not per loop condition check
          if (iso_div_charge == abs_charge - 1)
          {
            iso_div_lower = iso_div_this;
            iso_div_this = iso_div_upper;
            iso_div_upper = iso_div_by_mz / (abs_charge + 1);
          }
          else if (iso_div_charge != abs_charge)
          {
            iso_div_lower = abs_charge > 1 ? iso_div_by_mz / (abs_charge - 1) : 0;
            iso_div_this = iso_div_by_mz / abs_charge;
            iso_div_upper = iso_div_by_mz / (abs_charge + 1);
          }
          iso_div_charge = abs_charge;
          bool harmonic_found = false; // a harmonic peak between this peak and an isotope peak of it
          for (int d = 1; d >= last_d; d -= 2)
          {
            bool iso_exist = false;
            int next_iso_bin = 0;
            const int nib = (int)getBinNumber_(log_mz + (d > 0 ? iso_div_this : -iso_div_this), mz_bin_min_value_, bin_mul_factor);
            const int nibr
              = abs_charge > 1 ? (int)getBinNumber_(log_mz + (d > 0 ? iso_div_lower : -iso_div_lower), mz_bin_min_value_, bin_mul_factor) : 0;
            const int nibl = (int)getBinNumber_(log_mz + (d > 0 ? iso_div_upper : -iso_div_upper), mz_bin_min_value_, bin_mul_factor);
            if (std::abs(nib - nibr) < tol_div_factor
                || std::abs(nib - nibl)
                     < tol_div_factor) // if different charges are not distinguishable, we ignore. Not informative and the source of the errors.
              break;

            // in MSn the neighbouring isotope peak is looked for within +-2 bins (~0.8 x the input tolerance), not within +-1 bin
            // (~0.4 x, 2 ppm at 5 ppm): on Eclip_1480 the +-1 bin window lost 11% of the real fragment envelopes (fragments of the
            // proteoforms TopPIC found and FLASHTnT did not), mostly at charges 1-4
            const int iso_t = ms_level_ > 1 ? (int)tol_div_factor : 1;
            for (int t = -iso_t; t <= iso_t; t++)
            {
              int nibt = nib + t;
              if (std::abs(nibt - (int)mz_bin_index) >= tol_div_factor && nibt > 0 && nibt < (int)binned_log_mz_peaks_.size()
                  && binned_log_mz_peaks_[nibt])
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

                const int hdiff = (int)round((double)((int)next_iso_bin - (int)mz_bin_index)) / hc * (hc / 2);
                const int next_harmonic_iso_bin = (int)mz_bin_index + hdiff;
                // check if there are harmonic peaks between the current peak and the next isotope peak.

                // no perfect filtration. Just obvious ones are filtered out by checking if a peak is in the harmonic position and the intensity ratio
                // is within two folds from the current peak (specified by mz_bin_index)
                if (std::abs(next_harmonic_iso_bin - (int)mz_bin_index) >= tol_div_factor && next_harmonic_iso_bin >= 0
                    && next_harmonic_iso_bin < (int)binned_log_mz_peaks_.size() && binned_log_mz_peaks_[next_harmonic_iso_bin]
                    && mz_intensities[next_harmonic_iso_bin] > h_threshold / 2 && mz_intensities[next_harmonic_iso_bin] < h_threshold * 2)
                {
                  harmonic_cntr++;
                  harmonic_found = true;
                }

                if (harmonic_cntr > 0) { pass_first_check = false; }
              }
            }
            if (pass_first_check)
            {
              support_peak_intensity += mz_intensities[next_iso_bin];
            }
          }
          pass_first_check &= ! harmonic_found;
        }

        if (pass_first_check)
        {
          if (charge_continous) // check harmonic artifacts for high charge ranges
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
            int min_dis = (int)(tol_div_factor + 1);
            float max_harmonic_intensity = 0; // the largest of the harmonic intensities of the harmonic charges
            for (size_t k = 0; k < h_charge_size; k++)
            {
              if (ms_level_ > 1 && harmonic_charges_[k] * abs_charge > current_max_charge_) break;
              float harmonic_intensity = 0;
              for (int t = -(int)tol_div_factor; t <= (int)tol_div_factor; t++)
              {
                long hmz_bin_index = mass_bin_index - binned_harmonic_patterns.getValue(k, j) + t;
                if (hmz_bin_index > 0 && hmz_bin_index != (long)mz_bin_index && hmz_bin_index < (int)binned_log_mz_peaks_.size()
                    && binned_log_mz_peaks_[hmz_bin_index])
                {
                  float h_intensity = mz_intensities[hmz_bin_index];
                  if (h_intensity > low_threshold && h_intensity < high_threshold)
                  {
                    is_harmonic = true;
                    if (std::abs(t) < min_dis)
                    {
                      harmonic_intensity = h_intensity;
                      min_dis = abs(t);
                    }
                  }
                }
              }
              max_harmonic_intensity = std::max(max_harmonic_intensity, harmonic_intensity);
            }

            if (! is_harmonic) // if it is not harmonic
            {
              if (mass_intensities[target_bin] == 0) { written_mass_bins.push_back((Size)target_bin); }
              mass_intensities[target_bin] += intensity + support_peak_intensity;

              if (! binned_log_masses_[target_bin])
              {
                spc++;
                if (spc >= min_support_peak_count_ || spc >= abs_charge / 2) { binned_log_masses_[target_bin] = true; }
              }
            }
            else // if harmonic
            {
              if (mass_intensities[target_bin] == 0) { written_mass_bins.push_back((Size)target_bin); }
              mass_intensities[target_bin] -= max_harmonic_intensity;
              if (spc > 0) { spc--; }
            }
          }
          else if (abs_charge <= low_charge_ || ms_level_ > 1) // for low charge (every charge in MSn), include the mass if isotope is present
          {
            if (mass_intensities[target_bin] == 0) { written_mass_bins.push_back((Size)target_bin); }
            mass_intensities[target_bin] += intensity + support_peak_intensity;
            if (! binned_log_masses_[target_bin])
            {
              spc++;
              binned_log_masses_[target_bin] = true;
            }
          }
        }
        // the chain continues from this peak (its state is kept in the peak's own mass bin); the chain it continued cannot be continued again
        if (charge_continous && chain_bin != mass_bin_index) { prev_charges[chain_bin] = no_charge; }
        if (target_bin != mass_bin_index) { chain_redirects.emplace_back(mz_bin_index * current_max_charge_ + j, target_bin); }
        auto& chain_state = chain_states[mass_bin_index];
        chain_state.prev_intensity = intensity;
        prev_charges[mass_bin_index] = (unsigned short)(charge_tag_base + j);
        chain_state.support_peak_count = spc;
        chain_state.prev_residual = residual;
        chain_state.chain_offset = (short)(target_bin - mass_bin_index);
      }
    }
  }

  // Subfunction of updateMassBins_. If a peak corresponds to multiple masses, only one mass is selected for the peak based on intensities.
  // mass level harmonic check is also performed in this function
  // it also outputs the charge range of each mass bin
  std::vector<SpectralDeconvolution::MassBinChargeRange_> SpectralDeconvolution::filterMassBins_(const std::vector<float>& mass_intensities,
                                                                                                 std::vector<std::pair<Size, long>>& chain_redirects)
  {
    // chain_redirects in ascending order. updateCandidateMassBins_ adds them by m/z bin in descending order and, within one, by charge in
    // ascending order (unique keys): the runs of the m/z bins in reverse order are sorted, with no sort (std::sort for any other order)
    if (! std::is_sorted(chain_redirects.begin(), chain_redirects.end()))
    {
      std::reverse(chain_redirects.begin(), chain_redirects.end());
      for (Size run_start = 0; run_start < chain_redirects.size();)
      {
        Size run_end = run_start + 1; // the run is now descending
        while (run_end < chain_redirects.size() && chain_redirects[run_end].first < chain_redirects[run_end - 1].first) { run_end++; }
        std::reverse(chain_redirects.begin() + run_start, chain_redirects.begin() + run_end);
        run_start = run_end;
      }
      if (! std::is_sorted(chain_redirects.begin(), chain_redirects.end())) { std::sort(chain_redirects.begin(), chain_redirects.end()); }
    }
    Size redirect_index = 0;
    // (mass bin, charge index) per selection, reduced to a charge range per mass bin at the end: no table over all mass bins to fill
    std::vector<std::pair<Size, int>> selections;
    Size mz_bin_index = binned_log_mz_peaks_.find_first();
    long bin_size = (long)binned_log_masses_.size();

    // the mass bins selected by updateCandidateMassBins_; binned_log_masses_ (all clear) collects the ones selected here
    boost::dynamic_bitset<> candidate_bins(binned_log_masses_.size());
    candidate_bins.swap(binned_log_masses_);
    const bool targeted = ! target_mono_masses_.empty();

    const int select_top_N = 2; // select top N charges per peak. We allow up to 2 just to consider frequent coelution.
    std::vector<long> max_indices(select_top_N, -1);
    std::vector<int> max_intensity_abs_charge_ranges(select_top_N, -1);

    while (mz_bin_index != binned_log_mz_peaks_.npos)
    {
      std::fill(max_indices.begin(), max_indices.end(), -1);
      std::fill(max_intensity_abs_charge_ranges.begin(), max_intensity_abs_charge_ranges.end(), -1);

      float max_intensity = 0;

      for (int j = 0; j < current_max_charge_; j++)
      {
        long mass_bin_index = (long)mz_bin_index + binned_universal_pattern_[j];

        if (mass_bin_index < 0) { continue; }
        if (mass_bin_index >= bin_size) { break; }

        // a peak that continued a consecutive-charge chain counts for the mass bin of the chain (updateCandidateMassBins_)
        const Size key = mz_bin_index * current_max_charge_ + j;
        while (redirect_index < chain_redirects.size() && chain_redirects[redirect_index].first < key) { redirect_index++; }
        if (redirect_index < chain_redirects.size() && chain_redirects[redirect_index].first == key) { mass_bin_index = chain_redirects[redirect_index].second; }

        if (targeted && target_mass_bins_[mass_bin_index])
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
          if (! candidate_bins[mass_bin_index]) { continue; }

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
          selections.emplace_back((Size)max_index, max_intensity_abs_charge_range);
          binned_log_masses_[max_index] = true;
        }
      }
      mz_bin_index = binned_log_mz_peaks_.find_next(mz_bin_index);
    }

    // per selected mass bin (the set bits of binned_log_masses_, in increasing order) the min and max of its charge indices: the selections
    // sorted, as one 64-bit key each (bin, then charge index) when they fit
    std::vector<MassBinChargeRange_> mass_bin_charge_ranges;
    mass_bin_charge_ranges.reserve(selections.size());
    auto add_selection = [&mass_bin_charge_ranges](const Size bin, const int charge_index) {
      if (mass_bin_charge_ranges.empty() || mass_bin_charge_ranges.back().bin != bin) { mass_bin_charge_ranges.push_back({bin, charge_index, charge_index}); }
      else { mass_bin_charge_ranges.back().max_charge_index = charge_index; }
    };
    if ((std::uint64_t)bin_size <= ((std::uint64_t)1 << 32))
    {
      std::vector<std::uint64_t> keys;
      keys.reserve(selections.size());
      for (const auto& [bin, charge_index] : selections) { keys.push_back(((std::uint64_t)bin << 32) | (std::uint32_t)charge_index); }
      std::sort(keys.begin(), keys.end());
      for (const std::uint64_t key : keys) { add_selection((Size)(key >> 32), (int)(std::uint32_t)key); }
    }
    else
    {
      std::sort(selections.begin(), selections.end());
      for (const auto& [bin, charge_index] : selections) { add_selection(bin, charge_index); }
    }
    return mass_bin_charge_ranges;
  }

  // update mass bins which will be used to select peaks in the input spectrum...
  std::vector<SpectralDeconvolution::MassBinChargeRange_> SpectralDeconvolution::updateMassBins_(const std::vector<float>& mz_intensities)
  {
    std::vector<std::pair<Size, long>> chain_redirects;
    updateCandidateMassBins_(mz_intensities, chain_redirects); // universal pattern matching

    auto mass_bin_charge_ranges = filterMassBins_(mass_bin_scratch_.mass_intensities, chain_redirects);
    // mass_intensities all zero again for the next spectrum: the bins written (one cache line each), or all of them when they are more than the
    // cache lines of the array (Eclip_1480 at 5 ppm: ~32,000 of 3.8 million bins per MS1 spectrum, ~2,300 of 2.6 million per MS2 spectrum)
    auto& scratch = mass_bin_scratch_;
    if (scratch.written_mass_bins.size() < binned_log_masses_.size() / 16)
    {
      for (Size bin : scratch.written_mass_bins) { scratch.mass_intensities[bin] = .0f; }
    }
    else { std::fill(scratch.mass_intensities.begin(), scratch.mass_intensities.begin() + binned_log_masses_.size(), .0f); }
    scratch.written_mass_bins.clear();
    scratch.mass_intensities_clean = true;
    // a decoy instance (report_FDR) does not keep the state, which would be a second copy next to the one of its target instance
    if (target_decoy_type_ != PeakGroup::target) { mass_bin_scratch_ = MassBinScratch_ {}; }
    return mass_bin_charge_ranges;
  }

  // With the mass bins selected by updateMassBins_, select peaks from the same mass in the original input spectrum
  void SpectralDeconvolution::getCandidatePeakGroups_(const std::vector<MassBinChargeRange_>& mass_bin_charge_ranges)
  {
    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];
    double tol = tolerance_[ms_level_ - 1];
    int charge_range = current_max_charge_;
    Size mass_bin_size = binned_log_masses_.size();
    int log_mz_peak_size = (int)log_mz_peaks_.size();
    // this stores which peak is now being considered per charge. Per charge, peak is considered from left (lowest m/z) to right (highest m/z).
    auto current_peak_index = std::vector<int>(charge_range, 0);
    // the candidates go into the vector of the previous spectrum's candidates (pooled, empty), for its capacity
    candidate_pool_.clear();
    deconvolved_spectrum_.swapPeakGroups(candidate_pool_);
    candidate_pool_.clear(); // what the spectrum held before (nothing: generatePeakGroupsFromSpectrum_ cleared it) is dropped
    deconvolved_spectrum_.reserve(mass_bin_charge_ranges.size());
    Size range_index = 0; // the selected mass bin in mass_bin_charge_ranges
    auto& peak_bin_numbers = peak_bin_numbers_pool_; // pooled: every element is set below
    peak_bin_numbers.resize(log_mz_peak_size);
    // per peak, store the m/z bin number for fast processing
    for (int i = 0; i < log_mz_peak_size; i++)
    {
      peak_bin_numbers[i] = getBinNumber_(log_mz_peaks_[i].logMz, mz_bin_min_value_, bin_mul_factor);
    }

    // main iteration. mass_bin_charge_ranges gives the range of charges for each mass bin
    while (range_index < mass_bin_charge_ranges.size())
    {
      const MassBinChargeRange_& charge_range_of_bin = mass_bin_charge_ranges[range_index];
      const Size mass_bin_index = charge_range_of_bin.bin;
      double log_m = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_mul_factor);
      double mass = exp(log_m);

      // Peaks for this mass bin are staged in a per-worker scratch buffer; a PeakGroup with an
      // exact-size peak vector is only created when the bin yields a candidate.
      static thread_local std::vector<LogMzPeak> bin_peaks;
      bin_peaks.clear();
      // the range of isotope span. For a given peak the peaks within the span are searched.
      Size right_index = activeAveragine_().getRightCountFromApex(mass);
      Size left_index = activeAveragine_().getLeftCountFromApex(mass);

      // scan through charge - from mass to m/z
      for (size_t j = charge_range_of_bin.min_charge_index; j <= (size_t)charge_range_of_bin.max_charge_index; j++)
      {
        int max_peak_index = -1;
        size_t abs_charge = j + 1;
        int bin_offset = binned_universal_pattern_[j];

        if (mass_bin_index < (size_t)bin_offset) { continue; }

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
        if (max_peak_index < log_mz_peak_size - 1 && peak_bin_numbers[max_peak_index + 1] == b_index + 1
            && log_mz_peaks_[max_peak_index + 1].intensity > max_intensity)
        {
          continue;
        }

        // now we have a matching peak for this mass of charge  abs_charge. From here, isotope peaks are collected
        const double mz = log_mz_peaks_[max_peak_index].mz; // charged mz
        const double iso_delta = iso_da_distance_ / (double)abs_charge;
        const double inv_iso_delta = 1.0 / iso_delta;
        double mz_delta = std::min(max_mass_dalton_tolerance_ / (double)abs_charge, tol * mz); //

        double max_mz = mz;

        for (int peak_index = max_peak_index; peak_index < log_mz_peak_size; peak_index++)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          double mz_diff = observed_mz - mz;
          int tmp_i;
          if (! roundedQuotient(mz_diff, inv_iso_delta, tmp_i)) { tmp_i = (int)round(mz_diff / iso_delta); }

          if (observed_mz - max_mz > (double)right_index * iso_delta + mz_delta) { break; }

          if (std::abs(mz_diff - tmp_i * iso_delta) < mz_delta) // if peak is signal
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size && ! (bin < excluded_mass_bins_for_decoy_runs_.size() && excluded_mass_bins_for_decoy_runs_[bin]))
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = (int)abs_charge;
              p.isotopeIndex = tmp_i;
              bin_peaks.push_back(p);
            }
          }
        }

        for (int peak_index = max_peak_index - 1; peak_index >= 0; peak_index--)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          double mz_diff = observed_mz - mz;
          int tmp_i;
          if (! roundedQuotient(mz_diff, inv_iso_delta, tmp_i)) { tmp_i = (int)round(mz_diff / iso_delta); }

          if (max_mz - observed_mz > (float)left_index * iso_delta + mz_delta) { break; }
          if (std::abs(mz_diff - tmp_i * iso_delta) < mz_delta)
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size && ! (bin < excluded_mass_bins_for_decoy_runs_.size() && excluded_mass_bins_for_decoy_runs_[bin]))
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = (int)abs_charge;
              p.isotopeIndex = tmp_i;
              bin_peaks.push_back(p);
            }
          }
        }
      }

      if (! bin_peaks.empty())
      {
        double max_intensity = -1.0;
        double t_mass = .0;
        static thread_local std::vector<LogMzPeak> new_peaks;
        new_peaks.clear();
        for (const auto& p : bin_peaks)
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

        int apex_index = (int)activeAveragine_().getApexIndex(t_mass);
        for (auto& p : bin_peaks)
        {
          p.isotopeIndex = (int)round((p.getUnchargedMass() - t_mass) / iso_da_distance_);

          if (std::abs(t_mass - p.getUnchargedMass() + iso_da_distance_ * p.isotopeIndex) > iso_tolerance) { continue; }
          p.isotopeIndex += apex_index;
          new_peaks.push_back(p);
          min_off = min_off > p.isotopeIndex ? p.isotopeIndex : min_off;
          max_off = max_off < p.isotopeIndex ? p.isotopeIndex : max_off;
          max_charge = max_charge < p.abs_charge ? p.abs_charge : max_charge;
        }

        if (min_off != max_off)
        {
          PeakGroup pg(1, charge_range_of_bin.max_charge_index + 1, // make an empty peakGroup (mass)
                       is_positive_);
          pg.setIsotopeDaDistance(iso_da_distance_);
          std::vector<LogMzPeak> exact_peaks(new_peaks.begin(), new_peaks.end());
          pg.swap(exact_peaks);
          pg.updateMonoMassAndIsotopeIntensities(tol);

          if (pg.getMonoMass() < current_min_mass_ || pg.getMonoMass() > current_max_mass_) { continue; }
          pg.setScanNumber(deconvolved_spectrum_.getScanNumber());
          pg.sort();

          // fused: scored now (as scoreAndFilterPeakGroups_ would, in candidate order), and only kept if it passes
          if (! fused_scoring_ || scorePeakGroup_(pg)) { deconvolved_spectrum_.push_back(std::move(pg)); }
        }
      }
      range_index++;
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
    if (target_decoy_type == PeakGroup::signal_decoy)
      allowed_iso_error_ = param_.getValue("allowed_isotope_error");
    target_dspec_for_decoy_calculation_ = &target_dspec_for_decoy_calcualtion;
  }

  // spectral deconvolution main function
  void SpectralDeconvolution::generatePeakGroupsFromSpectrum_()
  {
    deconvolved_spectrum_.clear();
    int current_charge_range = current_max_charge_;
    int tmp_peak_cntr = current_charge_range - min_support_peak_count_;

    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;
    double mass_bin_max_value = std::min(log_mz_peaks_.back().logMz - universal_pattern_[tmp_peak_cntr],
                                         log(current_max_mass_ + (double)activeAveragine_().getRightCountFromApex(current_max_mass_) + 1.0));

    double bin_mul_factor = bin_mul_factors_[ms_level_ - 1];

    mass_bin_min_value_ = log(std::max(1.0, 50 - activeAveragine_().getAverageMassDelta(50)));
    mz_bin_min_value_ = log_mz_peaks_[0].logMz;

    double mz_bin_max_value = log_mz_peaks_.back().logMz;
    Size mass_bin_number = getBinNumber_(mass_bin_max_value, mass_bin_min_value_, bin_mul_factor) + 1;
    binned_universal_pattern_.clear();

    for (int i = 0; i < current_charge_range; i++)
    {
      binned_universal_pattern_.push_back((int)round((mz_bin_min_value_ - universal_pattern_[i] - mass_bin_min_value_) * bin_mul_factor));
    }

    binned_harmonic_patterns = Matrix<int>(harmonic_pattern_matrix_.rows(), harmonic_pattern_matrix_.cols(), 0);

    for (Size k = 0; k < (Size)harmonic_pattern_matrix_.rows(); k++)
    {
      for (Size i = 0; i < (Size)harmonic_pattern_matrix_.cols(); i++)
      {
        binned_harmonic_patterns.setValue(
          k, i, (int)round((mz_bin_min_value_ - harmonic_pattern_matrix_.getValue(k, i) - mass_bin_min_value_) * bin_mul_factor));
      }
    }

    Size mz_bin_number = getBinNumber_(mz_bin_max_value, mz_bin_min_value_, bin_mul_factor) + 1;
    auto& binned_log_mz_peak_intensities = mz_intensities_pool_; // pooled: refilled here
    binned_log_mz_peak_intensities.assign(mz_bin_number, .0f);

    // bin log mzs
    binLogMzPeaks_(mz_bin_number, binned_log_mz_peak_intensities);
    binned_log_masses_.resize(mass_bin_number);
    binned_log_masses_.reset();

    // for FDR estimation
    if (! excluded_peak_masses_for_decoy_runs_.empty())
    {
      excluded_mass_bins_for_decoy_runs_.resize(binned_log_masses_.size());
      excluded_mass_bins_for_decoy_runs_.reset();
      // always positive
      int bin_offset = (int)round(tol_div_factor);
      for (double m : excluded_peak_masses_for_decoy_runs_)
      {
        if (m <= 0) { continue; }
        Size j = getBinNumber_(log(m), mass_bin_min_value_, bin_mul_factors_[ms_level_ - 1]);
        if ((int)j >= bin_offset && j < excluded_mass_bins_for_decoy_runs_.size() - bin_offset)
        {
          for (int k = -bin_offset; k <= bin_offset; k++)
            excluded_mass_bins_for_decoy_runs_[j + k] = true;
        }
      }
    }

    // for targeted deconvolution
    if (! target_mono_masses_.empty())
    {
      target_mass_bins_.resize(binned_log_masses_.size());
      target_mass_bins_.reset();
      for (double& tm : target_mono_masses_)
      {
        for (int off = -1; off < 2; off++)
        {
          double m = tm + off * iso_da_distance_;
          double mass_delta = activeAveragine_().getMostAbundantMassDelta(m);

          Size j = getBinNumber_(log(m + mass_delta), mass_bin_min_value_, bin_mul_factors_[ms_level_ - 1]);
          if (j < 1) { continue; }

          if (j >= target_mass_bins_.size() - 2) { break; }

          target_mass_bins_[j - 1] = true;
          target_mass_bins_[j] = true;
          target_mass_bins_[j + 1] = true;
        }
      }
    }

    // main algorithm to generate mass bins
    const auto mass_bin_charge_ranges = updateMassBins_(binned_log_mz_peak_intensities);

    // detection limit of the spectrum: the intensity of its weakest peak (after the peak filter). Noise below it cannot be observed;
    // it floors the per-charge noise power of the peak groups (PeakGroup::updateChargeRange_)
    detection_limit_ = 0;
    for (const auto& p : deconvolved_spectrum_.getOriginalSpectrum())
    {
      if (p.getIntensity() <= 0) { continue; }
      if (detection_limit_ == 0 || p.getIntensity() < detection_limit_) { detection_limit_ = p.getIntensity(); }
    }
    // Inside an active parallel region (the spectrum workers of FLASHDeconvAlgorithm, or the two decoy runs of a spectrum) the candidates
    // are scored serially in candidate order (scoreAndFilterPeakGroups_). There each candidate is scored as soon as it is built and kept
    // only if it passes, so that a worker never holds all candidates of a spectrum: the scoring of a candidate reads only the candidate and
    // state fixed before the candidates are built, so the kept groups, their scores and their order are the same.
    fused_scoring_ = false;
#ifdef _OPENMP
    fused_scoring_ = omp_in_parallel();
#endif
    // main algorithm to generate peak groups
    getCandidatePeakGroups_(mass_bin_charge_ranges);
    scoreAndFilterPeakGroups_();
  }

  bool SpectralDeconvolution::isPeakGroupInExcludedMassForDecoyRuns_(const PeakGroup& peak_group, double tol, int offset) const
  {
    if (excluded_masses_for_decoy_runs_.empty()) return false;
    bool match = false;
    double mass = offset * iso_da_distance_ + peak_group.getMonoMass();
    double mass_da_tol = mass * tol;
    // the excluded masses at or below mass + tol, downwards while within tol (a mass at or above the heaviest excluded one was not checked:
    // the search started at end())
    auto it = std::upper_bound(excluded_masses_for_decoy_runs_.begin(), excluded_masses_for_decoy_runs_.end(), mass + mass_da_tol);
    while (it != excluded_masses_for_decoy_runs_.begin())
    {
      --it;
      double delta = mass - *it;
      if (delta > mass_da_tol) break;
      if (std::abs(delta) < mass_da_tol)
      {
        match = true;
        break;
      }
    }
    return match;
  }


  bool SpectralDeconvolution::scorePeakGroup_(PeakGroup& peak_group)
  {
    const double tol = tolerance_[ms_level_ - 1];
    const double snr_threshold = min_snr_[ms_level_ - 1];
    int offset = 0;
    peak_group.setTargetDecoyType(target_decoy_type_);
    peak_group.setDetectionLimit(detection_limit_);

    if (isPeakGroupInExcludedMassForDecoyRuns_(peak_group, tol)) { return false; }

    static thread_local std::vector<LogMzPeak> noisy_peaks; // per-worker scratch
    int num_iteration = 10;
    bool mass_determined = false;
    double prev_cos = 0;
    // isotope pattern matching and qscore update part. Isotope pattern matching is done multiple times until finding the maximum isotope cosine
    for (int k = 0; k < num_iteration; k++)
    {
      peak_group.recruitAllPeaksInSpectrum(deconvolved_spectrum_.getOriginalSpectrum(), tol, activeAveragine_(),
                                           peak_group.getMonoMass() + offset * iso_da_distance_, noisy_peaks, true, &original_log_mzs_,
                                           &original_nearest_peaks_);
      if (isPeakGroupInExcludedMassForDecoyRuns_(peak_group, tol)) { break; }
      // min cosine is checked in here. mono mass is also updated one last time. SNR, per charge SNR, and avg errors are updated here.
      const auto& [z1, z2] = peak_group.getAbsChargeRange();
      offset = peak_group.updateQscore(noisy_peaks, activeAveragine_(), min_isotope_cosine_[ms_level_ - 1], tol,
                                       (z1 + z2) < 2 * low_charge_, excluded_masses_for_decoy_runs_, false);

      if (prev_cos > peak_group.getIsotopeCosine() || offset == 0 || k >= num_iteration - 1) //
      {
        if (peak_group.getChargeSNR(peak_group.getRepAbsCharge()) < snr_threshold) // to speed up
          break;
        if (offset != 0) peak_group.recruitAllPeaksInSpectrum(deconvolved_spectrum_.getOriginalSpectrum(), tol, activeAveragine_(),
                                                              peak_group.getMonoMass() + offset * iso_da_distance_, noisy_peaks, true, &original_log_mzs_,
                                                              &original_nearest_peaks_);

        peak_group.updateQscore(noisy_peaks, activeAveragine_(), min_isotope_cosine_[ms_level_ - 1], tol,
                                (z1 + z2) < 2 * low_charge_, excluded_masses_for_decoy_runs_, true);
        mass_determined = true;
        break;
      }
      prev_cos = peak_group.getIsotopeCosine();
    }

    if (isPeakGroupInExcludedMassForDecoyRuns_(peak_group, tol))
    {
      return false;
    }

    if (! mass_determined || peak_group.empty() || peak_group.getQscore() <= 0 || peak_group.getMonoMass() < current_min_mass_
        || peak_group.getMonoMass() > current_max_mass_)
      return false;

    auto [z1, z2] = peak_group.getAbsChargeRange();

    // a high-charge mass must span several charge states in MS1; fragment ions in MSn mostly appear at one charge
    if (ms_level_ == 1 && z1 > low_charge_ && (z2 - z1) < min_support_peak_count_) { return false; }

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
    return true;
  }

  void SpectralDeconvolution::scoreAndFilterPeakGroups_()
  {
    double tol = tolerance_[ms_level_ - 1];
    const int num_peak_groups = static_cast<int>(deconvolved_spectrum_.size());
    auto selected = boost::dynamic_bitset<>(num_peak_groups);
    double snr_threshold = min_snr_[ms_level_ - 1];
    // detection limit of the spectrum (generatePeakGroupsFromSpectrum_)
    const float detection_limit = detection_limit_;

    // The peak groups are scored independently of each other, in parallel, except inside an active parallel region (the spectrum workers
    // of FLASHDeconvAlgorithm, which deconvolve whole spectra in parallel): there the loops run serially as plain loops, so that an
    // exception stays within the worker's own try/catch. The result is the same either way.
    bool parallel_peak_groups = true;
#ifdef _OPENMP
    parallel_peak_groups = ! omp_in_parallel();
#endif
    // whether peak group i is kept
    const auto score_peak_group = [&](const int i) -> bool { return scorePeakGroup_(deconvolved_spectrum_[i]); };
    if (fused_scoring_)
    {
      // the candidates were scored as they were built (getCandidatePeakGroups_), and only those that passed were kept
      selected.set();
      fused_scoring_ = false;
    }
    else if (parallel_peak_groups && num_peak_groups > 100) // a small spectrum (e.g. FLASHIda's, in real time) is scored serially
    {
  #pragma omp parallel for default(none) shared(selected, score_peak_group)
      for (int i = 0; i < (int)deconvolved_spectrum_.size(); i++)
      {
        if (! score_peak_group(i)) { continue; }
  #pragma omp critical
        selected[i] = true;
      }
    }
    else
    {
      for (int i = 0; i < (int)deconvolved_spectrum_.size(); i++)
      {
        if (score_peak_group(i)) { selected[i] = true; }
      }
    }

    Size selected_count = selected.count();
    if (selected_count == 0)
    {
      // no peak groups: the candidates are destroyed now, their vector (capacity) goes back to the pool
      deconvolved_spectrum_.swapPeakGroups(candidate_pool_);
      candidate_pool_.clear();
      return;
    }

    std::vector<PeakGroup> filtered_peak_groups;
    filtered_peak_groups.reserve(selected_count);

    // Candidates from different charge ranges can converge to the same monoisotopic mass, and masses move by isotope offsets during
    // scoring, so the groups are not sorted by mass here. Sort them, and of the groups with the same mass (within tolerance) keep the
    // one explaining the most signal: the first one in candidate order may cover only part of the ion's charge states.
    std::vector<Size> selected_indices;
    selected_indices.reserve(selected_count);
    for (Size index = selected.find_first(); index != selected.npos; index = selected.find_next(index))
    {
      selected_indices.push_back(index);
    }
    std::stable_sort(selected_indices.begin(), selected_indices.end(),
                     [this](Size a, Size b) { return deconvolved_spectrum_[a].getMonoMass() < deconvolved_spectrum_[b].getMonoMass(); });
    for (Size index : selected_indices) // each candidate is visited once, so it can be moved
    {
      auto& peak_group = deconvolved_spectrum_[index];
      if (! filtered_peak_groups.empty() && peak_group.getMonoMass() - filtered_peak_groups.back().getMonoMass() <= peak_group.getMonoMass() * tol)
      {
        if (peak_group.getIntensity() > filtered_peak_groups.back().getIntensity()) { filtered_peak_groups.back() = std::move(peak_group); }
        continue;
      }
      filtered_peak_groups.push_back(std::move(peak_group));
    }

    if (target_decoy_type_ == PeakGroup::noise_decoy)
    {
      for (const auto& pg2 : *target_dspec_for_decoy_calculation_)
      {
        filtered_peak_groups.push_back(pg2);
      }
    }

    deconvolved_spectrum_.swapPeakGroups(filtered_peak_groups);
    // the remaining (moved-from or unselected) candidates are destroyed now, their vector (capacity) goes back to the pool
    filtered_peak_groups.clear();
    candidate_pool_.swap(filtered_peak_groups);
    deconvolved_spectrum_.sort();

    removeChargeErrorPeakGroups_(deconvolved_spectrum_, target_decoy_type_);

    // final harmonics removal
    const int num_harmonics = static_cast<int>(deconvolved_spectrum_.size());
    selected = boost::dynamic_bitset<>(num_harmonics);
    filtered_peak_groups = std::vector<PeakGroup>();

    // whether peak group i is no harmonic mass
    const auto passes_harmonics = [&](const int i) -> bool {
      const auto& peak_group = deconvolved_spectrum_[i];
      bool pass = true;
      const auto& [z1, z2] = peak_group.getAbsChargeRange();

      for (auto hz : harmonic_charges_)
      {
        PeakGroup pg(std::min(current_max_charge_, z1 * hz), std::min(current_max_charge_, z2 * hz), is_positive_);
        pg.setIsotopeDaDistance(iso_da_distance_);
        pg.setDetectionLimit(detection_limit);
        const auto& [z1_, z2_] = pg.getAbsChargeRange();
        if (z2 - z1 > z2_ - z1_) break; // if harmonic charges are too high stop
        pg.setMonoisotopicMass(peak_group.getMonoMass() * hz);
        static thread_local std::vector<LogMzPeak> nps; // per-worker scratch
        pg.recruitAllPeaksInSpectrum(deconvolved_spectrum_.getOriginalSpectrum(), tol, activeAveragine_(), pg.getMonoMass(), nps, true, &original_log_mzs_,
                                     &original_nearest_peaks_);
        pg.updateQscore(nps, activeAveragine_(), min_isotope_cosine_[ms_level_ - 1], tol,
                        (z1 + z2) * hz < 2 * low_charge_, excluded_masses_for_decoy_runs_, true);

        if (pg.getQscore() > 0 && pg.getSNR() > peak_group.getSNR())
        {
          pass = false;
          break;
        }
      }
      return pass;
    };
    if (parallel_peak_groups && num_harmonics > 100)
    {
  #pragma omp parallel for default(none) shared(selected, passes_harmonics)
      for (int i = 0; i < (int)deconvolved_spectrum_.size(); i++)
      {
        if (! passes_harmonics(i)) { continue; }
  #pragma omp critical
        selected[i] = true;
      }
    }
    else
    {
      for (int i = 0; i < (int)deconvolved_spectrum_.size(); i++)
      {
        if (passes_harmonics(i)) { selected[i] = true; }
      }
    }

    selected_count = selected.count();
    if (selected_count == 0)
    {
      deconvolved_spectrum_.clear();
      return;
    }

    filtered_peak_groups.reserve(selected_count);
    Size index = selected.find_first();

    while (index != selected.npos)
    {
      auto& peak_group = deconvolved_spectrum_[index];
      if (filtered_peak_groups.empty() || filtered_peak_groups.back().getMonoMass() < peak_group.getMonoMass())
        filtered_peak_groups.push_back(std::move(peak_group));

      index = selected.find_next(index);
    }

    deconvolved_spectrum_.setPeakGroups(std::move(filtered_peak_groups));
    deconvolved_spectrum_.sort();

    removeOverlappingPeakGroups_(deconvolved_spectrum_, tol * 1.2, target_decoy_type_);

    filtered_peak_groups.clear();
    filtered_peak_groups.reserve(deconvolved_spectrum_.size());

    std::set<Size> indices; // the signal decoys that pass the checks below (none may pass)
    const bool signal_decoy_run = target_decoy_type_ == PeakGroup::signal_decoy;
    if (signal_decoy_run)
    {
      excluded_masses_for_decoy_runs_.clear();
      for (const auto& pg : *target_dspec_for_decoy_calculation_) // pg are the target peak groups from normal spectrum deconvolution
      {
        for (int j = -allowed_iso_error_; j <= allowed_iso_error_; j++)
        {
          excluded_masses_for_decoy_runs_.push_back(pg.getMonoMass() + j * iso_da_distance_);
        }
      }
      // isPeakGroupInExcludedMassForDecoyRuns_ binary-searches this list; with allowed_isotope_error > 0 the entries interleave
      std::sort(excluded_masses_for_decoy_runs_.begin(), excluded_masses_for_decoy_runs_.end());

      for (Size k = 0; k < deconvolved_spectrum_.size(); k++)
      {
        auto& pg = deconvolved_spectrum_[k];
        if (isPeakGroupInExcludedMassForDecoyRuns_(pg, tol)) continue;
        bool pass = true;
        for (const auto & pg2 : *target_dspec_for_decoy_calculation_)
        {
          if (std::abs(pg.getMonoMass() - pg2.getMonoMass()) < (3 + allowed_iso_error_) * iso_da_distance_ + .1) // if they are close enough
          {
            if (std::abs(pg2.getIsotopeCosine() - pg.getIsotopeCosine()) > .005 * (allowed_iso_error_ + 1))
            {
              pass = false;
              break;
            }
          }
        }
        if (!pass) continue;// || (is_isotope_error && qs == pg.getQscore())
        indices.insert(k);
      }
    }

    for (Size k = 0; k < deconvolved_spectrum_.size(); k++)
    {
      auto& pg = deconvolved_spectrum_[k];

      // a signal decoy must have passed the checks above (when none did, every one was kept)
      if (signal_decoy_run && indices.find(k) == indices.end())
      {
        continue;
      }
      if (pg.getQscore() >= min_qscore_ && pg.getChargeSNR(pg.getRepAbsCharge()) > snr_threshold)
      {
        filtered_peak_groups.push_back(std::move(pg));
      }
    }

    deconvolved_spectrum_.setPeakGroups(std::move(filtered_peak_groups));
    removeExcludedMasses_(deconvolved_spectrum_, excluded_masses_, tol);
    deconvolved_spectrum_.sort();
  }

  float SpectralDeconvolution::getIsotopeCosineAndIsoOffset(double mono_mass,
                                                            const std::vector<float>& per_isotope_intensities,
                                                            int& offset,
                                                            const PrecalculatedAveragine& avg,
                                                            const int iso_int_shift,
                                                            const int window_width,
                                                            const std::vector<double>& excluded_masses)
  {
    offset = 0;
    if ((int)per_isotope_intensities.size() < min_iso_size + iso_int_shift) { return .0; }
    const auto& iso = avg.getReference_(mono_mass);

    int right = (int)avg.getApexIndex(mono_mass) / 4 + 1;
    int left = right;

    right += iso_int_shift;
    left -= iso_int_shift;
    float max_cos = -1000;
    int max_isotope_index = (int)per_isotope_intensities.size(); // exclusive
    int min_isotope_index = -1;                                  // inclusive

    for (int i = 0; i < max_isotope_index; i++)
    {
      if (per_isotope_intensities[i] <= 0) { continue; }

      min_isotope_index = i; // the first one
      break;
    }
    if (max_isotope_index - min_isotope_index < min_iso_size) { return .0; }

    static thread_local std::vector<std::pair<int, float>> offset_cos; // per-worker scratch
    offset_cos.clear();
    offset_cos.reserve(right + left + 1);
    // getCosine(per_isotope_intensities, min_isotope_index, max_isotope_index, iso, offset, min_iso_size) for each offset, with the norm of
    // the observed intensities (the same for every offset, summed in the same order) computed once
    const int a_start = std::max(0, min_isotope_index);
    const int a_end = std::min((int)per_isotope_intensities.size(), max_isotope_index);
    double a_norm = .0;
    for (int j = a_start; j < a_end; j++)
    {
      a_norm += double(per_isotope_intensities[j]) * per_isotope_intensities[j];
    }
    const int iso_size = (int)iso.size(); // read once (IsotopeDistribution::size is not inline)
    // the norm of the averagine pattern, as in getCosine: a trimmed pattern need not have norm one (#10543)
    double b_norm = .0;
    for (const auto& peak : iso) b_norm += double(peak.getIntensity()) * peak.getIntensity();
    const auto cosine = [&](const int iso_offset) -> float {
      double n = .0;
      if (a_end - a_start < min_iso_size) { return 0; }
      for (int j = a_start; j < a_end; j++)
      {
        int i = j - iso_offset;
        if (i >= iso_size || i < 0 || iso[i].getIntensity() <= 0) { continue; }
        else { n += double(per_isotope_intensities[j]) * iso[i].getIntensity(); }
      }
      if (a_norm <= 0 || b_norm <= 0) { return 0; }
      return float(std::clamp(n / sqrt(a_norm * b_norm), 0.0, 1.0));
    };

    for (int tmp_offset = -left; tmp_offset <= right; tmp_offset++)
    {
      if (window_width >= 0 && abs(tmp_offset - iso_int_shift) > window_width)
        continue;
      if (!excluded_masses.empty())
      {
        bool exclude = false;
        for (auto em : excluded_masses)
        {
          if ( std::abs(mono_mass + (tmp_offset - iso_int_shift)* Constants::ISOTOPE_MASSDIFF_55K_U - em) < mono_mass * 1e-5) // tmp
          {
            exclude = true;
            break;
          }
        }
        if (exclude) continue;
      }
      float tmp_cos = cosine(tmp_offset);
      offset_cos.emplace_back(tmp_offset, tmp_cos);
    }
    if (offset_cos.empty()) return max_cos;

    // Only the first pair after the sort is used (all pairs pass the checks below). A unique largest cosine without NaN is first in any sort
    // order: it is moved to the front. Otherwise (equal largest cosines, whose order the sort decides, or NaN) the sort runs as before.
    Size max_pos = 0;
    Size max_count = 1;
    bool has_nan = std::isnan(offset_cos[0].second);
    for (Size k = 1; k < offset_cos.size(); k++)
    {
      const float c = offset_cos[k].second;
      has_nan = has_nan || std::isnan(c);
      if (c > offset_cos[max_pos].second)
      {
        max_pos = k;
        max_count = 1;
      }
      else if (c == offset_cos[max_pos].second) { max_count++; }
    }
    if (! has_nan && max_count == 1) { std::swap(offset_cos[0], offset_cos[max_pos]); }
    else
    {
      std::sort(offset_cos.begin(), offset_cos.end(),
                [](const std::pair<int, float>& p1, const std::pair<int, float>& p2) { return p1.second > p2.second; });
    }

    for (const auto& [o, c] : offset_cos)
    {
      if (o > right || o < -left) continue;
      if (window_width >= 0 && abs(o - iso_int_shift) > window_width) //
        continue;

      offset = o;
      max_cos = c;
      break;
    }

    max_cos = std::max(max_cos, .0f);
    offset -= iso_int_shift;

    return max_cos;
  }

  float SpectralDeconvolution::getCosine(const std::vector<float>& a, int a_start, int a_end, const IsotopeDistribution& b, int offset, int min_iso_len)
  {
    // Trimming and floating-point normalization can leave both RNA and peptide
    // averagine patterns with a norm that is not exactly one.
    // A cosine above one produces negative noise power in PeakGroup::updateSNR_.
    double n = 0.0, a_norm = 0.0, b_norm = 0.0;
    for (const auto& peak : b) b_norm += double(peak.getIntensity()) * peak.getIntensity();
    a_start = std::max(0, a_start);
    a_end = std::min((int)a.size(), a_end);

    if (a_end - a_start < min_iso_len) { return 0; }

    float max_intensity = 0;
    const int b_size = (int)b.size(); // read once (IsotopeDistribution::size is not inline)

    for (int j = a_start; j < a_end; j++)
    {
      int i = j - offset;
      a_norm += double(a[j]) * a[j];

      if (max_intensity < a[j]) { max_intensity = a[j]; }

      if (i >= b_size || i < 0 || b[i].getIntensity() <= 0) { continue; }
      else { n += double(a[j]) * b[i].getIntensity(); }
    }

    if (a_norm <= 0 || b_norm <= 0) { return 0; }

    return float(std::clamp(n / sqrt(a_norm * b_norm), 0.0, 1.0));
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

    // A peak shared by groups of different charges counts as overlapped for a group only by the part that the winning group's own
    // isotope envelope explains: its averagine pattern scaled on its peaks at that charge that no other group claims (without such
    // peaks the whole peak is its). A low-charge fragment mixed with a weaker high-charge envelope on the same peaks then keeps them.
    std::map<std::pair<int, int>, double> envelope_scales; // (group index, charge) -> scale of the averagine pattern
    auto explained_intensity = [&](int j, int z, double pmz, float pint) -> float {
      const auto& iso = activeAveragine_().getReference_(dspec[j].getMonoMass());
      auto it = envelope_scales.find({j, z});
      if (it == envelope_scales.end())
      {
        double num = 0, den = 0;
        for (const auto& p : dspec[j])
        {
          if (p.abs_charge != z || p.isotopeIndex < 0 || p.isotopeIndex >= (int)iso.size()) continue;
          auto shared = peak_to_pgs.find(p.mz);
          if (shared != peak_to_pgs.end() && shared->second.size() > 1) continue;
          const double e = iso[p.isotopeIndex].getIntensity();
          num += e * p.intensity;
          den += e * e;
        }
        it = envelope_scales.emplace(std::make_pair(j, z), den > 0 ? num / den : -1.0).first;
      }
      if (it->second < 0) return pint;
      for (const auto& p : dspec[j])
      {
        if (p.abs_charge != z || p.mz != pmz) continue;
        if (p.isotopeIndex < 0 || p.isotopeIndex >= (int)iso.size()) return pint;
        return std::min(pint, (float)(it->second * iso[p.isotopeIndex].getIntensity()));
      }
      return pint;
    };

    for (const auto& e : peak_to_pgs)
    {
      const auto& pg_is = e.second;
      double pmz = e.first;
      float pint = mz_to_intensities[pmz];

      if (pg_is.size() == 1) { continue; }

      for (auto i : pg_is)
      {
        bool is_overlap = false;
        float overlap_part = pint;
        double mass1 = dspec[i].getMonoMass();
        int repz1 = (int)round(mass1 / (pmz - FLASHHelperClasses::getChargeMass(is_positive_)));
        const double high_snr = 10;
        if (dspec[i].getSNR() > high_snr && dspec[i].getChargeSNR(repz1) > high_snr) // if signal is almost perfect, skip the filtration
          continue;

        for (auto j : pg_is)
        {
          if (i == j) { continue; }
          double mass2 = dspec[j].getMonoMass();
          int repz2 = (int)round(mass2 / (pmz - FLASHHelperClasses::getChargeMass(is_positive_)));
          if (repz1 == repz2) continue;
          if (dspec[i].getChargeSNR(repz1) > dspec[j].getChargeSNR(repz2)) { continue; }

          is_overlap = true;
          overlap_part = explained_intensity((int)j, repz2, pmz, pint);
          break;
        }
        if (is_overlap) overlap_intensity[i] += overlap_part;
      }
    }

    for (Size i = 0; i < dspec.size(); i++)
    {
      if (dspec[i].getTargetDecoyType() != target_decoy_type) { continue; }
      if ((ms_level_ == 1 && dspec[i].getRepAbsCharge() < min_abs_charge_) || dspec[i].getRepAbsCharge() > max_abs_charge_) { continue; }

      const double mul_factor = .5;
      if (! dspec[i].isTargeted() &&                                    // z1 != z2 &&
          overlap_intensity[i] >= dspec[i].getIntensity() * mul_factor) // If the overlapped intensity takes more than mul_factor * total intensity then
                                                                        // it is a peakgroup with a charge error. the smaller, the harsher
      {
        continue;
      }
      filtered_pg_vec.push_back(std::move(dspec[i]));
    }
    dspec.setPeakGroups(std::move(filtered_pg_vec));
  }

  void SpectralDeconvolution::removeOverlappingPeakGroups_(DeconvolvedSpectrum& dspec, double tol, PeakGroup::TargetDecoyType target_decoy_type)
  {
    if (dspec.empty()) { return; }
    // dspec.sort();
    std::vector<PeakGroup> filtered_pg_vec; //
    filtered_pg_vec.reserve(dspec.size());

    double start_mass = dspec[0].getMonoMass();
    float local_max_SNR = -1.0f;
    Size local_max_index = 0;
    Size last_local_max_index = dspec.size();

    for (Size i = 0; i < dspec.size(); i++)
    {
      double mass = dspec[i].getMonoMass();
      if (mass - start_mass > mass * tol)
      {
        if (! dspec[local_max_index].isTargeted()) // targeted ones were already push_backed.
        {
          if (dspec[local_max_index].getTargetDecoyType() == target_decoy_type
              && last_local_max_index != local_max_index)
            filtered_pg_vec.push_back(std::move(dspec[local_max_index]));
        }
        last_local_max_index = local_max_index;
        start_mass = mass;
        local_max_SNR = -1.0f;
      }

      float snr = dspec[i].getSNR();

      if (local_max_SNR < snr)
      {
        local_max_SNR = snr;
        local_max_index = i;
      }
      if (dspec[i].isTargeted())
      {
        if (dspec[i].getTargetDecoyType() == target_decoy_type)
          filtered_pg_vec.push_back(std::move(dspec[i]));
      }
    }

    if (local_max_SNR >= 0)
    {
      if (! dspec[local_max_index].isTargeted()) // targeted ones were already push_backed.
      {
        if (dspec[local_max_index].getTargetDecoyType() == target_decoy_type
            && last_local_max_index != local_max_index)
          filtered_pg_vec.push_back(std::move(dspec[local_max_index]));
      }
    }

    dspec.setPeakGroups(std::move(filtered_pg_vec));
  }

  void SpectralDeconvolution::removeExcludedMasses_(DeconvolvedSpectrum& dspec, std::vector<double> excluded_masses, double tol) const
  {
    if (excluded_masses.empty()) return;
    std::vector<PeakGroup> filtered_pg_vec; //
    filtered_pg_vec.reserve(dspec.size());
    // spec.sort();

    for (auto& peak_group : dspec)
    {
      if (peak_group.getTargetDecoyType() != target_decoy_type_) continue;
      double delta = peak_group.getMonoMass() * tol;
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
      filtered_pg_vec.push_back(std::move(peak_group));
    }
    dspec.setPeakGroups(std::move(filtered_pg_vec));
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
      if (excluded) { end = (int)(activeAveragine_().getApexIndex(m) + activeAveragine_().getRightCountFromApex(m)); }
      for (int j = start; j <= end + 1; j++)
      {
        if (excluded) excluded_masses_.push_back(m + iso_da_distance_ * j);
        else
          target_mono_masses_.push_back(m + iso_da_distance_ * j);
      }
    }
    // the lookups (upper_bound, then a walk down; the target bin loop stops at the first mass past the last bin) need ascending
    // order, which the isotope masses of two targets within 1 Da break
    if (excluded) std::sort(excluded_masses_.begin(), excluded_masses_.end());
    else std::sort(target_mono_masses_.begin(), target_mono_masses_.end());
  }

  void SpectralDeconvolution::setAveragine(const SpectralDeconvolution::PrecalculatedAveragine& avg)
  {
    avg_ = avg; // copy first: avg may be the shared averagine
    shared_avg_.reset();
  }

  void SpectralDeconvolution::setSharedAveragine_(std::shared_ptr<const PrecalculatedAveragine> avg)
  {
    if (! avg) { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The shared averagine must not be null"); }
    shared_avg_ = std::move(avg);
  }
  void SpectralDeconvolution::releasePools()
  {
    std::vector<float>().swap(mz_intensities_pool_);
    mass_bin_scratch_ = MassBinScratch_ {};
    std::vector<double>().swap(exact_mz_positions_pool_);
    std::vector<Size>().swap(mz_bin_index_reverse_pool_);
    std::vector<Size>().swap(peak_bin_numbers_pool_);
    std::vector<PeakGroup>().swap(candidate_pool_);
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