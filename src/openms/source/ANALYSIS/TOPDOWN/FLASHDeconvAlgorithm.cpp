// --------------------------------------------------------------------------
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

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  /// harmonic charge factors that will be considered for harmonic mass reduction.
  inline const std::vector<int> harmonic_charges_ {2, 3, 5, 7};
  /// high and low charges are differently deconvolved. This value determines the (inclusive) threshold for low charge.
  inline const int low_charge_ = 8; // 8 inclusive

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

  // The main function called from outside. precursor_map_for_FLASHIda is used to read FLASHIda information
  void FLASHDeconvAlgorithm::performSpectrumDeconvolution(const MSSpectrum& spec, const std::vector<DeconvolvedSpectrum>& survey_scans, const int scan_number,
                                                          const std::map<int, std::vector<std::vector<float>>>& precursor_map_for_FLASHIda)
  {
    // First prepare for decoy runs.
    iso_da_distance_ = decoy_flag_ != PeakGroup::noise_decoy ? Constants::ISOTOPE_MASSDIFF_55K_U : Constants::ISOTOPE_MASSDIFF_55K_U * sqrt(7.0);
    previously_deconved_mono_masses_for_decoy.clear();
    previously_deconved_mass_bins_for_decoy.clear();
    excluded_integer_mzs_.clear();

    if (decoy_flag_ == PeakGroup::charge_decoy) // charge decoy
    {
      for (auto& pg : targetFD_->deconvolved_spectrum_)
      {
        auto [m, M] = pg.getAbsChargeRange();
        int min_iso = -1, max_iso = 0;
        for (auto& p : pg)
        {
          previously_deconved_mono_masses_for_decoy.push_back(p.getUnchargedMass());
          min_iso = min_iso < 0 ? p.isotopeIndex : std::min(min_iso, p.isotopeIndex);
          max_iso = std::max(max_iso, p.isotopeIndex);
        }
        for (int i = -2 * M; i <= 0; i++)
        {
          previously_deconved_mono_masses_for_decoy.push_back(pg.getMonoMass() + (min_iso + i) * iso_da_distance_);
        }
        for (int i = 0; i <= 2 * M; i++)
        {
          previously_deconved_mono_masses_for_decoy.push_back(pg.getMonoMass() + (max_iso + i) * iso_da_distance_);
        }
        for (int i = -1; i <= 1; i++)
        {
          previously_deconved_mono_masses_for_decoy.push_back(pg.getMonoMass() + i * iso_da_distance_);
        }
      }
    }
    if (decoy_flag_ == PeakGroup::noise_decoy) // noise decoy
    {
      for (auto& pg : targetFD_->deconvolved_spectrum_)
      {
        for (auto& p : pg)
        {
          excluded_integer_mzs_.insert((int)p.mz);
        }
      }
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

    bin_width_.clear();
    tolerance_ = param_.getValue("tol");

    for (double& j : tolerance_)
    {
      j *= 1e-6;
      bin_width_.push_back(.5 / j);
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
    harmonic_filter_matrix_.clear();
    int charge_range = current_max_charge_;
    for (int i = 0; i < charge_range; i++)
    {
      filter_.push_back(-log(i + 1)); //+
    }

    harmonic_filter_matrix_.resize(harmonic_charges_.size(), charge_range);

    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      int hc = harmonic_charges_[k];
      int n = hc / 2;

      for (int i = 0; i < charge_range; i++)
      {
        double a = i > 0 ? exp(-filter_[i - 1]) : 0;
        double b = exp(-filter_[i]);
        harmonic_filter_matrix_.setValue(k, i, -log(b - (b - a) * n / hc));
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
      if (excluded_integer_mzs_.size() > 0 && excluded_integer_mzs_.find((int)peak.getMZ()) != excluded_integer_mzs_.end())
      {
        continue;
      }

      LogMzPeak log_mz_peak(peak, is_positive_);
      log_mz_peaks_.push_back(log_mz_peak);
    }
  }

  // from bin to raw value
  double FLASHDeconvAlgorithm::getBinValue_(const Size bin, const double min_value, const double bin_width)
  {
    return min_value + (double)bin / bin_width;
  }

  // from value to bin number
  Size FLASHDeconvAlgorithm::getBinNumber_(const double value, const double min_value, const double bin_width)
  {
    if (value < min_value)
    {
      return 0;
    }
    return (Size)(((value - min_value) * bin_width) + .5);
  }

  // From log mz to mz bins.
  void FLASHDeconvAlgorithm::updateMzBins_(const Size bin_number, std::vector<float>& mz_bin_intensities)
  {
    mz_bins_ = boost::dynamic_bitset<>(bin_number);
    double bin_width = bin_width_[ms_level_ - 1];
    for (auto& p : log_mz_peaks_)
    {
      Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_width);
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
    auto support_peak_count = std::vector<int>(mass_bins_.size(), 0); // per mass bin how many peaks are present

    // to calculate continuous charges, the previous charge value per mass should be stored
    auto prev_charges = std::vector<int>(mass_bins_.size(), charge_range + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prev_intensities = std::vector<float>(mass_bins_.size(), 1.0f);

    double bin_width = bin_width_[ms_level_ - 1];
    std::vector<float> sub_max_h_intensity(h_charge_size, .0f);

    for (int i = (int)mz_bin_index_reverse.size() - 1; i >= 0; i--)
    {
      mz_bin_index = mz_bin_index_reverse[i];
      float intensity = mz_intensities[mz_bin_index];
      double mz = -1.0, log_mz = 0;

      log_mz = getBinValue_(mz_bin_index, mz_bin_min_value_, bin_width);
      mz = exp(log_mz);
      // scan through charges
      for (int j = 0; j < charge_range; j++)
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

        if (previously_deconved_mono_masses_for_decoy.size() > 0 && previously_deconved_mass_bins_for_decoy[mass_bin_index])
        {
          continue;
        }

        auto& spc = support_peak_count[mass_bin_index];
        int abs_charge = (j + 1);
        float& prev_intensity = prev_intensities[mass_bin_index];
        int& prev_charge = prev_charges[mass_bin_index];
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
            Size next_iso_bin = getBinNumber_(log_mz + diff, mz_bin_min_value_, bin_width);
            if (next_iso_bin > 0 && next_iso_bin < mz_bins_.size() && mz_bins_[next_iso_bin])
            {
              iso_exist = true;
              pass_first_check = true;
            }

            // harmonic check
            if (iso_exist)
            {
              double h_threshold = std::min(intensity, mz_intensities[next_iso_bin]);
              for (size_t k = 0; k < h_charge_size; k++)
              {
                int hc = harmonic_charges_[k];
                if (hc * abs_charge > current_max_charge_)
                {
                  break;
                }
                double hdiff = diff / hc;
                {
                  Size next_harmonic_iso_bin = getBinNumber_(log_mz + hdiff, mz_bin_min_value_, bin_width);
                  float harmonic_intensity = mz_intensities[next_harmonic_iso_bin];
                  if (next_harmonic_iso_bin < mz_bins_.size() && mz_bins_[next_harmonic_iso_bin] && harmonic_intensity > h_threshold / 2.0 // no perfect filtration. Just obvious ones are filtered out
                  )
                  {
                    pass_first_check = false;
                    sub_max_h_intensity[k] = sub_max_h_intensity[k] < harmonic_intensity ? harmonic_intensity : sub_max_h_intensity[k];
                  }
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
            for (size_t k = 0; k < h_charge_size; k++)
            {
              long hmz_bin_index = mass_bin_index - harmonic_bin_offset_matrix_.getValue(k, j);
              if (hmz_bin_index > 0 && hmz_bin_index < (int)mz_bins_.size() && mz_bins_[hmz_bin_index])
              {
                float harmonic_intensity = mz_intensities[hmz_bin_index];
                if (harmonic_intensity > low_threshold && harmonic_intensity < high_threshold)
                {
                  max_h_intensity = std::max(max_h_intensity, harmonic_intensity);
                  is_harmonic = true;
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
              spc--;
            }
          }
          else if (abs_charge <= low_charge_) // for low charge, include the mass if isotope is present
          {
            mass_intensities[mass_bin_index] += intensity + support_peak_intensity;
            spc++;
            {
              mass_bins_[mass_bin_index] = true;
            }
          }
        }
        else if (abs_charge <= low_charge_) // if, for low charge, no isotope peak exists or is harmonic
        {
          mass_intensities[mass_bin_index] -= intensity + max_h_intensity;
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
    double bin_width = bin_width_[ms_level_ - 1];
    Matrix<int> abs_charge_ranges(2, mass_bins_.size(), INT_MAX);
    for (Size i = 0; i < mass_bins_.size(); i++)
    {
      abs_charge_ranges.setValue(1, (int)i, INT_MIN);
    }
    Size mz_bin_index = mz_bins_.find_first();
    long bin_size = (long)mass_bins_.size();

    auto to_skip = mass_bins_.flip();
    mass_bins_.reset();

    while (mz_bin_index != mz_bins_.npos)
    {
      long max_index = -1;
      float max_intensity = -1e11f;
      int max_intensity_abs_charge_range = 0;

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

        if (previously_deconved_mono_masses_for_decoy.size() > 0 && previously_deconved_mass_bins_for_decoy[mass_bin_index])
        {
          continue;
        }

        if (target_mono_masses_.size() > 0 && target_mass_bins_[mass_bin_index])
        {
          float t = mass_intensities[mass_bin_index];

          if (t <= 0)
          { // no signal
            continue;
          }
          max_intensity = 1e38f;
          max_index = mass_bin_index;
          max_intensity_abs_charge_range = j;
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

          int abs_charge = (j + 1);

          if (max_intensity < t)
          {
            bool artifact = false;

            double original_log_mass = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width);
            double mass = exp(original_log_mass);
            {
              double log_mass = log(mass);
              if (log_mass < 1)
              {
                continue;
              }

              for (int h : harmonic_charges_)
              {
                if (h * abs_charge > current_max_charge_)
                {
                  break;
                }
                for (int f : {-1, 1})
                {
                  double hmass = log_mass - log(h) * f;
                  Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_, bin_width);
                  if (hmass_index > 0 && hmass_index < mass_bins_.size() - 1)
                  {
                    if (mass_intensities[hmass_index] >= t)
                    {
                      artifact = true;
                      break;
                    }
                  }
                }
                if (artifact)
                {
                  break;
                }
              }
              // max_intensity_abs_charge off by one here
              if (!artifact)
              {
                for (int coff = 1; coff <= 2 && !artifact; coff++)
                {
                  for (int f : {-1, 1})
                  {
                    if (abs_charge + f * coff <= 0)
                    {
                      continue;
                    }
                    double hmass = log_mass - log(abs_charge) + log(abs_charge + f * coff);
                    Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_, bin_width);
                    if (hmass_index > 0 && hmass_index < mass_bins_.size() - 1)
                    {
                      if (mass_intensities[hmass_index] >= t)
                      {
                        artifact = true;
                        break;
                      }
                    }
                  }
                }
              }
            }

            if (!artifact)
            {
              max_intensity = t;
              max_index = mass_bin_index;
              max_intensity_abs_charge_range = j;
            }
            else
            {
              to_skip[mass_bin_index] = true;
              mass_bins_[mass_bin_index] = false;
            }
          }
        }
      }

      if (max_index >= 0 && max_index < bin_size)
      {
        abs_charge_ranges.setValue(0, max_index, std::min(abs_charge_ranges.getValue(0, max_index), max_intensity_abs_charge_range));
        abs_charge_ranges.setValue(1, max_index, std::max(abs_charge_ranges.getValue(1, max_index), max_intensity_abs_charge_range));
        mass_bins_[max_index] = true;
      }
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }

    return abs_charge_ranges;
  }

  // update mass bins which will be used to select peaks in the input spectrum...
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins_(const std::vector<float>& mz_intensities)
  {
    auto mass_intensities = std::vector<float>(mass_bins_.size(), 0);

    updateCandidateMassBins_(mass_intensities, mz_intensities);

    auto per_mass_abs_charge_ranges = filterMassBins_(mass_intensities);

    return per_mass_abs_charge_ranges;
  }

  // With mass_bins_ from updateMassBins_ function, select peaks from the same mass in the original input spectrum
  void FLASHDeconvAlgorithm::getCandidatePeakGroups_(const Matrix<int>& per_mass_abs_charge_ranges)
  {
    double bin_width = bin_width_[ms_level_ - 1];
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
      peak_bin_numbers[i] = getBinNumber_(log_mz_peaks_[i].logMz, mz_bin_min_value_, bin_width);
    }
    std::vector<double> total_harmonic_intensity(harmonic_charges_.size(), .0);
    std::vector<int> h_prev_iso(harmonic_charges_.size(), 0);
    std::vector<float> h_max_isotope_intensity(harmonic_charges_.size(), .0f);

    // main iteration. per_mass_abs_charge_ranges gives the range of charges for each mass bin
    while (mass_bin_index != mass_bins_.npos)
    {
      double log_m = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width);
      double mass = exp(log_m);

      PeakGroup pg(1, per_mass_abs_charge_ranges.getValue(1, mass_bin_index) + 1, // make an empty peakGroup (mass)
                   is_positive_);

      pg.reserve(charge_range * 128);
      pg.setIsotopeDaDistance(iso_da_distance_);
      // the range of isotope span. For a given peak the peaks within the span are searched.
      Size right_index = avg_.getRightCountFromApex(mass);
      Size left_index = avg_.getLeftCountFromApex(mass);

      double total_signal_intensity = 0;
      std::fill(total_harmonic_intensity.begin(), total_harmonic_intensity.end(), .0);

      // scan through charge - from mass to m/z
      for (size_t j = per_mass_abs_charge_ranges.getValue(0, mass_bin_index); j <= (size_t)per_mass_abs_charge_ranges.getValue(1, mass_bin_index); j++)
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
        const double mz = log_mz_peaks_[max_peak_index].mz;
        const double iso_delta = iso_da_distance_ / (double)abs_charge;
        double mz_delta = std::min(.32 / (double)abs_charge, tol * mz);

        double max_mz = mz;
        float max_peak_intensity = log_mz_peaks_[max_peak_index].intensity;
        float max_isotope_intensity = .0;
        int prev_iso = -1000;
        std::fill(h_prev_iso.begin(),h_prev_iso.end(),0);
        std::fill(h_max_isotope_intensity.begin(),h_max_isotope_intensity.end(),.0);

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
            if (bin < mass_bin_size && !(bin < previously_deconved_mass_bins_for_decoy.size() && previously_deconved_mass_bins_for_decoy[bin]))
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = (int)abs_charge;
              p.isotopeIndex = tmp_i;
              pg.push_back(p);
              if (max_peak_intensity < intensity)
              {
                max_peak_intensity = intensity;
                max_mz = observed_mz;
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
          else // if(abs_charge > low_charge_)
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
              if (tmp_hi / hc >= tmp_i + 1)
              {
                break;
              }

              if (abs(tmp_hi / hc - tmp_i) < 0.32)
              {
                continue;
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
                // total_harmonic_intensity[l] += std::min(max_peak_intensity, intensity) / (hc - 1); //  min(max_peak_intensity, intensity) to minimize effect of outlier peaks
              }
            }
          }
        }

        total_signal_intensity += max_isotope_intensity;
        for (Size l = 0; l < harmonic_charges_.size(); l++)
        {
          total_harmonic_intensity[l] += std::min(max_peak_intensity, h_max_isotope_intensity[l]);
        }

        max_isotope_intensity = .0;
        prev_iso = -1000;
        std::fill(h_prev_iso.begin(),h_prev_iso.end(),0);
        std::fill(h_max_isotope_intensity.begin(),h_max_isotope_intensity.end(),.0);

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
            if (bin < mass_bin_size && !(bin < previously_deconved_mass_bins_for_decoy.size() && previously_deconved_mass_bins_for_decoy[bin]))
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = (int)abs_charge;
              p.isotopeIndex = tmp_i;
              pg.push_back(p);
              // total_signal_intensity += intensity;
              if (max_peak_intensity < intensity)
              {
                max_peak_intensity = intensity;
                max_mz = observed_mz;
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
          else // if(abs_charge > low_charge_)
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
              if (tmp_hi / hc >= tmp_i - 1)
              {
                break;
              }
              if (abs(tmp_hi / hc - tmp_i) < 0.32)
              {
                continue;
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
          total_harmonic_intensity[l] += std::min(max_peak_intensity, h_max_isotope_intensity[l]);
        }
      }

      if (total_signal_intensity / 2.0 > *std::max_element(total_harmonic_intensity.begin(), total_harmonic_intensity.end()) &&
          !pg.empty())
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

        for (auto& p : pg)
        {
          p.isotopeIndex = (int)round((p.getUnchargedMass() - t_mass) / iso_da_distance_);
          if (abs(t_mass - p.getUnchargedMass() + iso_da_distance_ * p.isotopeIndex) > iso_tolerance)
          {
            continue;
          }
          new_peaks.push_back(p);
          min_off = min_off > p.isotopeIndex ? p.isotopeIndex : min_off;
          max_off = max_off < p.isotopeIndex ? p.isotopeIndex : max_off;
          max_charge = max_charge < p.abs_charge ? p.abs_charge : max_charge;
        }
        if (min_off == max_off) {}
        else
        {
          pg.swap(new_peaks);
          for (auto& p : pg)
          {
            p.isotopeIndex -= min_off;
          }
          pg.updateMonomassAndIsotopeIntensities();

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

  void FLASHDeconvAlgorithm::setDecoyFlag(PeakGroup::DecoyFlag flag, FLASHDeconvAlgorithm& targetFD)
  {
    decoy_flag_ = flag;
    targetFD_ = &targetFD;
  }

  // spectral deconvolution main function
  void FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_()
  {
    deconvolved_spectrum_.clear();
    int current_charge_range = current_max_charge_;
    int tmp_peak_cntr = current_charge_range - min_support_peak_count_;

    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;
    double mass_bin_max_value = std::min(log_mz_peaks_[log_mz_peaks_.size() - 1].logMz - filter_[tmp_peak_cntr], log(current_max_mass_ + (double)avg_.getRightCountFromApex(current_max_mass_) + 1.0));

    double bin_width = bin_width_[ms_level_ - 1];

    mass_bin_min_value_ = log(std::max(1.0, current_min_mass_ - avg_.getAverageMassDelta(current_min_mass_)));
    mz_bin_min_value_ = log_mz_peaks_[0].logMz;

    double mz_bin_max_value = log_mz_peaks_[log_mz_peaks_.size() - 1].logMz;
    Size mass_bin_number = getBinNumber_(mass_bin_max_value, mass_bin_min_value_, bin_width) + 1;
    bin_offsets_.clear();
    harmonic_bin_offset_matrix_.clear();

    for (int i = 0; i < current_charge_range; i++)
    {
      bin_offsets_.push_back((int)round((mz_bin_min_value_ - filter_[i] - mass_bin_min_value_) * bin_width));
    }

    harmonic_bin_offset_matrix_.resize(harmonic_charges_.size(), current_charge_range);
    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      for (int i = 0; i < current_charge_range; i++)
      {
        harmonic_bin_offset_matrix_.setValue(k, i, (int)round((mz_bin_min_value_ - harmonic_filter_matrix_.getValue(k, i) - mass_bin_min_value_) * bin_width));
      }
    }

    Size mz_bin_number = getBinNumber_(mz_bin_max_value, mz_bin_min_value_, bin_width) + 1;
    auto mz_bin_intensities = std::vector<float>(mz_bin_number, .0f);
    updateMzBins_(mz_bin_number, mz_bin_intensities);
    mass_bins_ = boost::dynamic_bitset<>(mass_bin_number);

    if (decoy_flag_ == PeakGroup::charge_decoy && previously_deconved_mono_masses_for_decoy.size() > 0)
    {
      std::sort(previously_deconved_mono_masses_for_decoy.begin(), previously_deconved_mono_masses_for_decoy.end());
      previously_deconved_mass_bins_for_decoy = boost::dynamic_bitset<>(mass_bins_.size());

      for (double m : previously_deconved_mono_masses_for_decoy)
      {
        if (m <= 0)
        {
          continue;
        }
        Size j = getBinNumber_(log(m), mass_bin_min_value_, bin_width_[ms_level_ - 1]);

        if (j > 0 && j < previously_deconved_mass_bins_for_decoy.size() - 1)
        {
          previously_deconved_mass_bins_for_decoy[j] = true;
        }
      }
    }

    if (target_mono_masses_.size() > 0)
    {
      target_mass_bins_.clear();
      target_mass_bins_ = boost::dynamic_bitset<>(mass_bins_.size());
      for (Size i = 0; i < target_mono_masses_.size(); i++)
      {
        double& tm = target_mono_masses_[i];
        for (int off = -1; off < 2; off++)
        {
          double m = tm + off * iso_da_distance_;
          double mass_delta = avg_.getMostAbundantMassDelta(m);

          Size j = getBinNumber_(log(m + mass_delta), mass_bin_min_value_, bin_width_[ms_level_ - 1]);
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

    if (decoy_flag_ != PeakGroup::isotope_decoy)
    {
      auto per_mass_abs_charge_ranges = updateMassBins_(mz_bin_intensities);
      getCandidatePeakGroups_(per_mass_abs_charge_ranges);
    }
    else
    {
      deconvolved_spectrum_ = targetFD_->deconvolved_spectrum_;
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
        auto& peak_group = deconvolved_spectrum_[i];
        int offset = 0;
        float cos = getIsotopeCosineAndDetermineIsotopeIndex(peak_group.getMonoMass(), peak_group.getIsotopeIntensities(), offset, avg_, -1, allowed_iso_error_, decoy_flag_);

        peak_group.setIsotopeCosine(cos);

        if (cos < .5)
        {
          continue;
        }

        peak_group.recruitAllPeaksInSpectrum(deconvolved_spectrum_.getOriginalSpectrum(), tol, avg_, peak_group.getMonoMass() + offset * iso_da_distance_, excluded_integer_mzs_);

        if (peak_group.empty() || peak_group.getMonoMass() < current_min_mass_ || peak_group.getMonoMass() > current_max_mass_)
        {
          continue;
        }

        auto cr = peak_group.getAbsChargeRange();

        if (std::get<0>(cr) > low_charge_ && (std::get<1>(cr) - std::get<0>(cr)) < min_support_peak_count_)
        {
          continue;
        }

        if (target_mono_masses_.size() > 0)
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

        // min cosine is checked in here. mono mass is also updated one last time. SNR, per charge SNR, and avg errors are updated here.
        peak_group.updateIsotopeCosineSNRAvgErrorAndQScore(avg_, peak_group.isTargeted() ? .5 : min_isotope_cosine_[ms_level_ - 1]);

        cr = peak_group.getAbsChargeRange();

        if (std::get<0>(cr) > low_charge_ && (std::get<1>(cr) - std::get<0>(cr)) < min_support_peak_count_)
        {
          continue;
        }

        if (peak_group.getQScore() <= 0)
        {
          continue;
        }

        if (decoy_flag_ == PeakGroup::charge_decoy && previously_deconved_mono_masses_for_decoy.size() > 0)
        {
          bool exclude = false;
          double delta = peak_group.getMonoMass() * tolerance_[ms_level_ - 1];
          auto upper = std::upper_bound(previously_deconved_mono_masses_for_decoy.begin(), previously_deconved_mono_masses_for_decoy.end(), peak_group.getMonoMass() + delta);

          while (!exclude)
          {
            if (upper != previously_deconved_mono_masses_for_decoy.end())
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
            if (upper == previously_deconved_mono_masses_for_decoy.begin())
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

        if (peak_group.getQScore() <= 0)
        {
          continue;
        }
        peak_group.setDecoyFlag(decoy_flag_);
        filtered_peak_groups_private.push_back(peak_group);
      }

#ifdef _OPENMP
#pragma omp for schedule(static) ordered
      for(int i=0; i<omp_get_num_threads(); i++) {
#pragma omp ordered
        filtered_peak_groups.insert(filtered_peak_groups.end(), filtered_peak_groups_private.begin(), filtered_peak_groups_private.end());
      }
#else
      filtered_peak_groups = filtered_peak_groups_private;
#endif
    }
    deconvolved_spectrum_.setPeakGroups(filtered_peak_groups);
    deconvolved_spectrum_.sort();

    removeOverlappingPeakGroups_(deconvolved_spectrum_);
    removeChargeErrorPeakGroups_(deconvolved_spectrum_); //
  }

  float FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(const double mono_mass, const std::vector<float>& per_isotope_intensities, int& offset, const PrecalculatedAveragine& avg,
                                                                       int window_width, int allowed_iso_error_for_second_best_cos, PeakGroup::DecoyFlag decoyFlag)
  {
    offset = 0;
    if (per_isotope_intensities.size() < min_iso_size_)
    {
      return .0;
    }
    auto iso = avg.get(mono_mass);

    int iso_size = (int)iso.size();
    int right = (int)avg.getApexIndex(mono_mass) / 4 + 2;
    int left = right;

    if (window_width >= 0)
    {
      right = std::min(right, window_width);
      left = std::min(left, window_width);
    }

    float max_cos = -1000;
    float second_max_cos = -1000;
    int second_max_offset = -1000;
    // int isotope_length = 0;
    int max_isotope_index = (int)per_isotope_intensities.size(); // exclusive
    int min_isotope_index = -1;                             // inclusive

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

    if (decoyFlag != PeakGroup::DecoyFlag::target)
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

        if (second_max_cos < tmp_cos)
        {
          second_max_cos = tmp_cos;
          second_max_offset = tmp_offset;
        }
      }
      max_cos = second_max_cos;
      offset = second_max_offset;
    }

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
        //n -= a[j] * b[0].getIntensity();
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

  void FLASHDeconvAlgorithm::removeChargeErrorPeakGroups_(DeconvolvedSpectrum& dspec)
  {
    std::map<double, std::set<int>> peak_to_pgs;
    std::set<int> to_remove_pgs;
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(dspec.size());

    for (Size i = 0; i < dspec.size(); i++)
    {
      auto& pg = dspec[i];
      for (auto& p : pg)
      {
        peak_to_pgs[p.mz].insert((int)i);
      }
    }

    for (auto& e : peak_to_pgs)
    {
      auto& pg_is = e.second;
      auto pmz = e.first;
      if (pg_is.size() == 1)
      {
        continue;
      }

      for (auto i : pg_is)
      {
        double mass1 = dspec[i].getMonoMass();
        // int iso_length = avg_.getRightCountFromApex(mass1) + avg_.getLeftCountFromApex(mass1) + 1;

        int repz1 = (int)round(mass1 / pmz);

        double snr1 = dspec[i].getSNR();
        for (auto j : pg_is)
        {
          if (i == j)
          {
            continue;
          }
          double mass2 = dspec[j].getMonoMass();
          int repz2 = (int)round(mass2 / pmz);
          if (repz2 != dspec[j].getRepAbsCharge())
          {
            continue;
          }
          double snr2 = dspec[j].getSNR();

          if (snr1 <= snr2)
          {
            continue;
          }
          bool charge_error = false; // if snr is already highly different, it is charge error..

          if (snr1 > snr2 * 4 || std::max(abs(mass2 - mass1 / repz1 * repz2), abs(mass1 - mass2 / repz2 * repz1)) <= 2 * iso_da_distance_)
          {
            charge_error = true;
          }

          if (charge_error)
          {
            to_remove_pgs.insert(j);
          }
        }

        // lastly check if it is excluded
        if (excluded_masses_.size() > 0)
        {
          auto& peak_group = dspec[i];
          double delta = peak_group.getMonoMass() * tolerance_[ms_level_ - 1] * 2;
          auto upper = std::upper_bound(excluded_masses_.begin(), excluded_masses_.end(), peak_group.getMonoMass() + delta);
          bool exclude = false;
          while (!exclude)
          {
            if (upper != excluded_masses_.end())
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
            if (upper == excluded_masses_.begin())
            {
              break;
            }
            --upper;
          }
          if (exclude)
          {
            to_remove_pgs.insert(i);
          }
        }
      }
    }

    for (Size i = 0; i < dspec.size(); i++)
    {
      if (!dspec[i].isTargeted() && to_remove_pgs.find((int)i) != to_remove_pgs.end())
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

  void FLASHDeconvAlgorithm::removeOverlappingPeakGroups_(DeconvolvedSpectrum& dspec)
  {
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(dspec.size());

    for (Size i = 0; i < dspec.size(); i++)
    {
      if (i > 0)
      {
        if (!dspec[i].isTargeted() && abs(dspec[i - 1].getMonoMass() - dspec[i].getMonoMass()) < 1e-3 && dspec[i - 1].getSNR() >= dspec[i].getSNR())
        {
          continue;
        }
      }
      filtered_pg_vec.push_back(dspec[i]);
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
      int start = -1;
      int end = 1;
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
              precursor_pg.setQScore(smap[2]);
              precursor_pg.setRepAbsCharge((int)smap[1]);
              precursor_pg.updateMonomassAndIsotopeIntensities();
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
  double FLASHDeconvAlgorithm::getMassFromMassBin_(Size mass_bin, double bin_width) const
  {
    return exp(getBinValue_(mass_bin, mass_bin_min_value_, bin_width));
  }

  double FLASHDeconvAlgorithm::getMzFromMzBin_(Size mass_bin, double bin_width) const
  {
    return exp(getBinValue_(mass_bin, mz_bin_min_value_, bin_width));
  }

}