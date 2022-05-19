// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>

namespace OpenMS
{
  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm() :
      DefaultParamHandler("FLASHDeconvAlgorithm")
  {
    prev_mass_bins_ms1_ = std::vector<std::vector<Size>>();
    prev_rts_ms1_ = std::vector<double>();

    prev_mass_bins_ms2_ = std::vector<std::map<int, std::vector<Size>>>();
    prev_rts_ms2_ = std::vector<double>();

    //prev_minbin_logmass_vector_ = std::vector<double>();
    defaults_.setValue("tol",
                       DoubleList{10.0, 10.0},
                       "ppm tolerance for MS1, 2, ... (e.g., -tol 10.0 5.0 to specify 10.0 and 5.0 ppm for MS1 and MS2, respectively)");

    defaults_.setValue("min_mass", 50.0, "minimum mass (Da)");
    defaults_.setValue("max_mass", 100000.0, "maximum mass (Da)");

    defaults_.setValue("min_charge", 2, "minimum charge state for MS1 spectra (can be negative for negative mode)");
    defaults_.setValue("max_charge", 100,
                       "maximum charge state for MS1 spectra (can be negative for negative mode)");

    defaults_.setValue("min_mz", -1.0, "if set to positive value, minimum m/z to deconvolve.");
    defaults_.setValue("max_mz", -1.0, "if set to positive value, maximum m/z to deconvolve.");
    defaults_.setValue("min_rt", -1.0, "if set to positive value, minimum RT to deconvolve.");
    defaults_.setValue("max_rt", -1.0, "if set to positive value, maximum RT to deconvolve.");

    defaults_.setValue("isolation_window", 5.0, "default isolation window with. If the input mzML file does not contain isolation window width information, this width will be used.");
    defaults_.addTag("isolation_window", "advanced");

    defaults_.setValue("min_isotope_cosine",
                       DoubleList{.85, .85},
                       "cosine threshold between avg. and observed isotope pattern for MS1, 2, ... (e.g., -min_isotope_cosine_ 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)");

    defaults_.setValue("min_qscore",
                       .0,
                       "minimum QScore threshold for precursors for MS2. QScore is the probability that a precursor is identified, learned by a logistic regression.");

    defaults_.setValue("min_peaks",
                       IntList{3, 3},
                       "minimum number of peaks of consecutive charge states per MS level."
                       "(e.g., -min_peaks 4 2 to specify 4 and 2 for MS1 and MS2, respectively). "
                       "This affects only for peaks of highly charged peaks (>8). "
                       "The peaks of low charges are detected based on m/z distance between isotopes.");

    defaults_.setValue("max_mass_count",
                       IntList{-1, -1},
                       "maximum mass count per spec for MS1, 2, ... (e.g., -max_mass_count_ 100 50 to specify 100 and 50 for MS1 and MS2, respectively. -1 specifies unlimited)");

    defaults_.setValue("min_intensity", 10.0, "intensity threshold");
    defaults_.setValue("rt_window", 180.0, "RT window for MS1 deconvolution");
    defaultsToParam_();
  }


  FLASHDeconvAlgorithm& FLASHDeconvAlgorithm::operator=(const FLASHDeconvAlgorithm& fd)
  {
    if (this ==& fd)
    {
      return *this;
    }
    //...
    return *this;
  }

  //Calcualte the nominla mass from double mass. Mutiply 0.999497 reduces the rounding error.
  int FLASHDeconvAlgorithm::getNominalMass(const double mass)
  {
    return (int) (mass * 0.999497 + .5);
  }


  //The main function called from outside. precursor_map_for_FLASHIda is used to read FLASHIda information
  const DeconvolvedSpectrum& FLASHDeconvAlgorithm::getDeconvolvedSpectrum(const MSSpectrum& spec,
                                                                     const std::vector<DeconvolvedSpectrum>& survey_scans,
                                                                      const int scan_number,
                                                                      const std::map<int, std::vector<std::vector<double>>>& precursor_map_for_FLASHIda)
  {
    ms_level_ = spec.getMSLevel();
    deconvolved_spectrum_ = DeconvolvedSpectrum(spec, scan_number);

    //for MSn (n>1) register precursor peak and peak group.
    registerPrecursor(survey_scans, precursor_map_for_FLASHIda);

    //rt range of analysis
    if (min_rt_ > 0&&  spec.getRT() < min_rt_)
    {
      return deconvolved_spectrum_;
    }
    if (max_rt_ > 0&&  spec.getRT() > max_rt_)
    {
      return deconvolved_spectrum_;
    }
    //based on MS level, adjust charge and mass ranges. Precursor charge and mass determine those.
    current_min_charge_ = ms_level_ == 1 ? min_abs_charge_ - 5 : 1;
    current_min_charge_ = current_min_charge_ < 1 ? 1 : current_min_charge_;
    current_max_charge_ = deconvolved_spectrum_.getCurrentMaxAbsCharge(max_abs_charge_);//
    current_max_mass_ = deconvolved_spectrum_.getCurrentMaxMass(max_mass_);
    current_min_mass_ = deconvolved_spectrum_.getCurrentMinMass(min_mass_);
    //set universal pattern filter and harmonic pattern filters
    setFilters_();
    //LogMzPeaks are generated from raw peaks
    updateLogMzPeaks_(&spec);
    if (log_mz_peaks_.empty())
    {
      return deconvolved_spectrum_;
    }
    //This is the main FLASHDeconv function in which deconvolution is performed.
    generatePeakGroupsFromSpectrum_();
    for (auto& pg : deconvolved_spectrum_)
    {
      pg.sort();
      //sort(pg.begin(), pg.end());
      pg.setScanNumber(scan_number);
    }
    return deconvolved_spectrum_;
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
    min_support_peak_count_ = param_.getValue("min_peaks");

    bin_width_.clear();
    tolerance_ = param_.getValue("tol");

    for (double& j : tolerance_)
    {
      j *= 1e-6;
      bin_width_.push_back(.5 / j);
    }

    min_isotope_cosine_ = param_.getValue("min_isotope_cosine");
    max_mass_count_ = param_.getValue("max_mass_count");
    rt_window_ = param_.getValue("rt_window");
  }

  const FLASHDeconvHelperStructs::PrecalculatedAveragine& FLASHDeconvAlgorithm::getAveragine()
  {
    return avg_;
  }

  void FLASHDeconvAlgorithm::calculateAveragine(const bool use_RNA_averagine)
  {
    auto generator = new CoarseIsotopePatternGenerator();

    auto iso = use_RNA_averagine ?
                   generator->estimateFromRNAWeight(max_mass_) :
                   generator->estimateFromPeptideWeight(max_mass_);
    iso.trimRight(0.01 * iso.getMostAbundant().getIntensity());

    generator->setMaxIsotope(iso.size());
    avg_ = FLASHDeconvHelperStructs::PrecalculatedAveragine(50, max_mass_, 25, generator, use_RNA_averagine);
    avg_.setMaxIsotopeIndex(iso.size() - 1);
  }

  // generate filters
  void FLASHDeconvAlgorithm::setFilters_()
  {
    filter_.clear();
    harmonic_filter_matrix_.clear();
    int charge_range = current_max_charge_ - current_min_charge_ + 1;
    for (int i = 0; i < charge_range; i++)
    {
      filter_.push_back(log(1.0 / (i + current_min_charge_)));//+
    }

    harmonic_filter_matrix_.resize(harmonic_charges_.size(), charge_range);

    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      int hc = harmonic_charges_[k];
      int n = hc / 2;

      for (int i = 0; i < charge_range; i++)
      {
        harmonic_filter_matrix_
            .setValue(k, i, log(1.0 / (-1.0 * n / hc + (i + current_min_charge_))));// + current_min_charge_
      }
    }
  }

  //Generate uncharged log mz transformated peaks
  void FLASHDeconvAlgorithm::updateLogMzPeaks_(const MSSpectrum* spec)
  {
    log_mz_peaks_.clear();
    double threshold = intensity_threshold_;

    if (spec->size() > max_peak_count_)
    {
      std::vector<double> intensities;
      intensities.reserve(spec->size());
      for (auto& peak : *spec)
      {
        if (min_mz_ > 0&&  peak.getMZ() < min_mz_)
        {
          continue;
        }
        if (max_mz_ > 0&&  peak.getMZ() > max_mz_)
        {
          break;
        }
        if (peak.getIntensity() <= threshold)//
        {
          continue;
        }
        intensities.push_back(peak.getIntensity());
      }

      if (intensities.size() > max_peak_count_)
      {
        std::sort(intensities.begin(), intensities.end());
        threshold = intensities[intensities.size() - 1 - max_peak_count_];
      }
    }

    log_mz_peaks_.reserve(spec->size());

    //threshold = threshold < min_intensity * 2 ? min_intensity * 2 : threshold;
    for (auto& peak : *spec)
    {
      if (min_mz_ > 0&&  peak.getMZ() < min_mz_)
      {
        continue;
      }
      if (max_mz_ > 0&&  peak.getMZ() > max_mz_)
      {
        break;
      }
      if (peak.getIntensity() <= threshold)//
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
    return min_value + bin / bin_width;
  }

  // from value to bin number
  Size FLASHDeconvAlgorithm::getBinNumber_(const double value, const double min_value, const double bin_width)
  {
    if (value < min_value)
    {
      return 0;
    }
    return (Size) (((value - min_value) * bin_width) + .5);
  }

  // From log mz to mz bins. To reduce edge effect, its adjacent bin is set when a bin is set.
  void FLASHDeconvAlgorithm::updateMzBins_(const Size bin_number,
                                           std::vector<float>& mz_bin_intensities)
  {
    //mz_bins_for_edge_effect_ = boost::dynamic_bitset<>(bin_number);
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
      //mz_bins_for_edge_effect_.set(bi);

      mz_bin_intensities[bi] += p.intensity;
    }
    mz_bins_for_edge_effect_ = mz_bins_;
    /*
    for (auto& p : log_mz_peaks_)
    {
      Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_width);
      double delta = (p.logMz - getBinValue_(bi, mz_bin_min_value_, bin_width));

      if (delta > 0)
      {
        if (bi < bin_number - 1&&  !mz_bins_for_edge_effect_[bi + 1])
        {
          mz_bins_for_edge_effect_.set(bi + 1);
          mz_bin_intensities[bi + 1] += p.intensity;
        }
      }
      else if (delta < 0)
      {
        if (bi > 0&&  !mz_bins_for_edge_effect_[bi - 1])
        {
          mz_bins_for_edge_effect_.set(bi - 1);
          mz_bin_intensities[bi - 1] += p.intensity;
        }
      }
    }*/
  }

  //take the mass bins from previous overlapping spectra and put them in the candidate mass bins. Only for MS1
  void FLASHDeconvAlgorithm::unionPrevMassBins_()
  {
    if (mass_bins_.empty())
    {
      return;
    }
    long shift = (long) (round((mass_bin_min_value_) *bin_width_[ms_level_ - 1]));

    if (ms_level_ == 1)
    {
      for (Size i = 0; i < prev_mass_bins_ms1_.size(); i++)
      {
        auto& pmb = prev_mass_bins_ms1_[i];
        if (pmb.empty())
        {
          continue;
        }

        for (Size& index : pmb)
        {
          long j = (long) index - shift;
          if (j < 1)
          {
            continue;
          }
          if ((Size) j >= mass_bins_.size() - 1)
          {
            break;
          }
          mass_bins_[j - 1] = true;
          mass_bins_[j] = true;
          mass_bins_[j + 1] = true;
        }
      }
    }
    else if (ms_level_ == 2)
    {
      if (deconvolved_spectrum_.getPrecursorPeakGroup().empty())
      {
        return;
      }

      int nominal_precursor_mass = getNominalMass(deconvolved_spectrum_.getPrecursorPeakGroup().getMonoMass());
      for (Size i = 0; i < prev_mass_bins_ms2_.size(); i++)
      {
        auto& pmb = prev_mass_bins_ms2_[i];
        if (pmb.empty() || pmb.find(nominal_precursor_mass) == pmb.end())
        {
          continue;
        }

        for (Size& index : pmb[nominal_precursor_mass])
        {
          long j = (long) index - shift;
          if (j < 1)
          {
            continue;
          }
          if ((Size) j >= mass_bins_.size() - 1)
          {
            break;
          }

          mass_bins_[j - 1] = true;
          mass_bins_[j] = true;
          mass_bins_[j + 1] = true;
        }
      }
    }

    target_mass_bins_.clear();
    if (target_masses_.size() > 0)
    {
      target_mass_bins_ = boost::dynamic_bitset<>(mass_bins_.size());
      for (Size i = 0; i < target_masses_.size(); i++)
      {
        double& tm = target_masses_[i];
        for (int off = -1; off < 2; off++)
        {
          double m = tm + off * Constants::ISOTOPE_MASSDIFF_55K_U;
          double mass_delta = avg_.getMostAbundantMassDelta(m);

          Size j = getBinNumber_(log(m + mass_delta), mass_bin_min_value_, bin_width_[ms_level_ - 1]);
          if (j < 1)
          {
            continue;
          }

          if (j >= target_mass_bins_.size() - 1)
          {
            break;
          }

          target_mass_bins_[j - 1] = true;
          target_mass_bins_[j] = true;
          target_mass_bins_[j + 1] = true;
        }
      }
    }
  }

  // Find candidate mass bins from the current spectrum. The runtime of FLASHDeconv is deteremined by this function..
  void FLASHDeconvAlgorithm::updateCandidateMassBins_(std::vector<float>& mass_intensities,
                                                      const std::vector<float>& mz_intensities)
  {
    int charge_range = current_max_charge_ - current_min_charge_ + 1;
    int h_charge_size = (int) harmonic_charges_.size();
    int min_peak_cntr = min_support_peak_count_[ms_level_ - 1] - 1;
    long bin_end = (long) mass_bins_.size();
    auto support_peak_count = std::vector<int>(mass_bins_.size(), 0);// per mass bin how many peaks are present

    Size mz_bin_index = mz_bins_.find_first();
    auto mz_bin_index_reverse = std::vector<Size>();
    mz_bin_index_reverse.reserve(mz_bins_.count());

    // invert mz bins so charges are counted from small to large given a mass
    while (mz_bin_index != mz_bins_.npos)
    {
      mz_bin_index_reverse.push_back(mz_bin_index);
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }

    // to calculate continuous charges, the previous charge value per mass should be stored
    auto prev_charges = std::vector<int>(mass_bins_.size(), charge_range + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prev_intensities = std::vector<float>(mass_bins_.size(), 1.0f);

    double bin_width = bin_width_[ms_level_ - 1];

    for (int i = mz_bin_index_reverse.size() - 1; i >= 0; i--)
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
        long mass_bin_index = mz_bin_index + bin_offsets_[j];

        if (mass_bin_index < 0)
        {
          continue;
        }
        if (mass_bin_index >= bin_end)
        {
          break;
        }

        auto& spc = support_peak_count[mass_bin_index];
        int abs_charge = (j + current_min_charge_);
        float& prev_intensity = prev_intensities[mass_bin_index];
        int& prev_charge = prev_charges[mass_bin_index];
        bool charge_not_continous = prev_charge - j != -1&&  (prev_charge <= charge_range);
        bool pass_first_check = false;

        // intensity ratio between consecutive charges should not exceed the factor.
        float factor = abs_charge <= low_charge_ ? 10.0f : 5.0f + 5.0f * low_charge_ / abs_charge;
        // intensity ratio between consecutive charges for possible harmonic should be within this factor
        float hfactor = factor / 2.0f;
        // intensity of previous charge
        // intensity ratio between current and previous charges
        float intensity_ratio = intensity / prev_intensity;
        intensity_ratio = intensity_ratio < 1 ? 1.0f / intensity_ratio : intensity_ratio;
        float support_peak_intensity = 0;
        // check if peaks of continuous charges are present

        // if charge not continous or intensity ratio is too high reset continuousChargePeakPairCount
        if (charge_not_continous || intensity_ratio > factor)
        {
          spc = 0;
        }
        else
        {
          pass_first_check = true;
          if (spc == 0)
          {
            support_peak_intensity = prev_intensity;
          }
        }

        // for low charges, check isotope peak presence.
        double max_h_intensity = 0;
        if (!pass_first_check&&  abs_charge <= low_charge_)
        {// for low charges
          for(int d=1;d>=-1;d-=2)
          {
            bool iso_exist = false;
            double diff = d * Constants::C13C12_MASSDIFF_U / abs_charge / mz;
            Size next_iso_bin = getBinNumber_(log_mz + diff, mz_bin_min_value_, bin_width);
            if (next_iso_bin > 0&&  next_iso_bin < mz_bins_for_edge_effect_.size()&&  mz_bins_for_edge_effect_[next_iso_bin])
            {
              iso_exist = true;
              pass_first_check = true;
            }

            //harmonic check
            if (iso_exist)
            {
              for (int hc : harmonic_charges_)
              {
                double hdiff = diff / hc * (hc / 2);
                {
                  Size next_harmonic_iso_bin = getBinNumber_(log_mz + hdiff, mz_bin_min_value_, bin_width);
                  double harmonic_intensity = mz_intensities[next_harmonic_iso_bin];
                  if (
                      next_harmonic_iso_bin < mz_bins_for_edge_effect_.size()&& 
                      mz_bins_for_edge_effect_[next_harmonic_iso_bin]&& 
                      intensity / 2.0 < harmonic_intensity )
                  {
                    max_h_intensity = max_h_intensity < harmonic_intensity ? harmonic_intensity : max_h_intensity;
                  }
                }
              }
              if(max_h_intensity > 0) // harmonic
              {
                pass_first_check = false;
              }
              else
              {
                support_peak_intensity += mz_intensities[next_iso_bin];
              }
              break;
            }
          }
        }

        if (pass_first_check)
        {

          if (prev_charge - j == -1)//prev_charge - j == -1)//check harmonic artifacts for high charge ranges
          {
            float max_intensity = intensity;
            float min_intensity = prev_intensity;
            if (prev_intensity <= 1.0)
            {
              max_intensity = intensity;//
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

            for (int k = 0; k < h_charge_size; k++)//
            {
              int hmz_bin_index = mass_bin_index - harmonic_bin_offset_matrix_.getValue(k, j);
              if (hmz_bin_index > 0&&  hmz_bin_index < mz_bins_for_edge_effect_.size()&& 
                  mz_bins_for_edge_effect_[hmz_bin_index])
              {
                float hintensity = mz_intensities[hmz_bin_index];
                if (hintensity > low_threshold&& 
                    hintensity < high_threshold)
                {
                  is_harmonic = true;
                  break;
                }
              }

              if (is_harmonic)
              {
                break;
              }
            }

            if (!is_harmonic)// if it is not harmonic
            {
              mass_intensities[mass_bin_index] += intensity + support_peak_intensity;
              spc++;
              if (spc >= min_peak_cntr || spc >= abs_charge / 2)//
              {
                //double tm = exp(getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width));
                mass_bins_[mass_bin_index] = true;
              }
            }
            else// if harmonic
            {
              mass_intensities[mass_bin_index] -= intensity;
            }
          }
          else if (abs_charge <= low_charge_)// for low charge, include the mass if isotope is present
          {
            mass_intensities[mass_bin_index] += intensity + support_peak_intensity;
            spc++;
            {
              //double tm = exp(getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width));
              mass_bins_[mass_bin_index] = true;
            }
          }
        }
        else if (abs_charge <= low_charge_)// if, for low charge, no isotope peak exists or is harmonic..
        {
          mass_intensities[mass_bin_index] -= intensity + max_h_intensity;
        }
        prev_intensity = intensity;
        prev_charge = j;
      }
    }
  }

  // Subfunction of updateMassBins_. If a peak corresponds to multiple masses, only one mass is selected for the peak based on intensities..
  // mass level harmonic check is also performed in this function
  // it also outputs the charge range of each mass bin
  Matrix<int> FLASHDeconvAlgorithm::filterMassBins_(const std::vector<float>& mass_intensities)
  {
    int charge_range = current_max_charge_ - current_min_charge_ + 1;
    double bin_width = bin_width_[ms_level_ - 1];
    Matrix<int> abs_charge_ranges(2, mass_bins_.size(), INT_MAX);
    for (int i = 0; i < mass_bins_.size(); i++)
    {
      abs_charge_ranges.setValue(1, i, INT_MIN);
    }
    Size mz_bin_index = mz_bins_.find_first();
    long bin_size = (long) mass_bins_.size();

    auto to_skip = mass_bins_.flip();
    mass_bins_.reset();

    while (mz_bin_index != mz_bins_.npos)
    {
      long max_index = -1;
      float max_intensity = -1e11;
      int max_intensity_abs_charge_range = 0;

      for (int j = 0; j < charge_range; j++)
      {
        long mass_bin_index = mz_bin_index + bin_offsets_[j];

        if (mass_bin_index < 0)
        {
          continue;
        }
        if (mass_bin_index >= bin_size)
        {
          break;
        }

        if (target_masses_.size() > 0&&  target_mass_bins_[mass_bin_index])
        {
          float t = mass_intensities[mass_bin_index];

          if (t <= 0)
          {// no signal
            continue;
          }
          max_intensity = 1e38;
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

/*
          std::set<int> unseen({1709, 3084,
                                3212,
                                3474,
                                16730,
                                16413,
                                15794,
                                15495});
          double tm = exp(getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width));//
          if (unseen.find(FLASHDeconvAlgorithm::getNominalMass(tm)) != unseen.end())
          {
            std::cout << "2:  " << tm << std::endl;
          }

*/

          if (t <= 0)
          {// no signal
            continue;
          }



          int abs_charge = (j + current_min_charge_);

          if (false)// for ms n, include all masses.
          {
            abs_charge_ranges
                .setValue(0, mass_bin_index,
                          std::min(abs_charge_ranges.getValue(0, mass_bin_index), j));
            abs_charge_ranges
                .setValue(1, mass_bin_index,
                          std::max(abs_charge_ranges.getValue(1, mass_bin_index), j));
            mass_bins_[mass_bin_index] = true;
            max_intensity = max_intensity < t ? t : max_intensity;
          }
          else if (max_intensity < t)
          {
            bool artifact = false;
            // mass level harmonic, charge off by n artifact removal
            if (abs_charge > low_charge_)
            {
              double original_log_mass = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width);
              double mass = exp(original_log_mass);
              //double diff = Constants::ISOTOPE_MASSDIFF_55K_U / mass;
              //for (int iso_off = 0; iso_off <= 0&&  !artifact; ++iso_off)
              {
                double log_mass = log(mass);// + iso_off * Constants::ISOTOPE_MASSDIFF_55K_U);//original_log_mass + diff * iso_off;
                if (log_mass < 1)
                {
                  continue;
                }

                for (int h : harmonic_charges_)
                {
                  for (int f = -1; f <= 1&&  !artifact; f += 2)//
                  {
                    double hmass = log_mass - log(h) * f;
                    Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_, bin_width);
                    if (hmass_index > 0&&  hmass_index < mass_bins_.size() - 1)
                    {
                      if (mass_intensities[hmass_index] >= t)
                      {
                        artifact = true;//
                        break;
                      }
                    }
                  }
                  if (artifact)
                  {
                    break;
                  }
                }
                //   max_intensity_abs_charge off by one here
                if (!artifact)//
                {
                  for (int coff = 1; coff <= 3&&  !artifact; coff++)
                  {
                    for (int f = -1; f <= 1&&  !artifact; f += 2)
                    {
                      if (abs_charge + f * coff <= 0)
                      {
                        continue;
                      }
                      double hmass = log_mass - log(abs_charge) + log(abs_charge + f * coff);
                      Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_, bin_width);
                      if (hmass_index > 0&&  hmass_index < mass_bins_.size() - 1)
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

      if (max_index >= 0&&  max_index < bin_size)
      {
        abs_charge_ranges
            .setValue(0, max_index,
                      std::min(abs_charge_ranges.getValue(0, max_index), max_intensity_abs_charge_range));
        abs_charge_ranges
            .setValue(1, max_index,
                      std::max(abs_charge_ranges.getValue(1, max_index), max_intensity_abs_charge_range));

        mass_bins_[max_index] = true;
      }
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }

    return abs_charge_ranges;
  }

  //update mass bins which will be used to select peaks in the input spectrum...
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins_(const std::vector<float>& mz_intensities)
  {
    auto mass_intensities = std::vector<float>(mass_bins_.size(), 0);
    updateCandidateMassBins_(mass_intensities, mz_intensities);
    long shift = (long) (round((mass_bin_min_value_) *bin_width_[ms_level_ - 1]));
    auto per_mass_abs_charge_ranges = filterMassBins_(mass_intensities);

    return per_mass_abs_charge_ranges;
  }

  //With mass_bins_ from updateMassBins_ function, select peaks from the same mass in the original input spectrum
  void FLASHDeconvAlgorithm::getCandidatePeakGroups_(const Matrix<int>& per_mass_abs_charge_ranges)
  {
    double bin_width = bin_width_[ms_level_ - 1];
    double tol = tolerance_[ms_level_ - 1];
    int charge_range = current_max_charge_ - current_min_charge_ + 1;
    Size mass_bin_size = mass_bins_.size();
    int log_mz_peak_size = (int) log_mz_peaks_.size();
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
    // main interation. per_mass_abs_charge_ranges gives the range of charges for each mass bin
    while (mass_bin_index != mass_bins_.npos)
    {
      double log_m = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width);
      double mass = exp(log_m);
      PeakGroup pg(current_min_charge_,
                   per_mass_abs_charge_ranges.getValue(1, mass_bin_index) + current_min_charge_,//
                   is_positive_);                                                               // make a empty peakGroup (mass)

      pg.reserve(charge_range * 128);
      // the range of isotope span. For a given peak the peaks within the span are searched.
      Size right_index = avg_.getRightCountFromApex(mass);
      Size left_index = avg_.getLeftCountFromApex(mass);

      double total_signal_pwr = 0;
      std::vector<double> total_harmonic_pwr(harmonic_charges_.size(), .0);

      // scan through charge - from mass to m/z
      for (int j = per_mass_abs_charge_ranges.getValue(0, mass_bin_index);
           j <= per_mass_abs_charge_ranges.getValue(1, mass_bin_index);
           j++)
      {
        int max_peak_index = -1;
        int abs_charge = j + current_min_charge_;//
        int& bin_offset = bin_offsets_[j];

        {
          int b_index = mass_bin_index - bin_offset;// m/z bin
          int& cpi = current_peak_index[j];         // in this charge which peak is to be considered?
          double max_intensity = -1;

          while (cpi < log_mz_peak_size - 1)// scan through peaks from cpi
          {
            if (peak_bin_numbers[cpi] == b_index)// if the peak of consideration matches to this mass with charge abs_charge
            {
              double intensity = log_mz_peaks_[cpi].intensity;
              if (intensity >
                  max_intensity)// compare with other matching peaks and select the most intense peak (in max_peak_index)
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
          if (max_peak_index > 0 && max_peak_index <= log_mz_peak_size &&
              peak_bin_numbers[max_peak_index - 1] == b_index - 1 &&
              log_mz_peaks_[max_peak_index - 1].intensity > max_intensity)
          {
            continue;
          }
          if (max_peak_index >= 0 && max_peak_index < log_mz_peak_size - 1 &&
              peak_bin_numbers[max_peak_index + 1] == b_index + 1 &&
              log_mz_peaks_[max_peak_index + 1].intensity > max_intensity)
          {
            continue;
          }
        }

        // now we have a mathcing peak for this mass of charge  abs_charge. From here, istope peaks are collected
        const double mz = log_mz_peaks_[max_peak_index].mz;
        const double iso_delta = Constants::ISOTOPE_MASSDIFF_55K_U / (abs_charge);
        double mz_delta = tol * mz;//

        double peak_pwr = .0;
        double max_noise_peak_intensity = .0;
        double max_mz = mz;
        double max_peak_intensity = log_mz_peaks_[max_peak_index].intensity;

        for (int peak_index = max_peak_index; peak_index < log_mz_peak_size; peak_index++)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          const double intensity = log_mz_peaks_[peak_index].intensity;
          double mz_diff = observed_mz - mz;
          int tmp_i = (int) round(mz_diff / iso_delta);

          if (observed_mz - max_mz > right_index * iso_delta + mz_delta)
          {
            break;
          }

          if (abs(mz_diff - tmp_i * iso_delta) < mz_delta)// if peak is signal
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size)
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = abs_charge;
              p.isotopeIndex = tmp_i;
              pg.push_back(p);
              peak_pwr += max_noise_peak_intensity * max_noise_peak_intensity;
              peak_pwr += intensity * intensity;
              total_signal_pwr += intensity * intensity;
              max_noise_peak_intensity = .0;
            }
            if (max_peak_intensity < intensity)
            {
              max_peak_intensity = intensity;
              max_mz = observed_mz;
            }
          }
          else// if(abs_charge <= low_charge_)
          {
            for (int l = 0; l < harmonic_charges_.size(); l++)
            {
              int hc = harmonic_charges_[l];
              double hiso_delta = iso_delta / (hc / (hc / 2));
              int tmp_hi = (int) round(mz_diff / hiso_delta);
              double err = abs(mz_diff - tmp_hi * hiso_delta);
              if (err < mz_delta)
              {
                total_harmonic_pwr[l] += intensity * intensity;
              }
            }
            max_noise_peak_intensity = max_noise_peak_intensity < intensity ? intensity : max_noise_peak_intensity;
          }
        }
        max_noise_peak_intensity = .0;

        for (int peak_index = max_peak_index - 1; peak_index >= 0; peak_index--)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          const double intensity = log_mz_peaks_[peak_index].intensity;
          //observedMz = mz + isof * i * d - d * mzDelta;
          double mz_diff = observed_mz - mz;
          int tmp_i = (int) round(mz_diff / iso_delta);

          if (max_mz - observed_mz > left_index * iso_delta + mz_delta)
          {
            break;
          }
          if (abs(mz_diff - tmp_i * iso_delta) < mz_delta)
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size)
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.abs_charge = abs_charge;
              p.isotopeIndex = tmp_i;
              pg.push_back(p);
              peak_pwr += max_noise_peak_intensity * max_noise_peak_intensity;
              peak_pwr += intensity * intensity;
              total_signal_pwr += intensity * intensity;
              max_noise_peak_intensity = .0;
            }
            if (max_peak_intensity < intensity)
            {
              max_peak_intensity = intensity;
              max_mz = observed_mz;
            }
          }
          else// if(abs_charge <= low_charge_)
          {
            for (int l = 0; l < harmonic_charges_.size(); l++)
            {
              int hc = harmonic_charges_[l];
              double hiso_delta = iso_delta / (hc / (hc / 2));
              int tmp_hi = (int) round(mz_diff / hiso_delta);
              double err = abs(mz_diff - tmp_hi * hiso_delta);

              if (err < mz_delta)
              {
                total_harmonic_pwr[l] += intensity * intensity;
              }
            }

            max_noise_peak_intensity = max_noise_peak_intensity < intensity ? intensity : max_noise_peak_intensity;
          }
        }
        pg.setChargePower(abs_charge, peak_pwr);
      }

      if (total_signal_pwr / 2.0 > *std::max_element(total_harmonic_pwr.begin(), total_harmonic_pwr.end())&& 
          !pg.empty())
      {
        double max_intensity = -1.0;
        //double sum_intensity = .0;
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
        double iso_tolerance = tol * t_mass;//
        int min_off = 10000;
        int max_off = -1;
        int max_charge = -1;
        std::vector<double> signal_power(current_max_charge_ + 1, .0);

        for (auto& p : pg)
        {
          p.isotopeIndex = round((p.getUnchargedMass() - t_mass) / Constants::ISOTOPE_MASSDIFF_55K_U);
          if (abs(t_mass - p.getUnchargedMass() + Constants::ISOTOPE_MASSDIFF_55K_U * p.isotopeIndex) >
              iso_tolerance)
          {
            continue;
          }
          new_peaks.push_back(p);
          min_off = min_off > p.isotopeIndex ? p.isotopeIndex : min_off;
          max_off = max_off < p.isotopeIndex ? p.isotopeIndex : max_off;
          max_charge = max_charge < p.abs_charge ? p.abs_charge : max_charge;
          signal_power[p.abs_charge] += p.intensity * p.intensity;
        }
        if (min_off == max_off)
        {
        }
        else
        {
          for (int abs_charge = 0; abs_charge < signal_power.size(); abs_charge++)
          {
            double sp = signal_power[abs_charge];
            if (sp <= 0)
            {
              continue;
            }
            pg.setChargeSignalPower(abs_charge, sp);
          }

          pg.swap(new_peaks);

          for (auto& p : pg)
          {
            p.isotopeIndex -= min_off;
          }
          pg.updateMassesAndIntensity();

          if (pg.getMonoMass() < current_min_mass_ || pg.getMonoMass() > current_max_mass_)
          {
            continue;
          }
          deconvolved_spectrum_.push_back(pg);//
        }
      }
      mass_bin_index = mass_bins_.find_next(mass_bin_index);
    }
  }

  //spectral deconvolution main function
  void FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_()
  {
    std::vector<PeakGroup> empty;
    deconvolved_spectrum_.swap(empty);
    int min_peak_cntr = min_support_peak_count_[ms_level_ - 1] - 1;
    int current_charge_range = current_max_charge_ - current_min_charge_ + 1;
    int tmp_peak_cntr = current_charge_range - min_peak_cntr;

    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;
    double mass_bin_max_value = std::min(
        log_mz_peaks_[log_mz_peaks_.size() - 1].logMz -
            filter_[tmp_peak_cntr],
        log(current_max_mass_ + avg_.getRightCountFromApex(current_max_mass_) + 1));

    double bin_width = bin_width_[ms_level_ - 1];
    tmp_peak_cntr = min_peak_cntr - 1;
    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;

    mass_bin_min_value_ = log(std::max(1.0, current_min_mass_ - avg_.getAverageMassDelta(current_min_mass_)));
    mz_bin_min_value_ = log_mz_peaks_[0].logMz;

    double mz_bin_max_value = log_mz_peaks_[log_mz_peaks_.size() - 1].logMz;
    Size mass_bin_number = getBinNumber_(mass_bin_max_value, mass_bin_min_value_, bin_width) + 1;
    bin_offsets_.clear();
    harmonic_bin_offset_matrix_.clear();
    for (int i = 0; i < current_charge_range; i++)
    {
      bin_offsets_.push_back((int) round((mz_bin_min_value_ - filter_[i] - mass_bin_min_value_) * bin_width));
    }

    harmonic_bin_offset_matrix_.resize(harmonic_charges_.size(), current_charge_range);
    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      for (int i = 0; i < current_charge_range; i++)
      {
        harmonic_bin_offset_matrix_
            .setValue(k,
                      i,
                      (int) round((mz_bin_min_value_ - harmonic_filter_matrix_.getValue(k, i) -
                                   mass_bin_min_value_) *
                                  bin_width));
      }
    }

    Size mz_bin_number = getBinNumber_(mz_bin_max_value, mz_bin_min_value_, bin_width) + 1;
    auto mz_bin_intensities = std::vector<float>(mz_bin_number, .0f);

    updateMzBins_(mz_bin_number, mz_bin_intensities);

    mass_bins_ = boost::dynamic_bitset<>(mass_bin_number);
    unionPrevMassBins_();

    auto per_mass_abs_charge_ranges = updateMassBins_(mz_bin_intensities);

    getCandidatePeakGroups_(per_mass_abs_charge_ranges);

    scoreAndFilterPeakGroups_();

    removeOverlappingPeakGroups_(tolerance_[ms_level_ - 1], 1);
    removeHarmonicsPeakGroups_();

    if (ms_level_ == 1)
    {
      while (!prev_rts_ms1_.empty()&& 
             deconvolved_spectrum_.getOriginalSpectrum().getRT() - prev_rts_ms1_[0] > rt_window_)//
      {
        prev_rts_ms1_.erase(prev_rts_ms1_.begin());
        prev_mass_bins_ms1_.erase(prev_mass_bins_ms1_.begin());
      }
      std::vector<Size> curr_mass_bin;
      curr_mass_bin.reserve(deconvolved_spectrum_.size());
      for (auto& pg : deconvolved_spectrum_)//filteredPeakGroups
      {
        pg.shrink_to_fit();
        LogMzPeak mzPeak;
        double max_int = 0;
        for (auto& p : pg)
        {
          if (max_int > p.intensity)
          {
            continue;
          }
          max_int = p.intensity;
          mzPeak = p;
        }
        if (max_int <= 0)
        {
          continue;
        }
        //double mass_delta = avg_.getAverageMassDelta(pg.getMonoMass());

        Size pg_bin = getBinNumber_(log(mzPeak.getUnchargedMass()), 0, bin_width);
        curr_mass_bin.push_back(pg_bin);
      }

      prev_rts_ms1_.push_back(deconvolved_spectrum_.getOriginalSpectrum().getRT());
      prev_mass_bins_ms1_.push_back(curr_mass_bin);//
      //prev_minbin_logmass_vector_.push_back(mass_bin_min_value_);
      prev_rts_ms1_.shrink_to_fit();
      prev_mass_bins_ms1_.shrink_to_fit();
    }
    else if (ms_level_ == 2&&  !deconvolved_spectrum_.getPrecursorPeakGroup().empty())
    {
      int nominal_precursor_mass = getNominalMass(deconvolved_spectrum_.getPrecursorPeakGroup().getMonoMass());

      while (!prev_rts_ms2_.empty()&& 
             deconvolved_spectrum_.getOriginalSpectrum().getRT() - prev_rts_ms2_[0] > rt_window_)//
      {
        prev_rts_ms2_.erase(prev_rts_ms2_.begin());
        prev_mass_bins_ms2_.erase(prev_mass_bins_ms2_.begin());
      }
      std::map<int, std::vector<Size>> curr_mass_bin;
      curr_mass_bin[nominal_precursor_mass].reserve(deconvolved_spectrum_.size());
      for (auto& pg : deconvolved_spectrum_)//filteredPeakGroups
      {
        pg.shrink_to_fit();

        LogMzPeak mzPeak;
        double max_int = 0;
        for (auto& p : pg)
        {
          if (max_int > p.intensity)
          {
            continue;
          }
          max_int = p.intensity;
          mzPeak = p;
        }
        if (max_int <= 0)
        {
          continue;
        }
        //double mass_delta = avg_.getAverageMassDelta(pg.getMonoMass());

        Size pg_bin = getBinNumber_(log(mzPeak.getUnchargedMass()), 0, bin_width);

        //double mass_delta = avg_.getAverageMassDelta(pg.getMonoMass());
        //Size pg_bin = getBinNumber_(log(pg.getMonoMass() + mass_delta), 0, bin_width);
        curr_mass_bin[nominal_precursor_mass].push_back(pg_bin);
      }

      prev_rts_ms2_.push_back(deconvolved_spectrum_.getOriginalSpectrum().getRT());
      prev_mass_bins_ms2_.push_back(curr_mass_bin);//
      //prev_minbin_logmass_vector_.push_back(mass_bin_min_value_);
      prev_rts_ms2_.shrink_to_fit();
      prev_mass_bins_ms2_.shrink_to_fit();
    }
  }


  double FLASHDeconvAlgorithm::getShapeDiff_(const std::vector<double>& a,
                                             const int& a_start,
                                             const int& a_end,
                                             const IsotopeDistribution& b,
                                             const int& b_size,
                                             const int max_b_index,
                                             const int offset)
  {
    double a_norm = .0, b_norm = .0;
    double diff = 0;

    for (int j = a_start; j <= a_end; j++)
    {
      int i = j - offset;

      if (i < 0 || i >= b_size || b[i].getIntensity() <= 0)
      {
        continue;
      }

      if (i < max_b_index - 2 || i > max_b_index + 2)
      {
        continue;
      }

      a_norm += a[j] * a[j];
      b_norm += b[i].getIntensity() * b[i].getIntensity();
    }

    if (a_norm <= 0 || b_norm <= 0)
    {
      return -1;
    }
    a_norm = sqrt(a_norm);
    b_norm = sqrt(b_norm);

    for (int j = a_start; j <= a_end; j++)
    {
      int i = j - offset;

      if (i < 0)
      {
        double t_diff = (a[j] / a_norm - b[0].getIntensity() / b_norm);
        diff -= (t_diff * t_diff);
      }
      else if (i >= b_size || b[i].getIntensity() <= 0)
      {
        continue;
      }

      double t_diff = (a[j] / a_norm - b[i].getIntensity() / b_norm);
      diff += (t_diff * t_diff);
    }

    return diff;
  }


  double FLASHDeconvAlgorithm::getCosine_(const std::vector<double>& a,
                                          const int& a_start,
                                          const int& a_end,
                                          const IsotopeDistribution& b,
                                          const int& b_size,
                                          const int offset)
  {
    double n = .0, a_norm = .0;
    //int c = 0;
    for (int j = a_start; j <= a_end; j++)
    {
      int i = j - offset;
      a_norm += a[j] * a[j];

      if (i < 0)
      {
        n -= a[j] * b[0].getIntensity();
      }
      else if (i >= b_size || b[i].getIntensity() <= 0)
      {
        continue;
      }
      else
      {
        n += a[j] * b[i].getIntensity();//
      }
    }
    if (a_norm <= 0)
    {
      return 0;
    }
    return n / sqrt(a_norm);
  }


  double FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(const double mono_mass,
                                                                        const std::vector<double>& per_isotope_intensities,
                                                                        int& offset,
                                                                        const PrecalculatedAveragine& avg,
                                                                        bool use_shape_diff)
  {
    auto iso = avg.get(mono_mass);

    int iso_size = (int) iso.size();
    int apex_index = avg.getApexIndex(mono_mass);
    int right = avg.getRightCountFromApex(mono_mass);
    int left = avg.getLeftCountFromApex(mono_mass);
    int iso_range = right + left;

    offset = 0;
    double min_diff = -1000;
    //int isotope_length = 0;
    int max_isotope_index = 0, min_isotope_index = -1;

    for (int i = 0; i < avg.getMaxIsotopeIndex(); i++)
    {
      if (per_isotope_intensities[i] <= 0)
      {
        continue;
      }

      //isotope_length++;
      max_isotope_index = i;
      if (min_isotope_index < 0)
      {
        min_isotope_index = i;
      }
    }

    for (int tmp_offset = -left; tmp_offset <= right; tmp_offset++)
    {
      double tmp_cos = getCosine_(per_isotope_intensities,
                                  min_isotope_index,
                                  max_isotope_index,
                                  iso,
                                  iso_size,//apex_index,
                                  tmp_offset);

      if (min_diff < tmp_cos)
      {
        min_diff = tmp_cos;
        offset = tmp_offset;
      }
    }

    min_diff = 1e10;
    int final_offset = offset;
    if (use_shape_diff)
    {
      for (int tmp_offset = offset - std::max(iso_range / 4, 2); tmp_offset <= offset + std::max(iso_range / 4, 2); tmp_offset++)
      {
        double tmp_diff = getShapeDiff_(per_isotope_intensities,
                                        min_isotope_index,
                                        max_isotope_index,
                                        iso,
                                        iso_size, apex_index,
                                        tmp_offset);
        if (tmp_diff < 0)
        {
          continue;
        }

        if (min_diff > tmp_diff)
        {
          min_diff = tmp_diff;
          final_offset = tmp_offset;
        }
      }
    }
    offset = final_offset;
    return getCosine_(per_isotope_intensities,
                      min_isotope_index,
                      max_isotope_index,
                      iso,
                      iso_size,
                      final_offset);
  }


  bool FLASHDeconvAlgorithm::checkChargeDistribution_(const std::vector<double>& per_charge_intensity)
  {
    if (min_support_peak_count_[ms_level_ - 1] <= 3)
    {
      return true;
    }
    double max_per_charge_intensity = .0;
    double sum_charge_intensity = .0;
    int cntr = 0;
    int non_zero_start = -1, non_zero_end = 0;
    int charge_range = current_max_charge_ - current_min_charge_ + 1;
    for (int i = 0; i < charge_range; i++)
    {
      if (per_charge_intensity[i] > 0)
      {
        sum_charge_intensity += per_charge_intensity[i];
        cntr++;
        max_per_charge_intensity = std::max(max_per_charge_intensity, per_charge_intensity[i]);
        if (non_zero_start < 0)
        {
          non_zero_start = i;
        }
        non_zero_end = i;
      }
    }
    if (cntr == 0)
    {
      return false;
    }

    int prev_charge = non_zero_start;
    int n_r = 0;
    float factor = 1;
    double int_threshold = (sum_charge_intensity / cntr) * factor;
    for (int k = prev_charge + 1; k <= non_zero_end; k++)
    {
      if (per_charge_intensity[k] < int_threshold)
      {
        continue;
      }

      if (k - prev_charge == 1)
      {
        n_r++;
      }
      //spc >= a_charge/2
      if (n_r >= min_support_peak_count_[ms_level_ - 1] - 2 || n_r >= cntr / 2)//
      {
        return true;
      }
      prev_charge = k;
    }
    return false;
  }

  void FLASHDeconvAlgorithm::scoreAndFilterPeakGroups_()
  {
    std::vector<PeakGroup> filtered_peak_groups;
    filtered_peak_groups.reserve(deconvolved_spectrum_.size());
    int charge_range = current_max_charge_ - current_min_charge_ + 1;

    Size max_c = max_mass_count_.size() > ms_level_ - 1 ? max_mass_count_[ms_level_ - 1] : -1;
    if (max_c > 0)
    {
      std::vector<double> intensities;
      intensities.reserve(deconvolved_spectrum_.size());

      for (auto& pg : deconvolved_spectrum_)
      {
        if (pg.getMonoMass() < current_min_mass_ || pg.getMonoMass() > current_max_mass_)
        {
          continue;
        }
        intensities.push_back(pg.getIntensity());
      }
      sort(intensities.begin(), intensities.end());
    }
    int avg_max_index = avg_.getMaxIsotopeIndex();
    for (auto& peak_group : deconvolved_spectrum_)
    {
      auto per_isotope_intensities = std::vector<double>(avg_max_index, 0);
      auto per_abs_charge_intensities = std::vector<double>(charge_range, 0);

      auto indices = calculatePerChargeIsotopeIntensity_(
          per_isotope_intensities, per_abs_charge_intensities,
          avg_max_index, peak_group);

      //bool is_charge_well_distributed = true;
      // is_charge_well_distributed = checkChargeDistribution_(per_abs_charge_intensities); //  could be permanently removed?
      //double tmp = getChargeFitScore_(per_abs_charge_intensities, charge_range);

      // if (!is_charge_well_distributed)
      // {
      //   continue;
      //}

      int offset = 0;
      auto tpg = peak_group[0];
      double cos = getIsotopeCosineAndDetermineIsotopeIndex(tpg.getUnchargedMass(),
                                                            per_isotope_intensities,
                                                            offset, avg_, false);
      peak_group.setIsotopeCosine(cos);
      peak_group.updateMassesAndIntensity(offset, avg_max_index);

      if (peak_group.getMonoMass() < current_min_mass_ || peak_group.getMonoMass() > current_max_mass_)
      {
        continue;
      }

      if (target_masses_.size() > 0)
      {
        double delta = peak_group.getMonoMass() * tolerance_[ms_level_ - 1] * 2;
        auto upper = std::upper_bound(target_masses_.begin(), target_masses_.end(), peak_group.getMonoMass() + delta);

        while (upper != target_masses_.begin()&&  !peak_group.isTargeted())
        {
          --upper;
          if (std::abs(*upper - peak_group.getMonoMass()) < delta)
          {
            peak_group.setTargeted();
          }
        }
      }

      if (peak_group.empty())
      {
        continue;
      }

      if (peak_group.isTargeted())
      {
        if (peak_group.getIsotopeCosine() <= 0.5)
        {
          continue;
        }
      }
      else if (peak_group.size() < 2 || peak_group.getIsotopeCosine() <= min_isotope_cosine_[ms_level_ - 1])
      {
        continue;
      }

      double cs = getChargeFitScore_(per_abs_charge_intensities, 1 + std::get<1>(peak_group.getAbsChargeRange()));

      peak_group.setChargeScore(cs);

      auto iso_dist = avg_.get(peak_group.getMonoMass());
      int iso_size = (int) iso_dist.size();

      auto current_charge_range = peak_group.getAbsChargeRange();
      for (int abs_charge = std::get<0>(current_charge_range);
           abs_charge <= std::get<1>(current_charge_range);
           abs_charge++)
      {
        int j = abs_charge - current_min_charge_;//current_min_charge_;
        if (per_abs_charge_intensities[j] <= 0)
        {
          continue;
        }
        auto current_per_isotope_intensities = std::vector<double>(avg_max_index, 0);

        int min_isotope_index = avg_max_index;
        int max_isotope_index = 0;

        double max_intensity = .0;

        for (auto& peak : peak_group)
        {
          if (peak.abs_charge != abs_charge)
          {
            continue;
          }

          if (peak.isotopeIndex > iso_size)
          {
            continue;
          }

          current_per_isotope_intensities[peak.isotopeIndex] += peak.intensity;
          //sumIntensity += p.intensity;
          min_isotope_index = min_isotope_index < peak.isotopeIndex ? min_isotope_index : peak.isotopeIndex;
          max_isotope_index = max_isotope_index < peak.isotopeIndex ? peak.isotopeIndex : max_isotope_index;

          if (max_intensity < peak.intensity)
          {
            max_intensity = peak.intensity;
            // perChargeMaxIntensity[j] = maxIntensity;
          }
          // sp += p.intensity * p.intensity;
        }
        if (max_intensity <= 0)
        {
          continue;
        }

        double cos_score = getCosine_(current_per_isotope_intensities,
                                      min_isotope_index,
                                      max_isotope_index,
                                      iso_dist,
                                      iso_size,
                                      0);

        // double cos_score_squared = cos_score * cos_score;

        peak_group.setChargeIsotopeCosine(abs_charge, cos_score);
        peak_group.setChargeIntensity(abs_charge, per_abs_charge_intensities[j]);
      }


      peak_group.setAvgPPMError(getAvgPPMError_(peak_group));

      peak_group.updateSNR();
      peak_group.setQScore(-10000);

      for (int abs_charge = std::get<0>(current_charge_range);
           abs_charge <= std::get<1>(current_charge_range);
           abs_charge++)
      {
        if (peak_group.getChargeIntensity(abs_charge) <= 0)
        {
          continue;
        }
        //        int j = abs_charge - current_min_charge_;//current_min_charge_;

        double q_score = QScore::getQScore(&peak_group, abs_charge);

        if (q_score <= peak_group.getQScore())
        {
          continue;
        }
        peak_group.setRepAbsCharge(abs_charge);
        peak_group.setQScore(q_score);
      }
      if (ms_level_ == 1&&  peak_group.getRepAbsCharge() < min_abs_charge_)
      {
        continue;
      }
      auto max_q_score_mz_range = peak_group.getMzRange(peak_group.getRepAbsCharge());
      if (std::get<0>(max_q_score_mz_range) > std::get<1>(max_q_score_mz_range))
      {
        continue;
      }
      peak_group.setMaxQScoreMzRange(std::get<0>(max_q_score_mz_range), std::get<1>(max_q_score_mz_range));
      filtered_peak_groups.push_back(peak_group);
    }

    deconvolved_spectrum_.swap(filtered_peak_groups);
    filterPeakGroupsByIsotopeCosine_(max_c);
  }

  void FLASHDeconvAlgorithm::filterPeakGroupsByIsotopeCosine_(const int current_max_mass_count)
  {
    if (current_max_mass_count <= 0 || deconvolved_spectrum_.size() <= (Size) current_max_mass_count)
    {
      return;
    }

    std::vector<double> iso_cos_scores;
    iso_cos_scores.reserve(deconvolved_spectrum_.size());
    for (auto& pg : deconvolved_spectrum_)
    {
      iso_cos_scores.push_back(pg.getIsotopeCosine());
    }

    sort(iso_cos_scores.begin(), iso_cos_scores.end());

    auto new_peak_groups = std::vector<PeakGroup>();
    new_peak_groups.reserve(deconvolved_spectrum_.size());
    double threshold = iso_cos_scores[iso_cos_scores.size() - current_max_mass_count];

    for (auto& pg : deconvolved_spectrum_)
    {
      if (new_peak_groups.size() > current_max_mass_count)
      {
        break;
      }
      if (pg.isTargeted())
      {
        new_peak_groups.push_back(pg);
      }
    }

    for (auto& pg : deconvolved_spectrum_)
    {
      if (new_peak_groups.size() > current_max_mass_count)
      {
        break;
      }
      if (pg.isTargeted())
      {
        continue;
      }
      if (pg.getIsotopeCosine() >= threshold)
      {
        new_peak_groups.push_back(pg);
      }
    }
    deconvolved_spectrum_.swap(new_peak_groups);
  }

  void FLASHDeconvAlgorithm::removeHarmonicsPeakGroups_()
  {
    std::map<double, std::set<int>> peak_to_pgs;
    std::set<int> to_remove_pgs;
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(deconvolved_spectrum_.size());

    for (int i = 0; i < deconvolved_spectrum_.size(); i++)
    {
      PeakGroup pg = deconvolved_spectrum_[i];
      for (auto& p : pg)
      {
        peak_to_pgs[p.mz].insert(i);
      }
    }

    for (auto& e : peak_to_pgs)
    {
      auto pg_is = e.second;
      if (pg_is.size() == 1)
      {
        continue;
      }

      for (auto i : pg_is)
      {
        double mass1 = deconvolved_spectrum_[i].getMonoMass();
        double snr1 = deconvolved_spectrum_[i].getSNR();
        for (auto j : pg_is)
        {
          double mass2 = deconvolved_spectrum_[j].getMonoMass();
          double snr2 = deconvolved_spectrum_[j].getSNR();
          if (mass1 <= mass2)
          {
            continue;
          }

          bool harmonic = false;
          for (int hc = 2; hc <= 10; hc++)
          {
            if (abs(hc - mass1 / mass2) < 0.01)
            {
              harmonic = true;
              break;
            }
          }

          if (harmonic)
          {
            if (snr1 < snr2)
            {
              to_remove_pgs.insert(i);
            }
            else
            {
              to_remove_pgs.insert(j);
            }
          }
        }
      }
    }

    for (int i = 0; i < deconvolved_spectrum_.size(); i++)
    {
      if (!deconvolved_spectrum_[i].isTargeted()&&  to_remove_pgs.find(i) != to_remove_pgs.end())
      {
        continue;
      }
      filtered_pg_vec.push_back(deconvolved_spectrum_[i]);
    }
    deconvolved_spectrum_.swap(filtered_pg_vec);
  }

  void FLASHDeconvAlgorithm::removeOverlappingPeakGroups_(const double tol, const int iso_length)
  {
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(deconvolved_spectrum_.size());
    deconvolved_spectrum_.sort();
    //sort(deconvolved_spectrum_.begin(), deconvolved_spectrum_.end());
    for (Size i = 0; i < deconvolved_spectrum_.size(); i++)
    {
      if (i > 0)
      {
        if (!deconvolved_spectrum_[i].isTargeted()&& 
            abs(deconvolved_spectrum_[i - 1].getMonoMass() - deconvolved_spectrum_[i].getMonoMass()) < 1e-3&& 
            deconvolved_spectrum_[i - 1].getIntensity() >= deconvolved_spectrum_[i].getIntensity())
        {
          continue;
        }
      }

      filtered_pg_vec.push_back(deconvolved_spectrum_[i]);
    }
    deconvolved_spectrum_.swap(filtered_pg_vec);
    std::vector<PeakGroup>().swap(filtered_pg_vec);
    filtered_pg_vec.reserve(deconvolved_spectrum_.size());

    for (Size i = 0; i < deconvolved_spectrum_.size(); i++)
    {
      bool select = true;
      auto& pg = (deconvolved_spectrum_)[i];
      if (pg.isTargeted())
      {
        filtered_pg_vec.push_back(pg);
        continue;
      }

      if (pg.getMonoMass() <= 0)
      {
        continue;
      }
      double mass_tolerance = pg.getMonoMass() * tol;

      int j = i + 1;
      for (int l = 0; l <= iso_length; l++)
      {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j < deconvolved_spectrum_.size(); j++)
        {
          auto& pgo = (deconvolved_spectrum_)[j];

          if (l != 0&&  pgo.getMonoMass() - pg.getMonoMass() < off - mass_tolerance)
          {
            continue;
          }

          if (pgo.getMonoMass() - pg.getMonoMass() > off + mass_tolerance)
          {
            break;
          }
          select &= pg.getSNR() >= pgo.getSNR();
          if (!select)
          {
            break;
          }
        }
        if (!select)
        {
          break;
        }
      }

      if (!select)
      {
        continue;
      }

      j = i - 1;
      for (int l = 0; l <= iso_length; l++)
      {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j >= 0; j--)
        {
          auto& pgo = (deconvolved_spectrum_)[j];

          if (l != 0&&  pg.getMonoMass() - pgo.getMonoMass() < off - mass_tolerance)
          {
            continue;
          }

          if (pg.getMonoMass() - pgo.getMonoMass() > off + mass_tolerance)
          {
            break;
          }
          select &= pg.getSNR() >= pgo.getSNR();
          if (!select)
          {
            break;
          }
        }
      }
      if (!select)
      {
        continue;
      }
      filtered_pg_vec.push_back(pg);
    }
    deconvolved_spectrum_.swap(filtered_pg_vec);
  }


  std::vector<int> FLASHDeconvAlgorithm::calculatePerChargeIsotopeIntensity_(
      std::vector<double>& per_isotope_intensity,
      std::vector<double>& per_charge_intensity,
      const int max_isotope_count,
      PeakGroup& pg)
  {
    int min_pg_charge = INT_MAX;
    int max_pg_charge = INT_MIN;
    //    double maxIntensity = -1;
    int max_intensity_charge_index = -1;
    //    double maxIntensity2 = -1;
    int max_intensity_iso_index = -1;

    for (auto& p : pg)
    {
      if (p.isotopeIndex < 0 || p.isotopeIndex >= max_isotope_count)
      {
        continue;
      }
      min_pg_charge = std::min(min_pg_charge, p.abs_charge);
      max_pg_charge = std::max(max_pg_charge, p.abs_charge);

      int index = p.abs_charge - current_min_charge_;//current_min_charge_;
      per_isotope_intensity[p.isotopeIndex] += p.intensity;
      per_charge_intensity[index] += p.intensity;
    }
    pg.setAbsChargeRange(min_pg_charge, max_pg_charge);
    //std::cout<<" * " << min_pg_charge << " " << max_pg_charge<<std::endl;
    return std::vector<int>{max_intensity_charge_index, max_intensity_iso_index};
  }

  double FLASHDeconvAlgorithm::getCosine_(const std::vector<double>& a, const std::vector<double>& b, const int off)
  {
    double n = .0, d1 = .0, d2 = .0;
    Size size = a.size();
    //int overlapCntr = 0;
    for (Size j = off; j < size - off; j++)
    {
      d1 += a[j] * a[j];
      d2 += b[j] * b[j];
      n += a[j] * b[j];
      //  if(a[j] > 0&&  b[j] > 0) overlapCntr++;
    }

    //if(overlapCntr < 2) return 0; //
    double d = (d1 * d2);
    if (d <= 0 || n <= 0)
    {
      return 0;
    }

    return n / sqrt(d);
  }

  double FLASHDeconvAlgorithm::getChargeFitScore_(const std::vector<double>& per_charge_intensity, Size len)
  {
    double max_per_charge_intensity = .0;
    double summed_intensity = .0;
    int max_index = -1;
    int first_index = -1;
    int last_index = -1;

    for (int i = 0; i < len; i++)
    {
      summed_intensity += per_charge_intensity[i];
      if (per_charge_intensity[i] <= 0)
      {
        if (first_index < 0)
        {
          first_index = i;
        }
        last_index = i;
      }

      if (max_per_charge_intensity > per_charge_intensity[i])
      {
        continue;
      }
      max_per_charge_intensity = per_charge_intensity[i];
      max_index = i;
    }

    if(summed_intensity <= 0)
    {
      return .0;
    }

    first_index = first_index < 0 ? 0 : first_index;

    double p = .0;
    for (int i = max_index; i < last_index - 1; i++)
    {
      double diff = per_charge_intensity[i + 1] - per_charge_intensity[i];
      //double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i + 1]);
      if (diff <= 0)
      {
        continue;
      }
      p += abs(diff);
    }

    for (int i = max_index; i > first_index; i--)
    {
      double diff = per_charge_intensity[i - 1] - per_charge_intensity[i];
      //      double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i - 1]);

      if (diff <= 0)
      {
        continue;
      }
      p += abs(diff);
    }
    return std::max(.0, 1.0 - p / summed_intensity);
  }

  float FLASHDeconvAlgorithm::getAvgPPMError_(PeakGroup pg)
  {
    std::vector<float> diffs;
    std::vector<std::vector<float>> per_isotope_masses;
    int isotope_end_index = 0;

    for (auto& p : pg)
    {
      isotope_end_index = isotope_end_index < p.isotopeIndex ? p.isotopeIndex : isotope_end_index;
    }
    per_isotope_masses = std::vector<std::vector<float>>(isotope_end_index + 1, std::vector<float>());
    for (auto& p : pg)
    {
      per_isotope_masses[p.isotopeIndex].push_back(p.getUnchargedMass());
    }
    diffs.reserve(pg.size());
    for (int i = 0; i < per_isotope_masses.size(); i++)
    {
      auto& v = per_isotope_masses[i];
      Size n = v.size();
      double average = n >= 2 ?
                           accumulate(v.begin(), v.end(), 0.0) / n :
                           pg.getMonoMass() + i * Constants::ISOTOPE_MASSDIFF_55K_U;//
      for (float& t : v)
      {
        diffs.push_back(pow(1e6 * (t - average) / average, 2.0));
      }
    }
    Size n = diffs.size();
    return n == 0 ? .0 : sqrt(accumulate(diffs.begin(), diffs.end(), 0.0) / n);
  }

  void FLASHDeconvAlgorithm::setTargetMasses(const std::vector<double>& masses)
  {
    target_masses_ = masses;
  }



  bool FLASHDeconvAlgorithm::registerPrecursor(const std::vector<DeconvolvedSpectrum>& survey_scans,
                                               const std::map<int, std::vector<std::vector<double>>>& precursor_map_for_real_time_acquisition)
  {

    if (ms_level_ == 1 ||  (survey_scans.empty() && precursor_map_for_real_time_acquisition.empty()))
    {
      return false;
    }
    deconvolved_spectrum_.setPrecursorIntensity(.0);//
    double start_mz = 0;
    double end_mz = 0;
    //int target_precursor_scan = -1;
    for (auto& precursor: deconvolved_spectrum_.getOriginalSpectrum().getPrecursors())
    {
      for (auto& activation_method: precursor.getActivationMethods())
      {
        deconvolved_spectrum_.setActivationMethod(Precursor::NamesOfActivationMethodShort[activation_method]);
        if (deconvolved_spectrum_.getActivationMethod().compare("HCID"))
        {
          deconvolved_spectrum_.setActivationMethod("HCD");
        }
      }
      deconvolved_spectrum_.setPrecursor(precursor);

      double loffset = precursor.getIsolationWindowLowerOffset();
      double uoffest = precursor.getIsolationWindowUpperOffset();
      loffset = loffset <=0? isolation_window_size_/2.0 : loffset;
      uoffest = uoffest <=0? isolation_window_size_/2.0 : uoffest;

      start_mz = loffset > 100.0 ?
                     loffset :
                     -loffset + precursor.getMZ();
      end_mz = uoffest > 100.0 ?
                   uoffest:
                   uoffest + precursor.getMZ();
    }
    if (!precursor_map_for_real_time_acquisition.empty()&&  deconvolved_spectrum_.getPrecursorPeakGroup().empty())
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
          for (auto& smap: map->second)
          {
            if (abs(start_mz - smap[3]) < .001&&  abs(end_mz - smap[4]) < .001)
            {
              LogMzPeak precursor_log_mz_peak(deconvolved_spectrum_.getPrecursor(), is_positive_);
              precursor_log_mz_peak.abs_charge = (int) smap[1];
              precursor_log_mz_peak.isotopeIndex = 0;
              precursor_log_mz_peak.mass = smap[0];
              precursor_log_mz_peak.intensity = smap[6];

              Precursor precursor(deconvolved_spectrum_.getPrecursor());
              precursor.setCharge(precursor_log_mz_peak.abs_charge);
              precursor.setIntensity(smap[5]);
              deconvolved_spectrum_.setPrecursor(precursor);

              PeakGroup precursor_pg(deconvolved_spectrum_.getPrecursorPeakGroup());
              precursor_pg.push_back(precursor_log_mz_peak);
              precursor_pg.setAbsChargeRange(smap[7], smap[8]);
              precursor_pg.setChargeIsotopeCosine(precursor_log_mz_peak.abs_charge, smap[9]);
              precursor_pg.setChargeSNR(precursor_log_mz_peak.abs_charge, smap[10]);//cnsr
              precursor_pg.setIsotopeCosine(smap[11]);
              precursor_pg.setSNR(smap[12]);
              precursor_pg.setChargeScore(smap[13]);
              precursor_pg.setAvgPPMError(smap[14]);
              precursor_pg.setQScore(smap[2]);
              precursor_pg.setRepAbsCharge((int) smap[1]);
              precursor_pg.updateMassesAndIntensity();
              precursor_pg.setScanNumber(map->first);
              deconvolved_spectrum_.setPrecursorPeakGroup(precursor_pg);

              for (int i = survey_scans.size() - 1; i >= 0; i--)
              {
                auto precursor_spectrum = survey_scans[i];
                if (precursor_spectrum.getScanNumber() != deconvolved_spectrum_.getPrecursorScanNumber())
                {
                  continue;
                }

                auto o_spec = precursor_spectrum.getOriginalSpectrum();
                auto spec_iterator = o_spec.PosBegin(start_mz);
                while (spec_iterator != o_spec.end()&&  spec_iterator->getMZ() < end_mz)
                {
                  double intensity = spec_iterator->getIntensity();
                  spec_iterator++;
                }

                double max_score = -1.0;
                for (auto& pg: precursor_spectrum)
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
                  const LogMzPeak *tmp_precursor = nullptr;
                  // double power = .0;
                  for (auto& tmp_peak: pg)
                  {
                    if (tmp_peak.mz < start_mz)
                    {
                      continue;
                    }
                    if (tmp_peak.mz > end_mz)
                    {
                      break;
                    }
                    //power += tmp_peak.intensity * tmp_peak.intensity;

                    if (tmp_peak.intensity < max_intensity)
                    {
                      continue;
                    }
                    max_intensity = tmp_peak.intensity;
                    tmp_precursor =& tmp_peak;
                    //sum_int += tmp_peak.intensity * tmp_peak.intensity;
                  }

                  if (tmp_precursor == nullptr)
                  {
                    continue;
                  }

                  auto score = pg.getChargeSNR(tmp_precursor->abs_charge);
                  // // most intense one should determine the mass
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

    for (int i = survey_scans.size() - 1; i >= 0; i--)
    {
      auto precursor_spectrum = survey_scans[i];

      if (precursor_spectrum.empty())
      {
        continue;
      }
      auto o_spec = precursor_spectrum.getOriginalSpectrum();
      auto spec_iterator = o_spec.PosBegin(start_mz);

      // double total_power = .0;
      while (spec_iterator != o_spec.end()&&  spec_iterator->getMZ() < end_mz)
      {
        //total_power += spec_iterator->getIntensity() * spec_iterator->getIntensity();
        spec_iterator++;
      }
      for (auto& pg: precursor_spectrum)
      {
        //std::sort(pg.begin(),pg.end());
        if (pg[0].mz > end_mz || pg[pg.size() - 1].mz < start_mz)
        {
          continue;
        }

        //double sum_int = 0;
        double max_intensity = .0;
        const LogMzPeak *tmp_precursor = nullptr;

        int c = int(.5 + pg.getMonoMass() / start_mz);
        //bool contained = true;
        //double power = 0;
        //double i_sum = 0;
        for (auto& tmp_peak: pg)
        {
          if (tmp_peak.abs_charge != c)
          {
            continue;
          }
          //sum += tmp_peak.mz * tmp_peak.intensity;

          if (tmp_peak.mz < start_mz || tmp_peak.mz > end_mz)
          {
            continue;
          }
          // power += tmp_peak.intensity * tmp_peak.intensity;

          if (tmp_peak.intensity < max_intensity)
          {
            continue;
          }
          max_intensity = tmp_peak.intensity;
          tmp_precursor =& tmp_peak;
        }

        if (tmp_precursor == nullptr)// || i_sum <= 0 || sum < start_mz || sum > end_mz)
        {
          continue;
        }
        //cos2θ ||V||2/(N1+sin2θ ||V||2).
        //auto t_cos = pg.getIsotopeCosine();
        //auto precursor_snr = t_cos * t_cos * power
        //                     / (total_power - power + (1 - t_cos) * (1 - t_cos) * power);

        auto score = pg.getChargeSNR(tmp_precursor->abs_charge); // most intense one should determine the mass
        if (score < max_score)
        {
          continue;
        }

        Precursor precursor_pg(deconvolved_spectrum_.getPrecursor());
        precursor_pg.setCharge(tmp_precursor->is_positive ? tmp_precursor->abs_charge
                                                    : -tmp_precursor->abs_charge);

        precursor_pg.setIntensity(tmp_precursor->intensity);
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


}// namespace OpenMS