// --------------------------------------------------------------------------
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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>
#include <Eigen/Dense>

namespace OpenMS
{
  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm() :
      DefaultParamHandler("FLASHDeconvAlgorithm")
  {
    prev_mass_bin_vector_ = std::vector<std::vector<Size>>();
    //prev_minbin_logmass_vector_ = std::vector<double>();
    defaults_.setValue("tol",
                       DoubleList{10.0, 5.0},
                       "ppm tolerance for MS1, 2, ... (e.g., -tol 10.0 5.0 to specify 10.0 and 5.0 ppm for MS1 and MS2, respectively)");

    defaults_.setValue("min_mass", 50.0, "minimum mass (Da)");
    defaults_.setValue("max_mass", 100000.0, "maximum mass (Da)");

    defaults_.setValue("min_charge", 1, "minimum charge state (can be negative for negative mode)");
    defaults_.setValue("max_charge", 100, "maximum charge state (can be negative for negative mode)");

    defaults_.setValue("min_mz", -1.0, "if set to positive value, minimum m/z to deconvolute.");
    defaults_.setValue("max_mz", -1.0, "if set to positive value, maximum m/z to deconvolute.");
    defaults_.setValue("min_rt", -1.0, "if set to positive value, minimum RT to deconvolute.");
    defaults_.setValue("max_rt", -1.0, "if set to positive value, maximum RT to deconvolute.");

    defaults_.setValue("min_isotope_cosine",
                       DoubleList{.75, .75},
                       "cosine threshold between avg. and observed isotope pattern for MS1, 2, ... (e.g., -min_isotope_cosine_ 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)");

    defaults_.setValue("min_peaks",
                       IntList{3, 1},
                       "minimum number of supporting peaks for MS1, 2, ...  (e.g., -min_peaks 3 2 to specify 3 and 2 for MS1 and MS2, respectively)");

    //defaults_.setValue("min_charge_score",
    //                   DoubleList{.0, .0},
    //                   "charge score threshold for MS1, 2, ... (e.g., -min_charge_score 0.7 0.3 to specify 0.7 and 0.3 for MS1 and MS2, respectively)");

    //defaults_.setValue("min_charge_cosine",
    //                   .5,
    //                   "cosine threshold between per-charge-intensity and fitted gaussian distribution (applies only to MS1)");
    defaults_.setValue("max_mass_count",
                       IntList{-1, -1},
                       "maximum mass count per spec for MS1, 2, ... (e.g., -max_mass_count_ 100 50 to specify 100 and 50 for MS1 and MS2, respectively. -1 specifies unlimited)");

    defaults_.setValue("min_mass_count",
                       IntList{-1, -1},
                       "minimum mass count per spec for MS1, 2, ... "
                       "this parameter is only for real time acquisition. "
                       "the parameter may not be satisfied in case spectrum quality is too poor. (e.g., -max_mass_count_ -1 2 to specify no min limit and 2 for MS1 and MS2, respectively. -1 specifies unlimited)");
    defaults_.addTag("min_mass_count", "advanced");

    defaults_.setValue("min_intensity", .0, "intensity threshold");
    defaults_.setValue("RT_window", 20.0, "RT window for MS1 deconvolution");
    defaultsToParam_();
  }

  FLASHDeconvAlgorithm::~FLASHDeconvAlgorithm()
  {
    if (log_mz_peaks_.empty())
    {
      return;
    }
    std::vector<LogMzPeak>().swap(log_mz_peaks_);
  }


  FLASHDeconvAlgorithm &FLASHDeconvAlgorithm::operator=(const FLASHDeconvAlgorithm &fd)
  {
    if (this == &fd)
    {
      return *this;
    }
    //...
    return *this;
  }

  ///Calcualte the nominla mass from double mass. Mutiply 0.999497 reduces the rounding error.
  int FLASHDeconvAlgorithm::getNominalMass(const double mass)
  {
    return (int) (mass * 0.999497 + .5);
  }


  //This function is the main function for the deconvolution. Takes empty DeconvolutedSpectrum and fill it up with peakGroups.
  // DeconvolutedSpectrum contains the recursor peak group for MSn.
  //A peakGroup is the collection of peaks from a single mass (monoisotopic mass). Thus it contains peaks from different charges and iostope indices.
  DeconvolutedSpectrum &FLASHDeconvAlgorithm::getDeconvolutedSpectrum(const MSSpectrum &spec,
                                                                      const std::vector<DeconvolutedSpectrum> &survey_scans,
                                                                      const int scan_number)
  {
    deconvoluted_spectrum_ = DeconvolutedSpectrum(spec, scan_number);

    if (!survey_scans.empty())
    {
      deconvoluted_spectrum_.registerPrecursor(survey_scans);
    }
    if (min_rt_ > 0 && spec.getRT() < min_rt_)
    {
      return deconvoluted_spectrum_;
    }
    if (max_rt_ > 0 && spec.getRT() > max_rt_)
    {
      return deconvoluted_spectrum_;
    }

    updateLogMzPeaks_(&spec);
    ms_level_ = spec.getMSLevel();
    //For MS2,3,.. max mass and charge ranges should be determined by precursors
    current_max_charge_ = deconvoluted_spectrum_.getCurrentMaxCharge(max_charge_); // TODO negative mode!!
    current_max_mass_ = deconvoluted_spectrum_.getCurrentMaxMass(max_mass_);
    //Perform deconvolution and fill in deconvoluted_spectrum_
    generatePeakGroupsFromSpectrum_();

    //deconvoluted_spectrum = deconvoluted_spectrum_;
    //Update peakGroup information after deconvolution
    for (auto &pg : deconvoluted_spectrum_)
    {
      sort(pg.begin(), pg.end());
      pg.setScanNumber(scan_number);
    }
    return deconvoluted_spectrum_;
  }

  void FLASHDeconvAlgorithm::updateMembers_()
  {
    min_mz_ = param_.getValue("min_mz");
    max_mz_ = param_.getValue("max_mz");

    min_rt_ = param_.getValue("min_rt");
    max_rt_ = param_.getValue("max_rt");

    min_charge_ = param_.getValue("min_charge");
    max_charge_ = param_.getValue("max_charge");

    if (min_charge_ > max_charge_)
    {
      int tmp = min_charge_;
      min_charge_ = max_charge_;
      max_charge_ = tmp;
    }

    max_mass_ = param_.getValue("max_mass");
    min_mass_ = param_.getValue("min_mass");

    intensity_threshold_ = param_.getValue("min_intensity");
    min_support_peak_count_ = param_.getValue("min_peaks");

    bin_width_.clear();
    tolerance_ = param_.getValue("tol");

    for (int j = 0; j < (int) tolerance_.size(); j++)
    {
      tolerance_[j] *= 1e-6;
      bin_width_.push_back(.5 / tolerance_[j]);
    }

    min_isotope_cosine_ = param_.getValue("min_isotope_cosine");
    //minChargeScore = param_.getValue("min_charge_score");
    //minChargeCosine = param_.getValue("min_charge_cosine");
    max_mass_count_ = param_.getValue("max_mass_count");
    min_mass_count_ = param_.getValue("min_mass_count");
    rt_window_ = param_.getValue("RT_window");
    setFilters_();
  }

  FLASHDeconvHelperStructs::PrecalculatedAveragine FLASHDeconvAlgorithm::getAveragine()
  {
    return avg_;
  }

  void FLASHDeconvAlgorithm::calculateAveragine(const bool use_RNA_averagine)
  {
    avg_ = FLASHDeconvHelperStructs::calculateAveragines(max_mass_, use_RNA_averagine);
  }


  // generate filters
  void FLASHDeconvAlgorithm::setFilters_()
  {
    filter_.clear();
    harmonic_filter_matrix_.clear();
    int charge_range = max_charge_ - min_charge_ + 1;
    for (int i = 0; i < charge_range; i++)
    {
      filter_.push_back(log(1.0 / abs(i + min_charge_)));
    }

    harmonic_filter_matrix_.resize(harmonic_charges_.size(), charge_range);

    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      int hc = harmonic_charges_[k];
      int n = hc / 2;

      for (int i = 0; i < charge_range; i++)
      {
        //auto harmonic_filter_matrix_ = log(1.0 / (i + n / hc + param.minCharge));
        //        harmonic_bin_offset_matrix_[i][k] = (long) round((filter[i] - harmonic_filter_matrix_) * param.bin_width_);

        harmonic_filter_matrix_.setValue(k, i, log(1.0 / (1.0 * n / hc + abs(i + min_charge_))));
      }
    }
  }

  //Generate uncharged log mz transformated peaks
  void FLASHDeconvAlgorithm::updateLogMzPeaks_(const MSSpectrum *spec)
  {
    std::vector<LogMzPeak>().swap(log_mz_peaks_);
    log_mz_peaks_.reserve(spec->size());
    //int index = 0;
    for (auto &peak : *spec)
    {
      if (min_mz_ > 0 && peak.getMZ() < min_mz_)
      {
        continue;
      }
      if (max_mz_ > 0 && peak.getMZ() > max_mz_)
      {
        break;
      }
      if (peak.getIntensity() <= intensity_threshold_)//
      {
        continue;
      }
      LogMzPeak log_mz_peak(peak, min_charge_ > 0);
      log_mz_peaks_.push_back(log_mz_peak);
    }
  }


  double FLASHDeconvAlgorithm::getBinValue_(const Size bin, const double min_value, const double bin_width)
  {
    return min_value + bin / bin_width;
  }

  Size FLASHDeconvAlgorithm::getBinNumber_(const double value, const double min_value, const double bin_width)
  {
    if (value < min_value)
    {
      return 0;
    }
    return (Size) (((value - min_value) * bin_width) + .5);

  }

  // From log mz to mz bins. To reduce edge effect, its adjacent bin is set when a bin is set.
  void FLASHDeconvAlgorithm::updateMzBins_(const Size &bin_number,
                                           std::vector<float> &mz_bin_intensities)
  {
    mz_bins_for_edge_effect_ = boost::dynamic_bitset<>(bin_number);
    mz_bins_ = boost::dynamic_bitset<>(bin_number);
    //std::fill(mz_bin_intensities.begin(), mz_bin_intensities.end(), .0);

    for (auto &p : log_mz_peaks_)
    {
      Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_width_[ms_level_ - 1]);
      if (bi >= bin_number)
      {
        continue;
      }
      mz_bins_.set(bi);
      mz_bins_for_edge_effect_.set(bi);
      mz_bin_intensities[bi] += p.intensity;
    }
    for (auto &p : log_mz_peaks_)
    {
      Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_width_[ms_level_ - 1]);
      double delta = (p.logMz - getBinValue_(bi, mz_bin_min_value_, bin_width_[ms_level_ - 1]));

      if (delta > 0)
      {
        if (bi < bin_number - 1
            && !mz_bins_for_edge_effect_[bi + 1]
            )
        {
          mz_bins_for_edge_effect_.set(bi + 1);
          mz_bin_intensities[bi + 1] += p.intensity;
        }
      }
      else if (delta < 0)
      {
        if (bi > 0
            && !mz_bins_for_edge_effect_[bi - 1]
            )
        {
          mz_bins_for_edge_effect_.set(bi - 1);
          mz_bin_intensities[bi - 1] += p.intensity;
        }
      }
    }
  }

  //take the mass bins from previous overlapping spectra and put them in the candidate mass bins.
  void FLASHDeconvAlgorithm::unionPrevMassBins_()
  {
    if (mass_bins_.empty())
    {
      return;
    }
    for (Size i = 0; i < prev_mass_bin_vector_.size(); i++)
    {
      auto &pmb = prev_mass_bin_vector_[i];
      //getBinNumber_(pg.getMonoMass() + mass_delta, 0, bin_width_[ms_level_ - 1]);
      if (pmb.empty())
      {
        continue;
      }
      long shift = (long) (round((mass_bin_min_value_) * bin_width_[ms_level_ - 1]));

      for (Size &index : pmb)
      {
        long j = (long) index - shift;
        if (j < 0)
        {
          continue;
        }
        if ((Size) j >= mass_bins_.size())
        {
          break;
        }
        mass_bins_[j] = true;
      }
    }
  }

  //Find candidate mass bins from the current spectrum. The runtime of FLASHDeconv is deteremined by this function..
  void FLASHDeconvAlgorithm::updateCandidateMassBins_(std::vector<float> &mass_intensitites,
                                                      const std::vector<float> &mz_intensities)
  {
    int charge_range = current_max_charge_ - min_charge_ + 1;
    int h_charge_size = (int) harmonic_charges_.size();
    int min_peak_cntr = min_support_peak_count_[ms_level_ - 1];
    long bin_end = (long) mass_bins_.size();
    //auto candidateMassBinsForThisSpectrum = boost::dynamic_bitset<>(mass_bins_.size());
    // how many peaks of continuous charges per mass
    auto support_peak_count = std::vector<int>(mass_bins_.size(), 0);

    Size mz_bin_index = mz_bins_.find_first();
    //std::fill(mass_intensitites.begin(), mass_intensitites.end(), 0);

    // to calculate continuous charges, the previous charge value per mass should be stored
    auto prev_charges = std::vector<int>(mass_bins_.size(), charge_range + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prev_intensities = std::vector<float>(mass_bins_.size(), 1.0f);

    double b_width = bin_width_[ms_level_ - 1];

    // intensity change ratio should not exceed the factor.
    float factor = 5.0;
    float hfactor = 1.1;

    while (mz_bin_index != mz_bins_.npos)
    {
      float intensity = mz_intensities[mz_bin_index];
      double mz = -1.0, log_mz = 0;
      if (ms_level_ > 1)
      {
        log_mz = getBinValue_(mz_bin_index, mz_bin_min_value_, b_width);
        mz = exp(log_mz);
      }

      // scan through charges
      for (int j = 0; j < charge_range; j++) // TDDO: charge should be from 0 to 100 or 0 to -100 // abs increasing.
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

        auto &spc = support_peak_count[mass_bin_index];
        int a_charge = abs(j + min_charge_); // abs charge

        if (ms_level_ == 1)
        {
          // intensity of previous charge
          float &prev_intensity = prev_intensities[mass_bin_index];
          // intensity ratio between current and previous charges
          float intensity_ratio = intensity / prev_intensity;
          intensity_ratio = intensity_ratio < 1 ? 1.0f / intensity_ratio : intensity_ratio;

          // check if peaks of continuous charges are present
          int &prev_charge = prev_charges[mass_bin_index];
          bool charge_not_continous = prev_charge - j != 1;
          // if charge not continous or intensity ratio is too high reset continuousChargePeakPairCount
          if (charge_not_continous || intensity_ratio > factor)
          {
            spc = 0;
          }
          else
          { // check harmonic artifacts
            float max_intensity = intensity;
            float min_intensity = prev_intensity;
            if (min_intensity > max_intensity)
            {
              float tmpi = min_intensity;
              min_intensity = max_intensity;
              max_intensity = tmpi;
            }

            float high_threshold = max_intensity * hfactor;
            float low_threshold = min_intensity / hfactor;// / factor;
            bool is_harmonic = false;
            for (int k = 0; k < h_charge_size; k++)//
            {
              for (int off = -1; off <= 1; off++)
              {
                int hmz_bin_index = off + mass_bin_index - harmonic_bin_offset_matrix_.getValue(k, j);
                if (hmz_bin_index > 0 && hmz_bin_index < mz_bins_for_edge_effect_.size() &&
                    mz_bins_for_edge_effect_[hmz_bin_index])
                {
                  float hintensity = mz_intensities[hmz_bin_index];
                  if (hintensity > low_threshold
                      &&
                      hintensity < high_threshold
                      )
                  {
                    is_harmonic = true;
                    break;//
                    //maxHcharge = k;
                  }
                }
              }
              if (is_harmonic)
              {
                break;
              }
            }
            if (!is_harmonic)
            {
              if (spc == 0)
              {
                mass_intensitites[mass_bin_index] += prev_intensity;
              }

              mass_intensitites[mass_bin_index] += intensity;
              if (++spc >= min_peak_cntr)
              {
                mass_bins_[mass_bin_index] = true;
              }
              else if (a_charge >= 3)
              {
                mass_bins_[mass_bin_index] = (spc >= a_charge / 2);
              }
            }
            else
            {
              mass_intensitites[mass_bin_index] -= intensity;
            }
          }
          prev_intensity = intensity;
          prev_charge = j;
        }
        else // for MS2,3,... iostopic peaks or water nh3 loss peaks are considered
        {
          bool support_peak_present = false;
          double iso_intensity = .0;
          double diff = Constants::ISOTOPE_MASSDIFF_55K_U / a_charge / mz;
          Size next_iso_bin = getBinNumber_(log_mz + diff, mz_bin_min_value_, b_width);

          if (next_iso_bin < mz_bins_for_edge_effect_.size() && mz_bins_for_edge_effect_[next_iso_bin])
          {
            iso_intensity = mz_intensities[next_iso_bin];
            support_peak_present = true;
            //spc++;
          }

          if (support_peak_present)
          {
            double water_loss_mz = log(mz - 18.010565 / a_charge); // 17.026549
            Size water_loss_bin = getBinNumber_(water_loss_mz, mz_bin_min_value_, b_width);

            if (water_loss_bin >= 0 && mz_bins_for_edge_effect_[water_loss_bin])
            {
              float water_loss_intensity = mz_intensities[water_loss_bin];
              if (water_loss_intensity < intensity)
              {
                iso_intensity += water_loss_intensity;
              }
            }

            double amonia_loss_mz = log(mz - 17.026549 / a_charge); // 17.026549
            Size amonia_loss_bin = getBinNumber_(amonia_loss_mz, mz_bin_min_value_, b_width);

            if (amonia_loss_bin >= 0 && mz_bins_for_edge_effect_[amonia_loss_bin])
            {
              float amonia_loss_intensity = mz_intensities[amonia_loss_bin];
              if (amonia_loss_intensity < intensity)
              {
                iso_intensity += amonia_loss_intensity;
              }
            }

            mass_intensitites[mass_bin_index] += intensity + iso_intensity;
            mass_bins_[mass_bin_index] = (++spc >= min_peak_cntr);
          }
          else
          {
            mass_intensitites[mass_bin_index] -= intensity;
            //spc = 0;
          }
        }
      }
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }
  }

  // Subfunction of updateMassBins_. If a peak corresponds to multiple masses, only one mass is selected based on intensities..
  Matrix<int> FLASHDeconvAlgorithm::filterMassBins_(const std::vector<float> &mass_intensities)
  {
    //int chargeRange = param.currentChargeRange;
    int charge_range = current_max_charge_ - min_charge_ + 1;
    Matrix<int> charge_ranges(2, mass_bins_.size(), INT_MAX);
    for (int i = 0; i < mass_bins_.size(); i++)
    {
      charge_ranges.setValue(1, i, INT_MIN);
    }
    Size mz_bin_index = mz_bins_.find_first();
    long bin_size = (long) mass_bins_.size();

    //massBinsForThisSpectrum = boost::dynamic_bitset<>(mass_bins_.size());

    auto to_skip = mass_bins_.flip();
    mass_bins_.reset();

    while (mz_bin_index != mz_bins_.npos)
    {
      long max_index = -1;
      float max_count = -1e11;
      int charge = 0;

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

        if (to_skip[mass_bin_index])
        {
          continue;
        }

        float t = mass_intensities[mass_bin_index];

        if (t <= 0)
        { // no signal
          continue;
        }

        if (max_count < t)
        {
          double logMass = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width_[ms_level_]);
          bool harmonic = false;
          for (int h = 2; h <= 6 && !harmonic; h++)
          {
            for (int f = -1; f <= 1 && !harmonic; f += 2) // only minus
            {
              double hmass = logMass - log(h) * f;
              Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_, bin_width_[ms_level_]);
              if (hmass_index > 0 && hmass_index < mass_bins_.size() - 1)
              {
                //for (int off = 0; off <= 0 && !harmonic; off++)
                //{
                if (mass_intensities[hmass_index] >= t)
                {
                  harmonic = true;
                  break;
                }
                //}
              }
            }
          }
          ///   charge off by one here
          if (!harmonic)
          {
            int acharge = abs(j + min_charge_);
            for (int coff = 1; coff <= 3 && !harmonic; coff++)
            {
              for (int f = -1; f <= 1 && !harmonic; f += 2)
              {
                if (acharge + f * coff <= 0)
                {
                  continue;
                }
                double hmass = logMass - log(acharge) + log(acharge + f * coff);
                Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_, bin_width_[ms_level_]);
                if (hmass_index > 0 && hmass_index < mass_bins_.size() - 1)
                {
                  //for (int off = 0; off <= 0 && !harmonic; off++)
                  //{
                  if (mass_intensities[hmass_index] >= t)
                  {
                    //mass_bins_[hmass_index + off] = false;
                    harmonic = true;
                    break;
                  }
                  //}
                }
              }
            }
          }

          if (!harmonic)
          {
            max_count = t;
            max_index = mass_bin_index;
            charge = j;
          }
          else
          {
            to_skip[mass_bin_index] = true;
            mass_bins_[mass_bin_index] = false;
          }
        }
      }

      if (max_index >= 0 && max_index < bin_size)
      {
        charge_ranges.setValue(0, max_index, std::min(charge_ranges.getValue(0, max_index), charge));
        charge_ranges.setValue(1, max_index, std::max(charge_ranges.getValue(1, max_index), charge));
        mass_bins_[max_index] = true;
      }
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }
    return charge_ranges;
  }

  //update mass bins which will be used to select peaks in the input spectrum...
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins_(//float *massIntensities,
      const std::vector<float> &mz_intensities)
  {
    auto mass_intensities = std::vector<float>(mass_bins_.size(), 0);
    updateCandidateMassBins_(mass_intensities, mz_intensities);
    auto per_mass_charge_ranges = filterMassBins_(mass_intensities);

    return per_mass_charge_ranges;
  }

  //With mass_bins_, select peaks from the same mass in the original input spectrum
  void FLASHDeconvAlgorithm::getCandidatePeakGroups_(const Matrix<int> &charge_ranges)
  {
    const int max_missing_isotope = 2;
    double b_width = bin_width_[ms_level_ - 1];
    double tol = tolerance_[ms_level_ - 1];
    int charge_range = current_max_charge_ - min_charge_ + 1;
    Size mass_bin_size = mass_bins_.size();
    int log_mz_peak_size = (int) log_mz_peaks_.size();
    auto current_peak_index = std::vector<int>(charge_range, 0);
    deconvoluted_spectrum_.reserve(mass_bins_.count());
    Size mass_bin_index = mass_bins_.find_first();
    auto peak_bin_numbers = std::vector<Size>(log_mz_peak_size);

    for (int i = 0; i < log_mz_peak_size; i++)
    {
      peak_bin_numbers[i] = getBinNumber_(log_mz_peaks_[i].logMz, mz_bin_min_value_, b_width);
    }

    while (mass_bin_index != mass_bins_.npos)
    {
      double log_m = getBinValue_(mass_bin_index, mass_bin_min_value_, b_width);
      double mass = exp(log_m);
      PeakGroup pg(min_charge_, charge_ranges.getValue(1, mass_bin_index) + min_charge_);

      pg.reserve(charge_range * 30);
      Size right_index = avg_.getIsotopeEndIndex(mass);
      Size left_index = avg_.getIsotopeStartIndex(mass);

      for (int j = charge_ranges.getValue(0, mass_bin_index); j <= charge_ranges.getValue(1, mass_bin_index); j++)
      {
        int &bin_offset = bin_offsets_[j];
        int b_index = mass_bin_index - bin_offset;

        double max_intensity = -1.0;
        int charge = j + min_charge_;
        int &cpi = current_peak_index[j];
        int max_peak_index = -1;

        while (cpi < log_mz_peak_size - 1)
        {
          if (peak_bin_numbers[cpi] == b_index)
          {
            double intensity = log_mz_peaks_[cpi].intensity;
            if (intensity > max_intensity)
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

        const double mz = log_mz_peaks_[max_peak_index].mz;
        const double iso_delta = Constants::ISOTOPE_MASSDIFF_55K_U / abs(charge);
        double mz_delta = tol * mz * 2; //

        double charge_snr = .0;

        int candidate_i = 0;
        // int peakcntr = 0;
        for (int peak_index = max_peak_index; peak_index < log_mz_peak_size; peak_index++)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          const double intensity = log_mz_peaks_[peak_index].intensity;
          //observedMz = mz + isof * i * d - d * mzDelta;
          double mz_diff = observed_mz - mz;

          int tmp_i = (int) (.5 + mz_diff / iso_delta);

          if (tmp_i > (int) right_index)
          {
            break;
          }

          if (tmp_i - candidate_i > max_missing_isotope)
          {
            break;
          }

          if (abs(mz_diff - tmp_i * iso_delta) >= mz_delta) // noise
          {
            charge_snr += intensity * intensity;
            //peakcntr++;
          }
          else
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size)
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.charge = charge;
              pg.push_back(p);
            }
            candidate_i = tmp_i;
          }
        }

        candidate_i = 0;

        for (int peak_index = max_peak_index - 1; peak_index >= 0; peak_index--)
        {
          const double observed_mz = log_mz_peaks_[peak_index].mz;
          const double intensity = log_mz_peaks_[peak_index].intensity;

          //observedMz = mz + isof * i * d - d * mzDelta;
          double mz_diff = mz - observed_mz;
          int tmp_i = (int) (.5 + mz_diff / iso_delta);

          if (tmp_i > (int) left_index)
          {
            break;
          }

          if (tmp_i - candidate_i > max_missing_isotope)
          {
            break;
          }

          if (abs(mz_diff - tmp_i * iso_delta) >= mz_delta)
          {
            charge_snr += intensity * intensity;
            //continue;
          }
          else
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;

            if (bin < mass_bin_size)
            {
              LogMzPeak p(log_mz_peaks_[peak_index]);
              p.charge = charge;
              pg.push_back(p);
            }

            candidate_i = tmp_i;
          }
        }
        if (charge_snr > 0)
        {
          pg.setChargeSNR(charge, charge_snr);
        }
      }

      if (!pg.empty())
      {
        double max_intensity = -1.0;
        double t_max_mass = .0;
        auto new_peaks = std::vector<LogMzPeak>();
        new_peaks.reserve(pg.size());
        for (auto &p : pg)
        {
          if (max_intensity < p.intensity)
          {
            max_intensity = p.intensity;
            t_max_mass = p.getUnchargedMass();
          }
        }
        double iso_tolerance = tol * t_max_mass;
        int min_off = 10000;
        for (auto &p : pg)
        {
          p.isotopeIndex = round((p.getUnchargedMass() - t_max_mass) / Constants::ISOTOPE_MASSDIFF_55K_U);
          if (abs(t_max_mass - p.getUnchargedMass() + Constants::ISOTOPE_MASSDIFF_55K_U * p.isotopeIndex) >
              iso_tolerance)
          {
            continue;
          }
          new_peaks.push_back(p);
          min_off = min_off > p.isotopeIndex ? p.isotopeIndex : min_off;
        }

        pg.swap(new_peaks);
        //std::vector<LogMzPeak>().swap(newPeaks);

        for (auto &p : pg)
        {
          p.isotopeIndex -= min_off;
        }
        pg.updateMassesAndIntensity();
        deconvoluted_spectrum_.push_back(pg); //
      }
      mass_bin_index = mass_bins_.find_next(mass_bin_index);
    }
  }

  //spectral deconvolution main function
  void FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_()
  {
    std::vector<PeakGroup> empty;
    deconvoluted_spectrum_.swap(empty);
    int min_peak_cntr = min_support_peak_count_[ms_level_ - 1];
    int current_charge_range = current_max_charge_ - min_charge_ + 1;
    int tmp_peak_cntr = current_charge_range - min_peak_cntr;

    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;
    double mass_bin_max_value = std::min(
        log_mz_peaks_[log_mz_peaks_.size() - 1].logMz -
        filter_[tmp_peak_cntr],
        log(current_max_mass_ + avg_.getAverageMassDelta(current_max_mass_) + 1));

    double b_width = bin_width_[ms_level_ - 1];
    tmp_peak_cntr = min_peak_cntr - 1;
    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;
    //filter.push_back(log(1.0 / abs(i + minCharge)));

    mass_bin_min_value_ = log(std::max(1.0, min_mass_ - avg_.getAverageMassDelta(min_mass_)));
    mz_bin_min_value_ = log_mz_peaks_[0].logMz;

    double mz_bin_max_value = log_mz_peaks_[log_mz_peaks_.size() - 1].logMz;
    Size mass_bin_number = getBinNumber_(mass_bin_max_value, mass_bin_min_value_, b_width) + 1;

    for (int i = 0; i < current_charge_range; i++)
    {
      bin_offsets_.push_back((int) round((mz_bin_min_value_ - filter_[i] - mass_bin_min_value_) * b_width));
    }

    //long massBinIndex = mzBinIndex + bin_offsets_[j];

    harmonic_bin_offset_matrix_.resize(harmonic_charges_.size(), current_charge_range);
    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      //      std::vector<int> _h_bin_offsets;
      for (int i = 0; i < current_charge_range; i++)
      {
        harmonic_bin_offset_matrix_
            .setValue(k,
                      i,
                      (int) round((mz_bin_min_value_ - harmonic_filter_matrix_.getValue(k, i) - mass_bin_min_value_) *
                                  b_width));
      }
    }

    Size mz_bin_number = getBinNumber_(mz_bin_max_value, mz_bin_min_value_, b_width) + 1;
    auto mz_bin_intensities = std::vector<float>(mz_bin_number, .0f);

    updateMzBins_(mz_bin_number, mz_bin_intensities);

    mass_bins_ = boost::dynamic_bitset<>(mass_bin_number);
    //massBinsForThisSpectrum = boost::dynamic_bitset<>(massBinNumber);

    if (ms_level_ == 1)
    {
      unionPrevMassBins_();
    }
    auto per_mass_charge_ranges = updateMassBins_(//massIntensities,
        mz_bin_intensities);

    getCandidatePeakGroups_(//massIntensities,
        per_mass_charge_ranges);

    scoreAndFilterPeakGroups_();

    removeHarmonicPeakGroups_(tolerance_[ms_level_ - 1]); //

    if (ms_level_ == 1)
    {
      removeOverlappingPeakGroups_(tolerance_[ms_level_ - 1]);
    }
    else
    {
      removeOverlappingPeakGroupsWithNominalMass_();
    }

    if (ms_level_ == 1)
    {
      while (!prev_rt_vector_.empty() &&
             deconvoluted_spectrum_.getOriginalSpectrum().getRT() - prev_rt_vector_[0] > rt_window_)//
      {
        prev_rt_vector_.erase(prev_rt_vector_.begin());
        prev_mass_bin_vector_.erase(prev_mass_bin_vector_.begin());
        //prev_minbin_logmass_vector_.erase(prev_minbin_logmass_vector_.begin());
      }
      std::vector<Size> curr_mass_bin;
      curr_mass_bin.reserve(deconvoluted_spectrum_.size());
      for (auto &pg : deconvoluted_spectrum_)//filteredPeakGroups
      {
        pg.shrink_to_fit();
        //if (massBinsForThisSpectrum[pg.massBinIndex])
        //{
        double mass_delta = avg_.getAverageMassDelta(pg.getMonoMass());
        Size pg_bin = getBinNumber_(pg.getMonoMass() + mass_delta, 0, bin_width_[ms_level_ - 1]);
        curr_mass_bin.push_back(pg_bin);
        //  }
      }

      prev_rt_vector_.push_back(deconvoluted_spectrum_.getOriginalSpectrum().getRT());
      prev_mass_bin_vector_.push_back(curr_mass_bin); //
      //prev_minbin_logmass_vector_.push_back(mass_bin_min_value_);
      prev_rt_vector_.shrink_to_fit();
      prev_mass_bin_vector_.shrink_to_fit();
      // prev_minbin_logmass_vector_.shrink_to_fit();
    }
  }

  double FLASHDeconvAlgorithm::getCosine_(const std::vector<double> &a,
                                          const int &a_start,
                                          const int &a_end,
                                          const IsotopeDistribution &b,
                                          const int &b_size,
                                          const double &b_norm,
                                          const int offset)
  {
    double n = .0, d1 = .0;
    //int c = 0;
    for (int j = a_start; j <= a_end; j++)
    {
      d1 += a[j] * a[j];
      int i = j - offset;
      if (i < 0 || i >= b_size)
      {
        continue;
      }
      n += a[j] * b[i].getIntensity(); //
    }
    double d = (d1 * b_norm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  /*
  double FLASHDeconvAlgorithm::getCosine_(const double *a, double *b, Size size)
  {
    double n = .0, d1 = .0, d2 = .0;
    //int overlapCntr = 0;
    for (Size j = 0; j < size; j++)
    {
      d1 += a[j] * a[j];
      d2 += b[j] * b[j];
      n += a[j] * b[j];
      //  if(a[j] > 0 && b[j] > 0) overlapCntr++;
    }

    //if(overlapCntr < 2) return 0; //
    double d = (d1 * d2);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }*/


  double FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(const double mono_mass,
                                                                        const std::vector<double> &per_isotope_intensities,
                                                                        int &offset,
                                                                        const PrecalculatedAveragine &avg)
  {
    auto iso = avg.get(mono_mass);
    double iso_norm = avg.getNorm(mono_mass);

    int iso_size = (int) iso.size();

    offset = 0;
    double max_cosine = -1;
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

    int max_cntr = 0;
    for (int tmp_offset = -iso_size - min_isotope_index; tmp_offset <= max_isotope_index; tmp_offset++)
    {
      double tmp_cos = getCosine_(per_isotope_intensities,
                                  min_isotope_index,
                                  max_isotope_index,
                                  iso,
                                  iso_size,
                                  iso_norm,
                                  tmp_offset);

      if (max_cosine <= tmp_cos)
      {
        if (max_cosine == tmp_cos)
        {
          max_cntr++;
          offset += tmp_offset;
        }
        else
        {
          max_cosine = tmp_cos;
          max_cntr = 1;
          offset = tmp_offset;
        }
      }
    }
    offset /= max_cntr;
    return max_cosine;
  }


  bool FLASHDeconvAlgorithm::checkChargeDistribution_(const std::vector<double> &per_charge_intensity)
  {
    if (min_support_peak_count_[ms_level_ - 1] <= 1)
    {
      return true;
    }
    double max_per_charge_intensity = .0;
    double sum_charge_intensity = .0;
    int cntr = 0;
    int non_zero_start = -1, non_zero_end = 0;
    int charge_range = current_max_charge_ - min_charge_ + 1;
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
      if (n_r >= min_support_peak_count_[ms_level_ - 1])//
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
    filtered_peak_groups.reserve(deconvoluted_spectrum_.size());
    double min_threshold = std::numeric_limits<double>::max();
    int charge_range = current_max_charge_ - min_charge_ + 1;

    Size max_c = max_mass_count_.size() > ms_level_ - 1 ? max_mass_count_[ms_level_ - 1] : -1;
    Size min_c = min_mass_count_.size() > ms_level_ - 1 ? min_mass_count_[ms_level_ - 1] : -1;

    if (max_c > 0 || min_c > 0)
    {
      std::vector<double> intensities;
      intensities.reserve(deconvoluted_spectrum_.size());

      for (auto &pg : deconvoluted_spectrum_)
      {
        if (pg.getMonoMass() < min_mass_ || pg.getMonoMass() > max_mass_)
        {
          continue;
        }
        intensities.push_back(pg.getIntensity());
      }
      sort(intensities.begin(), intensities.end());

      if (intensities.size() > (Size) min_c)
      {
        min_threshold = intensities[intensities.size() - min_c];
      }
    }

    for (auto &peak_group : deconvoluted_spectrum_)
    {
      bool pass = false;
      if (peak_group.getIntensity() >= min_threshold)
      {
        pass = true; //
      }

      auto per_isotope_intensities = std::vector<double>(avg_.getMaxIsotopeIndex(), 0);
      auto per_charge_intensities = std::vector<double>(charge_range, 0);

      auto indices = calculatePerChargeIsotopeIntensity_(
          per_isotope_intensities, per_charge_intensities,
          avg_.getMaxIsotopeIndex(), peak_group);

      double cs = getChargeFitScore(per_charge_intensities, charge_range);

      peak_group.setChargeScore(cs);

      if (ms_level_ == 1)
      {
        bool is_charge_well_distributed = checkChargeDistribution_(per_charge_intensities);
        //double tmp = getChargeFitScore_(per_charge_intensities, charge_range);

        if (!is_charge_well_distributed)
        {
          if (!pass)
          {
            continue;
          }
        }
      }

      int offset = 0;
      double cos = getIsotopeCosineAndDetermineIsotopeIndex(peak_group[0].getUnchargedMass(),
                                                            per_isotope_intensities,
                                                            offset, avg_);
      peak_group.setIsotopeCosine(cos);

      if (peak_group.empty() ||
          (peak_group.getIsotopeCosine() <=
           min_isotope_cosine_[ms_level_ -
                               1]))// (msLevel <= 1 ? param.minIsotopeCosineSpec : param.minIsotopeCosineSpec2)))
      {
        if (!pass)
        {
          continue;
        }
      }

      peak_group.updateMassesAndIntensity(offset, avg_.getMaxIsotopeIndex());
      if (peak_group.getMonoMass() < min_mass_ || peak_group.getMonoMass() > max_mass_)
      {
        continue;
      }
      auto iso_dist = avg_.get(peak_group.getMonoMass());
      double iso_norm = avg_.getNorm(peak_group.getMonoMass());
      int iso_size = (int) iso_dist.size();
      float total_noise = .0;
      float total_signal = .0;
      //auto perChargeMaxIntensity = std::vector<double>(chargeRange);

      auto current_charge_range = peak_group.getChargeRange();
      for (int charge = std::get<0>(current_charge_range); charge <= std::get<1>(current_charge_range); charge++)
      {
        int j = charge - min_charge_;
        if (per_charge_intensities[j] <= 0)
        {
          continue;
        }
        auto current_per_isotope_intensities = std::vector<double>(avg_.getMaxIsotopeIndex(), 0);

        int min_isotope_index = avg_.getMaxIsotopeIndex();
        int max_isotope_index = 0;

        double max_intensity = .0;
        //double sumIntensity = .0;
        double summed_intensity_squares = .0;

        for (auto &peak: peak_group)
        {
          if (peak.charge != charge)
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

          //min_mz = min_mz < p.mz ? min_mz : p.mz;
          // max_mz = max_mz > p.mz ? max_mz : p.mz;
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

        for (int k = min_isotope_index; k <= max_isotope_index; ++k)
        {
          if (k > iso_size)
          {
            break;
          }
          summed_intensity_squares += current_per_isotope_intensities[k] * current_per_isotope_intensities[k];
        }

        double cos_score = getCosine_(current_per_isotope_intensities,
                                      min_isotope_index,
                                      max_isotope_index,
                                      iso_dist,
                                      iso_size,
                                      iso_norm,
                                      0);

        double cos_score_squared = cos_score * cos_score;

        peak_group.setChargeIsotopeCosine(charge, cos_score);
        peak_group.setChargeIntensity(charge, per_charge_intensities[j]);

        double signal = (1 - cos_score_squared) * summed_intensity_squares + peak_group.getChargeSNR(charge) + 1;
        double noise = cos_score_squared * summed_intensity_squares + 1;

        peak_group.setChargeSNR(charge, noise / signal);

        total_noise += signal;
        total_signal += noise;
      }

      peak_group.setSNR(total_signal / total_noise);
      peak_group.setQScore(-10000);

      for (int charge = std::get<0>(current_charge_range); charge <= std::get<1>(current_charge_range); charge++)
      {
        if (peak_group.getChargeIntensity(charge) <= 0)
        {
          continue;
        }
        int j = charge - min_charge_;

        double q_score = QScore::getQScore(&peak_group, charge);

        if (q_score <= peak_group.getQScore())
        {
          continue;
        }
        peak_group.setRepCharge(charge);
        peak_group.setQScore(q_score);
      }

      double max_q_score_mz_start = peak_group.getMonoMass() * 2;
      double max_q_score_mz_end = .0;
      for (auto &tmp_p:peak_group)
      {
        if (tmp_p.charge != peak_group.getRepCharge())
        {
          continue;
        }
        if (tmp_p.isotopeIndex > iso_size)
        {
          continue;
        }

        max_q_score_mz_start = max_q_score_mz_start < tmp_p.mz ? max_q_score_mz_start : tmp_p.mz;
        max_q_score_mz_end = max_q_score_mz_end > tmp_p.mz ? max_q_score_mz_end : tmp_p.mz;
      }
      if (max_q_score_mz_start > max_q_score_mz_end)
      {
        continue;
      }
      peak_group.setMaxQScoreMzRange(max_q_score_mz_start, max_q_score_mz_end);
      filtered_peak_groups.push_back(peak_group);
    }
    deconvoluted_spectrum_.swap(filtered_peak_groups);

    if (ms_level_ > 1)
    {
      filterPeakGroupsByIsotopeCosine_(max_c);
    }
    else
    {
      filterPeakGroupsByQScore_(max_c);
    }
  }

  void FLASHDeconvAlgorithm::filterPeakGroupsByIsotopeCosine_(const int current_max_mass_count)
  {
    if (current_max_mass_count <= 0 || deconvoluted_spectrum_.size() <= (Size) current_max_mass_count)
    {
      return;
    }

    std::vector<double> iso_cos_scores;
    iso_cos_scores.reserve(deconvoluted_spectrum_.size());
    for (auto &pg : deconvoluted_spectrum_)
    {
      iso_cos_scores.push_back(pg.getIsotopeCosine());
    }

    sort(iso_cos_scores.begin(), iso_cos_scores.end());

    auto new_peak_groups = std::vector<PeakGroup>();
    new_peak_groups.reserve(deconvoluted_spectrum_.size());
    double threshold = iso_cos_scores[iso_cos_scores.size() - current_max_mass_count];
    for (auto &pg : deconvoluted_spectrum_)
    {
      if (new_peak_groups.size() > current_max_mass_count)
      {
        break;
      }

      if (pg.getIsotopeCosine() >= threshold)
      {
        new_peak_groups.push_back(pg);
      }
    }
    deconvoluted_spectrum_.swap(new_peak_groups);
  }

  void FLASHDeconvAlgorithm::filterPeakGroupsByQScore_(const int current_max_mass_count)
  {
    if (current_max_mass_count <= 0 || deconvoluted_spectrum_.size() <= (Size) current_max_mass_count)
    {
      return;
    }

    Size mc = (Size) current_max_mass_count;
    std::vector<double> q_scores;
    q_scores.reserve(deconvoluted_spectrum_.size());
    for (auto &pg : deconvoluted_spectrum_)
    {
      q_scores.push_back(pg.getQScore());
    }

    sort(q_scores.begin(), q_scores.end());

    auto new_peak_groups = std::vector<PeakGroup>();
    new_peak_groups.reserve(deconvoluted_spectrum_.size());
    double threshold = q_scores[q_scores.size() - mc];
    for (auto &pg : deconvoluted_spectrum_)
    {
      if (new_peak_groups.size() > mc)
      {
        break;
      }

      if (pg.getQScore() >= threshold)
      {
        new_peak_groups.push_back(pg);
      }
    }
    deconvoluted_spectrum_.swap(new_peak_groups);
  }


  void FLASHDeconvAlgorithm::removeHarmonicPeakGroups_(const double tol)
  {
    sort(deconvoluted_spectrum_.begin(), deconvoluted_spectrum_.end());
    std::vector<PeakGroup> merged_pg_vec;
    merged_pg_vec.reserve(deconvoluted_spectrum_.size());
    std::vector<double> pg_masses;
    pg_masses.reserve(deconvoluted_spectrum_.size());
    for (auto &peakGroup : deconvoluted_spectrum_)
    {
      pg_masses.push_back(peakGroup.getMonoMass());
    }
    for (auto &pg : deconvoluted_spectrum_)
    {
      bool select = true;
      for (int h = 2; h <= 3; h++)
      {
        for (int k = 0; k < 2; k++)
        {
          for (int i = -2; i <= 2; ++i)
          {
            double omass = pg.getMonoMass() + i * Constants::ISOTOPE_MASSDIFF_55K_U;
            double harmonic_mass = k == 0 ? omass * h : omass / h;
            double mass_tolerance = 2 * std::max(omass, harmonic_mass) * tol;
            auto iter = std::lower_bound(pg_masses.begin(), pg_masses.end(), harmonic_mass - mass_tolerance);
            Size j = iter - pg_masses.begin();

            if (j >= 0 && j < deconvoluted_spectrum_.size())
            {
              for (; j < deconvoluted_spectrum_.size(); j++)
              {
                auto &pgo = (deconvoluted_spectrum_)[j];
                if (harmonic_mass - pgo.getMonoMass() > mass_tolerance)
                {
                  continue;
                }

                if (!select || pgo.getMonoMass() - harmonic_mass > mass_tolerance)
                {
                  break;
                }
                select &= pg.getIntensity() >= pgo.getIntensity();
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
          }
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
      merged_pg_vec.push_back(pg);
    }
    deconvoluted_spectrum_.swap(merged_pg_vec);
  }

  void FLASHDeconvAlgorithm::removeOverlappingPeakGroups_(const double tol)
  {
    int iso_length = 1; // inclusive
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(deconvoluted_spectrum_.size());
    sort(deconvoluted_spectrum_.begin(), deconvoluted_spectrum_.end());

    for (Size i = 0; i < deconvoluted_spectrum_.size(); i++)
    {
      bool select = true;
      auto &pg = (deconvoluted_spectrum_)[i];

      if (pg.getMonoMass() <= 0)
      {
        continue;
      }
      double mass_tolerance = pg.getMonoMass() * tol * 2;

      int j = i + 1;
      for (int l = 0; l <= iso_length; l++)
      {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j < deconvoluted_spectrum_.size(); j++)
        {
          auto &pgo = (deconvoluted_spectrum_)[j];
          if (l != 0 && pgo.getMonoMass() - pg.getMonoMass() < off - mass_tolerance)
          {
            continue;
          }

          if (!select || pgo.getMonoMass() - pg.getMonoMass() > off + mass_tolerance)
          {
            break;
          }
          select &= pg.getIsotopeCosine() > pgo.getIsotopeCosine();
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
          auto &pgo = (deconvoluted_spectrum_)[j];

          if (l != 0 && pg.getMonoMass() - pgo.getMonoMass() < off - mass_tolerance)
          {
            continue;
          }

          if (!select || pg.getMonoMass() - pgo.getMonoMass() > off + mass_tolerance)
          {
            break;
          }
          select &= pg.getIsotopeCosine() > pgo.getIsotopeCosine();
        }
      }
      if (!select)
      {
        continue;
      }
      filtered_pg_vec.push_back(pg);
    }
    deconvoluted_spectrum_.swap(filtered_pg_vec);
  }

  void FLASHDeconvAlgorithm::removeOverlappingPeakGroupsWithNominalMass_()
  {
    int iso_length = 1; // inclusive
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(deconvoluted_spectrum_.size());
    sort(deconvoluted_spectrum_.begin(), deconvoluted_spectrum_.end());

    for (Size i = 0; i < deconvoluted_spectrum_.size(); i++)
    {
      bool select = true;
      auto &pg = (deconvoluted_spectrum_)[i];

      if (pg.getMonoMass() <= 0)
      {
        continue;
      }
      int j = i + 1;
      for (int l = 0; l <= iso_length; l++)
      {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j < deconvoluted_spectrum_.size(); j++)
        {
          auto &pgo = (deconvoluted_spectrum_)[j];
          if (l != 0 && getNominalMass(pgo.getMonoMass()) < getNominalMass(pg.getMonoMass() + off))
          {
            continue;
          }

          if (!select || getNominalMass(pgo.getMonoMass()) > getNominalMass(pg.getMonoMass() + off))
          {
            break;
          }
          select &= pg.getIsotopeCosine() > pgo.getIsotopeCosine();
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
          auto &pgo = (deconvoluted_spectrum_)[j];

          if (l != 0 && getNominalMass(pg.getMonoMass()) < getNominalMass(pgo.getMonoMass() + off))
          {
            continue;
          }

          if (!select || getNominalMass(pg.getMonoMass()) > getNominalMass(pgo.getMonoMass() + off))
          {
            break;
          }
          select &= pg.getIsotopeCosine() > pgo.getIsotopeCosine();
        }
      }
      if (!select)
      {
        continue;
      }
      filtered_pg_vec.push_back(pg);
    }
    deconvoluted_spectrum_.swap(filtered_pg_vec);
  }

  void FLASHDeconvAlgorithm::reassignPeaksinPeakGroups_()
  {
    /*
    //auto maxSNR = std::vector<double>(log_mz_peaks_.size() , 0);
    auto max_intensity = std::vector<double>(log_mz_peaks_.size(), 0);
    //auto max_intensity_charge = std::vector<int>(log_mz_peaks_.size() , 0);

    for (auto& pg : *deconvoluted_spectrum_)
    {

      auto intensity = pg.intensity / (pg.maxCharge - pg.minCharge + 1);

      for (auto& p: pg)
      {
        auto idx = p.index;
        //if (maxSNR[idx] < snr)
        //{
        //  maxSNR[idx] = snr;
        //}
        if (max_intensity[idx] < intensity)
        {
          max_intensity[idx] = intensity;
          //max_intensity_charge[idx] = p.charge;
        }
      }
    }

    for (auto& pg : *deconvoluted_spectrum_)
    {
      //auto snr = pg.total_snr_;
      auto intensity = pg.intensity / (pg.maxCharge - pg.minCharge + 1);

      std::vector<LogMzPeak> tmp;
      tmp.swap(pg);
      pg.reserve(tmp.size());
      for (auto& p: tmp)
      {
        auto idx = p.index;
        if (//maxSNR[idx] * .5 > snr
          //||
            max_intensity[idx] > intensity
          //max_intensity_charge[idx] != p.charge
            )
        {
          continue;
        }
        pg.push_back(p);
      }
    }*/
  }


  std::vector<int> FLASHDeconvAlgorithm::calculatePerChargeIsotopeIntensity_(
      std::vector<double> &per_isotope_intensity,
      std::vector<double> &per_charge_intensity,
      const int max_isotope_count,
      PeakGroup &pg)
  {
    int min_pg_charge = INT_MAX;
    int max_pg_charge = INT_MIN;
    //    double maxIntensity = -1;
    int max_intensity_charge_index = -1;
    //    double maxIntensity2 = -1;
    int max_intensity_iso_index = -1;

    for (auto &p : pg)
    {
      if (p.isotopeIndex < 0 || p.isotopeIndex >= max_isotope_count)
      {
        continue;
      }
      min_pg_charge = std::min(min_pg_charge, p.charge);
      max_pg_charge = std::max(max_pg_charge, p.charge);

      int index = p.charge - min_charge_;
      per_isotope_intensity[p.isotopeIndex] += p.intensity;
      per_charge_intensity[index] += p.intensity;
    }
    pg.setChargeRange(min_pg_charge, max_pg_charge);

    return std::vector<int>{max_intensity_charge_index, max_intensity_iso_index};
  }

  double FLASHDeconvAlgorithm::getCosine_(const std::vector<double> &a, const std::vector<double> &b, const int off)
  {
    double n = .0, d1 = .0, d2 = .0;
    Size size = a.size();
    //int overlapCntr = 0;
    for (Size j = off; j < size - off; j++)
    {
      d1 += a[j] * a[j];
      d2 += b[j] * b[j];
      n += a[j] * b[j];
      //  if(a[j] > 0 && b[j] > 0) overlapCntr++;
    }

    //if(overlapCntr < 2) return 0; //
    double d = (d1 * d2);
    if (d <= 0 || n <= 0)
    {
      return 0;
    }

    return n / sqrt(d);
  }


  double FLASHDeconvAlgorithm::getChargeFitScore_(std::vector<double> &per_charge_intensity, int charge_range)
  {
    double maxPerChargeIntensity = .0;
    int max_index = 0;
    std::vector<double> xs;
    std::vector<double> ys;

    xs.reserve(+2);
    ys.reserve(charge_range + 2);

    for (int i = 0; i < charge_range; i++)
    {
      if (maxPerChargeIntensity > per_charge_intensity[i])
      {
        continue;
      }
      maxPerChargeIntensity = per_charge_intensity[i];
      max_index = i;
    }

    double th = maxPerChargeIntensity * .02;// as recommended in the original paper...
    int first = -1, last = 0;

    for (int i = max_index - 1; i >= 0; i--)
    {
      if (per_charge_intensity[i] <= th)
      {
        break;
      }
      first = i;
    }

    for (int i = max_index + 1; i < charge_range; i++)
    {
      if (per_charge_intensity[i] <= th)
      {
        break;
      }
      last = i;
    }

    for (int i = first; i <= last; i++)
    {
      xs.push_back(i);
      ys.push_back((1 + per_charge_intensity[i]));
    }

    if (xs.size() <= 3)
    {
      return 0.5;
    }

    Eigen::Matrix3d m;
    Eigen::Vector3d v;

    double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    double t0 = 0, t1 = 0, t2 = 0;

    for (Size i = 0; i < xs.size(); i++)
    {
      auto &x = xs[i];
      auto y = log(ys[i]);
      s0++;
      s1 += x;
      s2 += x * x;
      s3 += x * x * x;
      s4 += x * x * x * x;
      t0 += y;
      t1 += y * x;
      t2 += y * x * x;
    }
    m(0, 0) = s0;
    m(1, 0) = m(0, 1) = s1;
    m(2, 0) = m(1, 1) = m(0, 2) = s2;
    m(2, 1) = m(1, 2) = s3;
    m(2, 2) = s4;

    auto im = m.inverse();
    v(0) = t0;
    v(1) = t1;
    v(2) = t2;
    //cout<<v<<endl;
    auto abc = im * v;
    //cout<<abc<<endl;
    double mu = -abc(1) / abc(2) / 2;
    double omega = -1 / abc(2) / 2;

    if (omega <= 0)
    {
      return 0;
    }
    std::vector<double> tys;

    for (Size i = 0; i < ys.size(); i++)
    {
      double ty = exp(-(xs[i] - mu) * (xs[i] - mu) / 2 / omega);
      tys.push_back(ty);
    }
    return getCosine_(ys, tys);
  }

  double FLASHDeconvAlgorithm::getChargeFitScore(const std::vector<double> &per_charge_intensity,
                                                 const int charge_range)
  {
    double max_per_charge_intensity = .0;
    double summed_intensity = .0;
    int max_index = -1;
    int first_index = -1;
    int last_index = charge_range - 1;

    for (int i = 0; i < charge_range; i++)
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
    first_index = first_index < 0 ? 0 : first_index;

    double p = .0;
    for (int i = max_index; i < last_index - 1; i++)
    {
      double diff = per_charge_intensity[i + 1] - per_charge_intensity[i];
      double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i + 1]);
      if (diff <= 0 && ratio < 5.0)
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
}
