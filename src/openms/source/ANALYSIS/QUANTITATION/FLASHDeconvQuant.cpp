// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/FLASHDeconvQuant.h>
#include <queue>

// test purpose
#include <iostream>
#include <fstream>

namespace OpenMS
{
  FLASHDeconvQuant::FLASHDeconvQuant():
      ProgressLogger(),
      DefaultParamHandler("FLASHDeconvQuant")
  {
    defaults_.setValue("local_rt_range", 15.0, "RT range where to look for coeluting mass traces", ListUtils::create<String>("advanced"));
    defaults_.setValue("local_mz_range", 6.5, "MZ range where to look for isotopic mass traces", ListUtils::create<String>("advanced")); // (-> decides size of isotopes =(local_mz_range_ * lowest_charge))
    defaults_.setValue("charge_lower_bound", 5, "Lowest charge state to consider");
    defaults_.setValue("charge_upper_bound", 50, "Highest charge state to consider");
    defaults_.setValue("min_mass", 10000, "minimim mass");
    defaults_.setValue("max_mass", 70000, "maximum mass");
    defaults_.setValue("mz_tol", 20, "ppm tolerance for m/z");

    defaults_.setValue("use_smoothed_intensities", "true", "Use LOWESS intensities instead of raw intensities.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("use_smoothed_intensities", ListUtils::create<String>("false,true"));

    defaultsToParam_();

    this->setLogType(CMD);
  }

  FLASHDeconvQuant::~FLASHDeconvQuant(){}

  void FLASHDeconvQuant::updateMembers_()
  {
    local_rt_range_ = (double)param_.getValue("local_rt_range");
    local_mz_range_ = (double)param_.getValue("local_mz_range");

    charge_lower_bound_ = (Size)param_.getValue("charge_lower_bound");
    charge_upper_bound_ = (Size)param_.getValue("charge_upper_bound");
    charge_range_ = charge_upper_bound_-charge_lower_bound_ + 1;

    min_mass_ = (double)param_.getValue("min_mass");
    max_mass_ = (double)param_.getValue("max_mass");

    mz_tolerance_ = (double)param_.getValue("mz_tol");
    mz_tolerance_ *= 1e-6;
    mz_bin_width_ = .5 / mz_tolerance_;

    use_smoothed_intensities_ = param_.getValue("use_smoothed_intensities").toBool();
  }

  void FLASHDeconvQuant::writeFeatureGroupsInFile(std::vector<FeatureGroup>& feat)
  {
    ofstream out;
    OPENMS_LOG_INFO << "writing output..." << outfile_path << endl;
    out.open(outfile_path, ios::out);

    // header
    out << "mono_mass\tquant_value\tcharge_score\tiso_cosine\tcharges\tfwhm_start\tfwhm_end\trt_start\trt_end\n";

    for (auto& m : feat)
    {
      double q = m.getIntensity();
//      if (q==0.0)
//      {
//        for(auto& idx : m.second.feature_idx)
//        {
//          q += hypotheses[idx].computeQuant();
//        }
//      }

      out << to_string(m.getMonoisotopicMass()) << "\t"
          << to_string(q) << "\t"
          << to_string(m.getChargeScore()) << "\t"
          << to_string(m.getIsotopeCosine()) << "\t";

      std::string charges = to_string(std::get<0>(m.getChargeRange())) + ","
                    + to_string(std::get<1>(m.getChargeRange()));

      // rt range
      double min_rt = std::numeric_limits<double>::max();
      double max_rt = 0;
      for (auto& l : m)
      {
        double tmp_min = l.getMassTrace()->getConvexhull().getBoundingBox().minX();
        double tmp_max = l.getMassTrace()->getConvexhull().getBoundingBox().maxX();

        if (tmp_min < min_rt)
          min_rt = tmp_min;
        if (tmp_max > max_rt)
          max_rt = tmp_max;
      }

      // fwhm range
//      m.getFwhmRange()

      out << charges << "\t"
      << to_string(m.getFwhmRange().first) << "\t"
      << to_string(m.getFwhmRange().second) << "\t"
      << to_string(min_rt) << "\t"
      << to_string(max_rt) << "\n";
    }

    out.close();
  }

  void FLASHDeconvQuant::run(std::vector<MassTrace> &input_mtraces, FeatureMap &output_featmap)
  {
    // *********************************************************** //
    // Step 1 deconvolute mass traces
    // *********************************************************** //
    vector<LogMassTrace> log_mtraces;
    logTransformMassTraces_(input_mtraces, log_mtraces);
    setFilters_();
    setAveragineModel_();
    if (Constants::C13C12_MASSDIFF_U * max_nr_traces_ / charge_lower_bound_ < local_mz_range_ )
    {
      local_mz_range_ = Constants::C13C12_MASSDIFF_U * max_nr_traces_ / charge_lower_bound_;
    }

    std::vector<FeatureGroup> features;
    features.reserve(log_mtraces.size());
    buildMassTraceGroups_(log_mtraces, features);
    features.shrink_to_fit();

    // *********************************************************** //
    // Step 2 mass artifact removal & post processing...
    // *********************************************************** //
    refineFeatureGroups_(features);
    writeFeatureGroupsInFile(features);

    // *********************************************************** //
    // Step 3 clustering features
    // *********************************************************** //
    std::vector<std::vector<Size>> shared_m_traces(input_mtraces.size(), std::vector<Size>());
    clusterFeatureGroups_(features, shared_m_traces);
  }

  void FLASHDeconvQuant::logTransformMassTraces_(std::vector<MassTrace> &input_mtraces, std::vector<LogMassTrace> &log_mtraces)
  {
    // shortest mass trace length? -> update local_rt_range_;
    double shortest_mt_length = local_rt_range_;
    double min_mz = std::numeric_limits<double>::max();
    double max_mz = 0;
    Size index = 0;

    log_mtraces.reserve(input_mtraces.size());
    for (auto iter=input_mtraces.begin(); iter!=input_mtraces.end(); ++iter, ++index)
    {
      LogMassTrace tmp_lmt(*iter);
      tmp_lmt.setTraceIndex(index);
      log_mtraces.push_back(tmp_lmt);
      if (iter->getFWHM() < shortest_mt_length)
      {
        shortest_mt_length = iter->getFWHM();
      }
      if (iter->getCentroidMZ() > max_mz)
      {
        max_mz = iter->getCentroidMZ();
      }
      else if(iter->getCentroidMZ() < min_mz)
      {
        min_mz = iter->getCentroidMZ();
      }
    }
    local_rt_range_ = shortest_mt_length;
    upper_bound_mz_ = max_mz;
    lower_bound_mz_ = min_mz;

    // sort input mass traces in RT
    std::sort(log_mtraces.begin(), log_mtraces.end(), CmpLogMassTraceByRT());
  }



  // TODO: replace this with FLASHDeconv's
  void FLASHDeconvQuant::setFilters_()
  {
    filter_.clear();
    harmonic_filter_matrix_.clear();
    for (int i = 0; i < charge_range_; i++)
    {
      filter_.push_back(log(1.0 / (i + charge_lower_bound_)));
    }

    harmonic_filter_matrix_.resize(harmonic_charges_.size(), charge_range_);

    for (Size k = 0; k < harmonic_charges_.size(); k++)
    {
      int hc = harmonic_charges_[k];
      int n = hc / 2;

      for (int i = 0; i < charge_range_; i++)
      {
        harmonic_filter_matrix_.setValue(k, i, log(1.0 / (1.0 * n / hc + (i + charge_lower_bound_))));
      }
    }
  }

  void FLASHDeconvQuant::setAveragineModel_()
  {
    auto generator = new CoarseIsotopePatternGenerator();
    auto maxIso = generator->estimateFromPeptideWeight(max_mass_);
    maxIso.trimRight(0.01 * maxIso.getMostAbundant().getIntensity());

    generator->setMaxIsotope(maxIso.size());
    iso_model_ = PrecalculatedAveragine(50, max_mass_, 25, generator);
    iso_model_.setMaxIsotopeIndex(maxIso.size() - 1);

    max_nr_traces_ = iso_model_.getMaxIsotopeIndex();
  }

  double FLASHDeconvQuant::getBinValue_(const Size &bin, const double &min_value, const double &bin_width) const
  {
    return min_value + bin / bin_width;
  }

  Size FLASHDeconvQuant::getBinNumber_(const double &value, const double &min_value, const double &bin_width) const
  {
    if (value < min_value)
    {
      return 0;
    }
    return (Size) (((value - min_value) * bin_width) + .5);
  }

  // From log mz to mz bins. To reduce edge effect, its adjacent bin is set when a bin is set.
  void FLASHDeconvQuant::updateMzBins_(std::vector<LogMassTrace*> &local_traces, const Size &bin_number,
                                       const double& mz_bin_min, std::vector<float> &mz_bin_intensities)
  {
    mz_bins_for_edge_effect_ = boost::dynamic_bitset<>(bin_number);
    mz_bins_ = boost::dynamic_bitset<>(bin_number);
    //std::fill(mz_bin_intensities.begin(), mz_bin_intensities.end(), .0);
    for (auto &trace : local_traces)
    {
      Size bi = getBinNumber_(trace->getLogCentroidMz(), mz_bin_min, mz_bin_width_);
      if (bi >= bin_number)
      {
        continue;
      }
      mz_bins_.set(bi);
      mz_bins_for_edge_effect_.set(bi);
      mz_bin_intensities[bi] += trace->getIntensity();
    }
    for (auto &trace : local_traces)
    {
      Size bi = getBinNumber_(trace->getLogCentroidMz(), mz_bin_min, mz_bin_width_);
      double delta = (trace->getLogCentroidMz() - getBinValue_(bi, mz_bin_min, mz_bin_width_));

      if (delta > 0)
      {
        if (bi < bin_number - 1
            && !mz_bins_for_edge_effect_[bi + 1]
            )
        {
          mz_bins_for_edge_effect_.set(bi + 1);
          mz_bin_intensities[bi + 1] += trace->getIntensity();
        }
      }
      else if (delta < 0)
      {
        if (bi > 0
            && !mz_bins_for_edge_effect_[bi - 1]
            )
        {
          mz_bins_for_edge_effect_.set(bi - 1);
          mz_bin_intensities[bi - 1] += trace->getIntensity();
        }
      }
    }
  }

  //take the mass bins from previous overlapping spectra and put them in the candidate mass bins.
  void FLASHDeconvQuant::unionPrevMassBins_()
  {
    if (mass_bins_.empty())
    {
      return;
    }
    long shift = (long) (round((mass_bin_min_value_) * mz_bin_width_));

    for (Size i = 0; i < prev_mass_bin_vector_.size(); i++)
    {
      auto &pmb = prev_mass_bin_vector_[i];
      //getBinNumber_(pg.getMonoMass() + mass_delta, 0, bin_width_[ms_level_ - 1]);
      if (pmb.empty())
      {
        continue;
      }

      for (Size &index : pmb)
      {
        long j = (long) index - shift;
        if (j < 1)
        {
          continue;
        }
        //std::cout << index << " " << shift << " " << j << " " << mass_bins_.size() << std::endl;
        if ((Size) j >= mass_bins_.size()-1)
        {
          break;
        }
        mass_bins_[j-1] = true;
        mass_bins_[j] = true;
        mass_bins_[j+1] = true;
      }
    }
    for (Size &index : target_mass_bins_) {
      long j = (long) index - shift;
      if (j < 1) {
        continue;
      }
      if ((Size) j >= mass_bins_.size()-1) {
        break;
      }
      mass_bins_[j-1] = true;
      mass_bins_[j] = true;
      mass_bins_[j+1] = true;
    }
  }

  //update mass bins which will be used to select peaks in the input spectrum...
  Matrix<int> FLASHDeconvQuant::updateMassBins_(const std::vector<float> &mz_intensities)
  {
    auto mass_intensities = std::vector<float>(mass_bins_.size(), 0);
    updateCandidateMassBins_(mass_intensities, mz_intensities);
    auto per_mass_abs_charge_ranges = filterMassBins_(mass_intensities);

    return per_mass_abs_charge_ranges;
  }

  //With mass_bins_, select peaks from the same mass in the original input spectrum
  void FLASHDeconvQuant::getCandidatePeakGroups_(const std::vector<LogMassTrace*> &log_mtraces,
                                                 const Matrix<int> &per_mass_abs_charge_ranges,
                                                 std::vector<FeatureGroup> &fgroup)
  {
    //const int max_missing_isotope = 3;
    Size mass_bin_size = mass_bins_.size();
    int log_mz_peak_size = (int) log_mtraces.size();
    auto current_peak_index = std::vector<int>(charge_range_, 0);
    fgroup.reserve(mass_bins_.count());
    Size mass_bin_index = mass_bins_.find_first();
    auto peak_bin_numbers = std::vector<Size>(log_mz_peak_size);

    for (int i = 0; i < log_mz_peak_size; i++)
    {
      peak_bin_numbers[i] = getBinNumber_(log_mtraces[i]->getLogCentroidMz(), mz_bin_min_value_, mz_bin_width_);
    }

    while (mass_bin_index != mass_bins_.npos)
    {
      double log_m = getBinValue_(mass_bin_index, mass_bin_min_value_, mz_bin_width_);
      double mass = exp(log_m);
      FeatureGroup fg(charge_lower_bound_,
                   per_mass_abs_charge_ranges.getValue(1, mass_bin_index) + lower_bound_mz_);

      fg.reserve(charge_range_ * 30);
      Size right_index = iso_model_.getLeftCountFromApex(mass);
      Size left_index = iso_model_.getRightCountFromApex(mass);

      for (int j = per_mass_abs_charge_ranges.getValue(0, mass_bin_index);
           j <= per_mass_abs_charge_ranges.getValue(1, mass_bin_index);
           j++)
      {
        int &bin_offset = bin_offsets_[j];
        int b_index = mass_bin_index - bin_offset;

        double max_intensity = -1.0;
        int abs_charge = j + charge_lower_bound_;
        int &cpi = current_peak_index[j];
        int max_peak_index = -1;

        while (cpi < log_mz_peak_size - 1)
        {
          if (peak_bin_numbers[cpi] == b_index)
          {
            double intensity = log_mtraces[cpi]->getIntensity();
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

        const double mz = log_mtraces[max_peak_index]->getCentroidMz();
        const double iso_delta = Constants::ISOTOPE_MASSDIFF_55K_U / (abs_charge);
        double mz_delta = mz_tolerance_ * mz * 2; //

        double charge_snr = .0;

        int candidate_i = 0;
        // int peakcntr = 0;
        for (int peak_index = max_peak_index; peak_index < log_mz_peak_size; peak_index++)
        {
          const double observed_mz = log_mtraces[peak_index]->getCentroidMz();
          const double intensity = log_mtraces[peak_index]->getIntensity();
          //observedMz = mz + isof * i * d - d * mzDelta;
          double mz_diff = observed_mz - mz;

          int tmp_i = (int) (.5 + mz_diff / iso_delta);

          if (tmp_i > (int) right_index)
          {
            break;
          }

          if (abs(mz_diff - tmp_i * iso_delta) >= mz_delta) // noise   max_intensity  vs   intensity
          {
            charge_snr += intensity * intensity;
            //peakcntr++;
          }
          else
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size)
            {
              LogMassTrace p(*log_mtraces[peak_index]);
              p.setCharge(abs_charge);
              p.setIsotopeIndex(tmp_i);
              fg.push_back(p);
            }
            candidate_i = tmp_i;
          }
        }

        candidate_i = 0;

        for (int peak_index = max_peak_index - 1; peak_index >= 0; peak_index--)
        {
          const double observed_mz = log_mtraces[peak_index]->getCentroidMz();
          const double intensity = log_mtraces[peak_index]->getIntensity();

          //observedMz = mz + isof * i * d - d * mzDelta;
          double mz_diff = mz - observed_mz;
          int tmp_i = (int) (.5 + mz_diff / iso_delta);

          if (tmp_i > (int) left_index)
          {
            break;
          }

          //if (tmp_i - candidate_i > max_missing_isotope)
          //{
          //  break;
          //}

          if (abs(mz_diff - tmp_i * iso_delta) >= mz_delta)
          {
            charge_snr += intensity * intensity;
          }
          else
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size)
            {
              LogMassTrace p(*log_mtraces[peak_index]);
              p.setCharge(abs_charge);
              p.setIsotopeIndex(tmp_i);
              fg.push_back(p);
            }
            candidate_i = tmp_i;
          }
        }
        // TODO : revive charge snr?
//        if (charge_snr > 0)
//        {
//          pg.setChargeSNR(abs_charge, charge_snr);
//        }
      }

      if (!fg.empty())
      {

        double max_intensity = -1.0;
        //double sum_intensity = .0;
        double t_mass = .0;
        auto new_peaks = std::vector<LogMassTrace>();
        new_peaks.reserve(fg.size());

        //std::unordered_map<int, double> isotope_intensity;

        for (auto &p : fg)
        {
          //sum_intensity += p.intensity;
          // t_mass += p.intensity * (p.getUnchargedMass() - p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U);
          if (max_intensity < p.getIntensity())
          {
            max_intensity = p.getIntensity();
            t_mass = p.getUnchargedMass();
          }
        }

        double iso_tolerance = mz_tolerance_ * t_mass;
        int min_off = 10000;
        std::vector<double> noise_power(charge_upper_bound_ + 1, .0);
        for (auto &p : fg)
        {
          p.setIsotopeIndex(round((p.getUnchargedMass() - t_mass) / Constants::ISOTOPE_MASSDIFF_55K_U));
          if (abs(t_mass - p.getUnchargedMass() + Constants::ISOTOPE_MASSDIFF_55K_U * p.getIsotopeIndex()) >
              iso_tolerance)
          {
            noise_power[p.getCharge()] += p.getIntensity() * p.getIntensity();
            continue;
          }
          new_peaks.push_back(p);
          min_off = min_off > p.getIsotopeIndex() ? p.getIsotopeIndex() : min_off;
        }
        for (int abs_charge = 0; abs_charge < noise_power.size(); abs_charge++)
        {
          double np = noise_power[abs_charge];
          if (np <= 0)
          {
            continue;
          }
          // TODO : revive charge snr?
//          pg.setChargeSNR(abs_charge, fg.getChargeSNR(abs_charge) + np);
        }

        fg.swap(new_peaks);

        for (auto &p : fg)
        {
          p.setIsotopeIndex(p.getIsotopeIndex()-min_off);
        }
        fg.updateMassesAndIntensity();

        fgroup.push_back(fg); //
      }
      mass_bin_index = mass_bins_.find_next(mass_bin_index);
    }
  }

  //Find candidate mass bins from the current spectrum. The runtime of FLASHDeconv is deteremined by this function..
  void FLASHDeconvQuant::updateCandidateMassBins_(std::vector<float> &mass_intensitites,
                                                  const std::vector<float> &mz_intensities)
  {
    int charge_range = charge_upper_bound_ - charge_lower_bound_ + 1;
    int h_charge_size = (int) harmonic_charges_.size();
    int min_peak_cntr = min_nr_mtraces_;
    long bin_end = (long) mass_bins_.size();
    // how many peaks of continuous charges per mass
    auto support_peak_count = std::vector<int>(mass_bins_.size(), 0);

    Size mz_bin_index = mz_bins_.find_first();

    // to calculate continuous charges, the previous charge value per mass should be stored
    auto prev_charges = std::vector<int>(mass_bins_.size(), charge_range + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prev_intensities = std::vector<float>(mass_bins_.size(), 1.0f);

    double bin_width = mz_bin_width_;

    // intensity change ratio should not exceed the factor.
    const float factor = 5.0;
    const float hfactor = 1.1;
    const int low_charge = 6;
    while (mz_bin_index != mz_bins_.npos)
    {
      float intensity = mz_intensities[mz_bin_index];
      double mz = -1.0, log_mz = 0;
      log_mz = getBinValue_(mz_bin_index, mz_bin_min_value_, bin_width);
      mz = exp(log_mz);

      // scan through charges
      for (int j = 0; j < charge_range; j++) //  increasing.
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
        int abs_charge = (j + charge_lower_bound_);
        float &prev_intensity = prev_intensities[mass_bin_index];
        int &prev_charge = prev_charges[mass_bin_index];

        bool pass_first_check = false;
        if (abs_charge <= low_charge)
        { // for low charges
          double diff = Constants::C13C12_MASSDIFF_U / abs_charge / mz;
          Size next_iso_bin = getBinNumber_(log_mz + diff, mz_bin_min_value_, bin_width);

          if (next_iso_bin < mz_bins_for_edge_effect_.size() && mz_bins_for_edge_effect_[next_iso_bin])
          {
            if (mz_intensities[next_iso_bin] < intensity)
            {
              pass_first_check = true;
            }
          }
        }
        else
        {
          // intensity of previous charge
          // intensity ratio between current and previous charges
          float intensity_ratio = intensity / prev_intensity;
          intensity_ratio = intensity_ratio < 1 ? 1.0f / intensity_ratio : intensity_ratio;

          // check if peaks of continuous charges are present
          bool charge_not_continous = prev_charge - j != 1;
          // if charge not continous or intensity ratio is too high reset continuousChargePeakPairCount
          if (charge_not_continous || intensity_ratio > factor)
          {
            spc = 0;
            //mass_intensitites[mass_bin_index] -= intensity; //
          }
          else
          {
            pass_first_check = true;
          }
        }
        if (pass_first_check)
        { // check harmonic artifacts
          float max_intensity = intensity;
          float min_intensity = prev_intensity;
          if (abs_charge <= low_charge && prev_intensity <= 1.0)
          {
            max_intensity = intensity * factor; // no consecutive intensity. Thus, allow factor
            min_intensity = intensity / factor;
          }
          else if (min_intensity > max_intensity)
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
            if (spc == 0 && abs_charge > low_charge)
            {
              mass_intensitites[mass_bin_index] += prev_intensity;
            }
            mass_intensitites[mass_bin_index] += intensity;
            if (++spc >= min_peak_cntr || abs_charge <= low_charge) //
            {
              mass_bins_[mass_bin_index] = true;
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
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }
  }

  // Subfunction of updateMassBins_. If a peak corresponds to multiple masses, only one mass is selected based on intensities..
  Matrix<int> FLASHDeconvQuant::filterMassBins_(const std::vector<float> &mass_intensities)
  {
    //int chargeRange = param.currentChargeRange;
    int charge_range = charge_upper_bound_ - charge_lower_bound_ + 1;
    double bin_width = mz_bin_width_;
    Matrix<int> abs_charge_ranges(2, mass_bins_.size(), INT_MAX);
    for (int i = 0; i < mass_bins_.size(); i++)
    {
      abs_charge_ranges.setValue(1, i, INT_MIN);
    }
    Size mz_bin_index = mz_bins_.find_first();
    long bin_size = (long) mass_bins_.size();

    //massBinsForThisSpectrum = boost::dynamic_bitset<>(mass_bins_.size());

    auto to_skip = mass_bins_.flip();
    mass_bins_.reset();

    while (mz_bin_index != mz_bins_.npos)
    {
      long max_index = -1;
      float max_intensity = -1e11;
      int max_intensity_abs_charge = 0;

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

        if (max_intensity < t)
        {
          bool artifact = false;

          double original_log_mass = getBinValue_(mass_bin_index, mass_bin_min_value_, bin_width);
          double mass = exp(original_log_mass);
          double diff = Constants::C13C12_MASSDIFF_U / mass;
          for (int iso_off = -2; iso_off <= 2 && !artifact; ++iso_off)
          {
            double log_mass = original_log_mass + diff * iso_off;
            if (log_mass < 1)
            {
              continue;
            }
            for (int h = 2; h <= 6 && !artifact; h++)
            {
              for (int f = -1; f <= 1 && !artifact; f += 2) //
              {
                double hmass = log_mass - log(h) * f;
                Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_, bin_width);
                if (hmass_index > 0 && hmass_index < mass_bins_.size() - 1)
                {
                  //for (int off = 0; off <= 0 && !artifact; off++)
                  //{
                  if (mass_intensities[hmass_index] >= t)
                  {
                    artifact = true;
                    break;
                  }
                  //}
                }
              }
            }
            //   max_intensity_abs_charge off by one here
            if (!artifact)
            {
              int abs_charge = (j + charge_lower_bound_);
              for (int coff = 1; coff <= 2 && !artifact; coff++)
              {
                for (int f = -1; f <= 1 && !artifact; f += 2)
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
            max_intensity_abs_charge = j;
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
        abs_charge_ranges
            .setValue(0, max_index, std::min(abs_charge_ranges.getValue(0, max_index), max_intensity_abs_charge));
        abs_charge_ranges
            .setValue(1, max_index, std::max(abs_charge_ranges.getValue(1, max_index), max_intensity_abs_charge));
        mass_bins_[max_index] = true;
      }
      mz_bin_index = mz_bins_.find_next(mz_bin_index);
    }

    return abs_charge_ranges;
  }

  bool FLASHDeconvQuant::checkChargeDistribution_(const std::vector<double> &per_charge_intensity) {

    double max_per_charge_intensity = .0;
    double sum_charge_intensity = .0;
    int cntr = 0;
    int non_zero_start = -1, non_zero_end = 0;
    for (int i = 0; i < charge_range_; i++) {
      if (per_charge_intensity[i] > 0) {
        sum_charge_intensity += per_charge_intensity[i];
        cntr++;
        max_per_charge_intensity = std::max(max_per_charge_intensity, per_charge_intensity[i]);
        if (non_zero_start < 0) {
          non_zero_start = i;
        }
        non_zero_end = i;
      }
    }
    if (cntr == 0) {
      return false;
    }

    int prev_charge = non_zero_start;
    int n_r = 0;
    float factor = 1;
    double int_threshold = (sum_charge_intensity / cntr) * factor;
    for (int k = prev_charge + 1; k <= non_zero_end; k++) {
      if (per_charge_intensity[k] < int_threshold) {
        continue;
      }

      if (k - prev_charge == 1) {
        n_r++;
      }
      //spc >= a_charge/2
      if (n_r >= min_nr_mtraces_ - 1 || n_r >= cntr / 2)//
      {
        return true;
      }
      prev_charge = k;
    }

    return false;
  }

  std::vector<int> FLASHDeconvQuant::calculatePerChargeIsotopeIntensity_(
      std::vector<double> &per_isotope_intensity,
      std::vector<double> &per_charge_intensity,
      const int max_isotope_count,
      FeatureGroup &fg) const {
    int min_pg_charge = INT_MAX;
    int max_pg_charge = INT_MIN;
    //    double maxIntensity = -1;
    int max_intensity_charge_index = -1;
    //    double maxIntensity2 = -1;
    int max_intensity_iso_index = -1;

    for (auto &p : fg) {
      if (p.getIsotopeIndex() < 0 || p.getIsotopeIndex() >= max_isotope_count) {
        continue;
      }
      min_pg_charge = std::min(min_pg_charge, p.getCharge());
      max_pg_charge = std::max(max_pg_charge, p.getCharge());

      int index = p.getCharge() - charge_lower_bound_;//current_min_charge_;
      per_isotope_intensity[p.getIsotopeIndex()] += p.getIntensity();
      per_charge_intensity[index] += p.getIntensity();
    }
    fg.setChargeRange(min_pg_charge, max_pg_charge);

    return std::vector<int>{max_intensity_charge_index, max_intensity_iso_index};
  }

  double FLASHDeconvQuant::getChargeFitScore_(const std::vector<double> &per_charge_intensity) const {
    double max_per_charge_intensity = .0;
    double summed_intensity = .0;
    int max_index = -1;
    int first_index = -1;
    int last_index = charge_range_ - 1;

    for (int i = 0; i < charge_range_; i++) {
      summed_intensity += per_charge_intensity[i];
      if (per_charge_intensity[i] <= 0) {
        if (first_index < 0) {
          first_index = i;
        }
        last_index = i;
      }

      if (max_per_charge_intensity > per_charge_intensity[i]) {
        continue;
      }
      max_per_charge_intensity = per_charge_intensity[i];
      max_index = i;
    }
    first_index = first_index < 0 ? 0 : first_index;

    double p = .0;
    for (int i = max_index; i < last_index - 1; i++) {
      double diff = per_charge_intensity[i + 1] - per_charge_intensity[i];
      double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i + 1]);
      if (diff <= 0 && ratio < 5.0) {
        continue;
      }
      p += abs(diff);
    }

    for (int i = max_index; i > first_index; i--) {
      double diff = per_charge_intensity[i - 1] - per_charge_intensity[i];
      //      double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i - 1]);

      if (diff <= 0) {
        continue;
      }
      p += abs(diff);
    }
    return std::max(.0, 1.0 - p / summed_intensity);
  }

  double FLASHDeconvQuant::getCosine_(const std::vector<double> &a,
                                          const int &a_start,
                                          const int &a_end,
                                          const IsotopeDistribution &b,
                                          const int &b_size,
                                          const int offset) const {
    double n = .0, a_norm = .0;
    //int c = 0;
    for (int j = a_start; j <= a_end; j++) {
      int i = j - offset;
      if (i < 0 || i >= b_size || b[i].getIntensity() <= 0) {
        continue;
      }
      a_norm += a[j] * a[j];
      n += a[j] * b[i].getIntensity(); //
    }
    if (a_norm <= 0) {
      return 0;
    }
    return n / sqrt(a_norm);
  }

  double FLASHDeconvQuant::getShapeDiff_(const std::vector<double> &a,
                                             const int &a_start,
                                             const int &a_end,
                                             const IsotopeDistribution &b,
                                             const int &b_size,
                                             const int max_b_index,
                                             const int offset) const{
    double a_norm = .0, b_norm = .0;
    double diff = 0;

    for (int j = a_start; j <= a_end; j++) {
      int i = j - offset;
      if (i < 0 || i >= b_size || b[i].getIntensity() <= 0) {
        continue;
      }

      if (i < max_b_index - 2 || i > max_b_index + 2) {
        continue;
      }
      a_norm += a[j] * a[j];
      b_norm += b[i].getIntensity() * b[i].getIntensity();
    }

    if (a_norm <= 0 || b_norm <= 0) {
      return -1;
    }
    a_norm = sqrt(a_norm);
    b_norm = sqrt(b_norm);

    for (int j = a_start; j <= a_end; j++) {
      int i = j - offset;
      if (i < 0 || i >= b_size || b[i].getIntensity() <= 0) {
        continue;
      }
      double t_diff = (a[j] / a_norm - b[i].getIntensity() / b_norm);
      diff += (t_diff * t_diff);
    }

    return diff;
  }

  double FLASHDeconvQuant::getIsotopeCosineAndDetermineIsotopeIndex(const double mono_mass,
                                                                        const std::vector<double> &per_isotope_intensities,
                                                                        int &offset,
                                                                        const PrecalculatedAveragine &avg) const {
    auto iso = avg.get(mono_mass);
    //double iso_norm = avg.getNorm(mono_mass);

    int iso_size = (int) iso.size();
    int apex_index = avg.getApexIndex(mono_mass);
    int iso_range = avg.getRightCountFromApex(mono_mass) + avg.getLeftCountFromApex(mono_mass);

    offset = 0;
    double min_diff = -1;
    //int isotope_length = 0;
    int max_isotope_index = 0, min_isotope_index = -1;

    for (int i = 0; i < avg.getMaxIsotopeIndex(); i++) {
      if (per_isotope_intensities[i] <= 0) {
        continue;
      }
      //isotope_length++;
      max_isotope_index = i;
      if (min_isotope_index < 0) {
        min_isotope_index = i;
      }
    }

    //        double norm = .0;
    //       for (int j = min_isotope_index; j <= max_isotope_index; j++) {
    //           norm += per_isotope_intensities[j] * per_isotope_intensities[j];
    //       }

    for (int tmp_offset = -apex_index - 1; tmp_offset <= -apex_index + max_isotope_index + 1; tmp_offset++) {
      double tmp_diff = getCosine_(per_isotope_intensities,
                                   min_isotope_index,
                                   max_isotope_index,
                                   iso,
                                   iso_size, //apex_index,
                                   tmp_offset);
      //if (tmp_diff < 0) {
      //    continue;
      //}

      if (min_diff < tmp_diff) {//min_diff < 0 ||
        min_diff = tmp_diff;
        offset = tmp_offset;
      }
    }

    min_diff = -1;
    int final_offset = offset;
    if (iso_range / 4 > 0) {
      for (int tmp_offset = offset - iso_range / 4; tmp_offset <= offset + iso_range / 4; tmp_offset++) {
        double tmp_diff = getShapeDiff_(per_isotope_intensities,
                                        min_isotope_index,
                                        max_isotope_index,
                                        iso,
                                        iso_size, apex_index,
                                        tmp_offset);
        if (tmp_diff < 0) {
          continue;
        }

        if (min_diff < 0 || min_diff > tmp_diff) {//||
          min_diff = tmp_diff;
          final_offset = tmp_offset;
        }
      }
    }

    return getCosine_(per_isotope_intensities,
                      min_isotope_index,
                      max_isotope_index,
                      iso,
                      iso_size,
                      final_offset);
  }

  float FLASHDeconvQuant::getAvgPPMError_(FeatureGroup pg) const{
    std::vector<float> diffs;
    std::vector<std::vector<float>> per_isotope_masses;
    int isotope_end_index = 0;

    for (auto &p : pg) {
      isotope_end_index = isotope_end_index < p.getIsotopeIndex() ? p.getIsotopeIndex() : isotope_end_index;
    }
    per_isotope_masses = std::vector<std::vector<float>>(isotope_end_index + 1, std::vector<float>());
    for (auto &p : pg) {
      per_isotope_masses[p.getIsotopeIndex()].push_back(p.getUnchargedMass());
    }
    diffs.reserve(pg.size());
    for (int i = 0; i < per_isotope_masses.size(); i++) {
      auto &v = per_isotope_masses[i];
      Size n = v.size();
      //if (n < 2)
      //{
      //  continue;
      //}
      double average = n >= 2 ?
                       accumulate(v.begin(), v.end(), 0.0) / n :
                       pg.getMonoisotopicMass() + i * Constants::ISOTOPE_MASSDIFF_55K_U;  //
      for (float &t:v) {
        diffs.push_back(pow(1e6 * (t - average) / average, 2.0));
      }
    }
    Size n = diffs.size();
    return n == 0 ? .0 : sqrt(accumulate(diffs.begin(), diffs.end(), 0.0) / n);

  }



  void FLASHDeconvQuant::scoreAndFilterPeakGroups_(std::vector<FeatureGroup> &local_fgroup) {
    std::vector<FeatureGroup> filtered_peak_groups;
    filtered_peak_groups.reserve(local_fgroup.size());

    for (auto &feature_group : local_fgroup) {
      auto per_isotope_intensities = std::vector<double>(iso_model_.getMaxIsotopeIndex(), 0);
      auto per_abs_charge_intensities = std::vector<double>(charge_range_, 0);

      auto indices = calculatePerChargeIsotopeIntensity_(
          per_isotope_intensities, per_abs_charge_intensities,
          iso_model_.getMaxIsotopeIndex(), feature_group);

      // TODO : remove? does charge distribution makes big difference?
      double cs = getChargeFitScore_(per_abs_charge_intensities);
      feature_group.setChargeScore(cs);

      // if (ms_level_ == 1)
      {
        bool is_charge_well_distributed = checkChargeDistribution_(per_abs_charge_intensities);
        //double tmp = getChargeFitScore_(per_abs_charge_intensities, charge_range);

        if (!is_charge_well_distributed) {
            continue;
        }
      }
      /// TODO_ends_here

      int offset = 0;
      double cos = getIsotopeCosineAndDetermineIsotopeIndex(feature_group[0].getUnchargedMass(),
                                                            per_isotope_intensities,
                                                            offset, iso_model_);
      feature_group.setIsotopeCosine(cos);

      if (feature_group.empty() ||
          (feature_group.getIsotopeCosine() <= min_isotope_cosine_))// (msLevel <= 1 ? param.minIsotopeCosineSpec : param.minIsotopeCosineSpec2)))
      {
          continue;
      }

      feature_group.updateMassesAndIntensity(offset, iso_model_.getMaxIsotopeIndex());
      if (feature_group.getMonoisotopicMass() < min_mass_ || feature_group.getMonoisotopicMass() > max_mass_) {
        continue;
      }
      auto iso_dist = iso_model_.get(feature_group.getMonoisotopicMass());
      int iso_size = (int) iso_dist.size();
      float total_noise = .0;
      float total_signal = .0;
      //auto perChargeMaxIntensity = std::vector<double>(chargeRange);

      auto current_charge_range = feature_group.getChargeRange();
      for (int abs_charge = std::get<0>(current_charge_range);
           abs_charge <= std::get<1>(current_charge_range);
           ++abs_charge) {
        int j = abs_charge - charge_lower_bound_;//current_min_charge_;
        if (per_abs_charge_intensities[j] <= 0) {
          continue;
        }
        auto current_per_isotope_intensities = std::vector<double>(iso_model_.getMaxIsotopeIndex(), 0);

        int min_isotope_index = iso_model_.getMaxIsotopeIndex();
        int max_isotope_index = 0;

        double max_intensity = .0;
        //double sumIntensity = .0;
        //double summed_intensity_squares = .0;

        for (auto &peak: feature_group) {
          if (peak.getCharge() != abs_charge) {
            continue;
          }

          if (peak.getIsotopeIndex() > iso_size) {
            continue;
          }

          current_per_isotope_intensities[peak.getIsotopeIndex()] += peak.getIntensity();
          //sumIntensity += p.intensity;
          min_isotope_index = min_isotope_index < peak.getIsotopeIndex() ? min_isotope_index : peak.getIsotopeIndex();
          max_isotope_index = max_isotope_index < peak.getIsotopeIndex() ? peak.getIsotopeIndex() : max_isotope_index;

          //min_mz = min_mz < p.mz ? min_mz : p.mz;
          // max_mz = max_mz > p.mz ? max_mz : p.mz;
          if (max_intensity < peak.getIntensity()) {
            max_intensity = peak.getIntensity();
            // perChargeMaxIntensity[j] = maxIntensity;
          }
          // sp += p.intensity * p.intensity;
        }
        if (max_intensity <= 0) {
          continue;
        }

        for (int k = min_isotope_index; k <= max_isotope_index; ++k) {
          if (k > iso_size) {
            break;
          }
          //summed_intensity_squares += current_per_isotope_intensities[k] * current_per_isotope_intensities[k];
        }
        //double norm = .0;
        //for (int j = min_isotope_index; j <= max_isotope_index; j++) {
        //    norm += current_per_isotope_intensities[j] * current_per_isotope_intensities[j];
        //}

        double cos_score = getCosine_(current_per_isotope_intensities,
                                      min_isotope_index,
                                      max_isotope_index,
                                      iso_dist,
                                      iso_size,
            // norm,
            //1,
                                      0);

        // double cos_score_squared = cos_score * cos_score;

        feature_group.setChargeIsotopeCosine(abs_charge, cos_score);
        feature_group.setChargeIntensity(abs_charge, per_abs_charge_intensities[j]);

        // double noise = (1 - cos_score_squared) * summed_intensity_squares + peak_group.getChargeSNR(abs_charge) + 1;
        // double signal = cos_score_squared * summed_intensity_squares + 1;

        //peak_group.setChargeSNR(abs_charge, noise / signal);

      }

      //peak_group.setSNR(total_signal / total_noise);

      feature_group.setAvgPPMError(getAvgPPMError_(feature_group));

      //if (ms_level_==1 &&  peak_group.getSNR() < 2.0) // tmp
      //{//
      //     continue;
      //}

      // TODO revive here?
//      feature_group.updateSNR();
//      feature_group.setQScore(-10000);

//      for (int abs_charge = std::get<0>(current_charge_range);
//           abs_charge <= std::get<1>(current_charge_range);
//           abs_charge++) {
//        if (feature_group.getChargeIntensity(abs_charge) <= 0) {
//          continue;
//        }
//        int j = abs_charge - charge_lower_bound_;//current_min_charge_;
//
//        double q_score = QScore::getQScore(&peak_group, abs_charge);
//
//        if (q_score <= feature_group.getQScore()) {
//          continue;
//        }
//        feature_group.setRepAbsCharge(abs_charge);
//        feature_group.setQScore(q_score);
//      }
//      if (ms_level_ == 1 && peak_group.getRepAbsCharge() < min_abs_charge_) {
//        continue;
//      }

//      auto max_q_score_mz_range = feature_group.getMzRange(feature_group.getRepAbsCharge());
//      if (std::get<0>(max_q_score_mz_range) > std::get<1>(max_q_score_mz_range)) {
//        continue;
//      }
//      feature_group.setMaxQScoreMzRange(std::get<0>(max_q_score_mz_range), std::get<1>(max_q_score_mz_range));
      filtered_peak_groups.push_back(feature_group);
    }
    local_fgroup.swap(filtered_peak_groups);

    // TODO : revive here?
//    filterPeakGroupsByIsotopeCosine_(max_c);
  }

  void FLASHDeconvQuant::removeOverlappingPeakGroups_(std::vector<FeatureGroup> &local_fgroup,
                                                      const double tol,
                                                      const int iso_length) const{
    std::vector<FeatureGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(local_fgroup.size());
    sort(local_fgroup.begin(), local_fgroup.end());

    for (Size i = 0; i < local_fgroup.size(); i++) {
      if (i > 0) {
        if (abs(local_fgroup[i - 1].getMonoisotopicMass() - local_fgroup[i].getMonoisotopicMass()) < 1e-3
            &&
            local_fgroup[i - 1].getIntensity() >= local_fgroup[i].getIntensity()) {
          continue;
        }
      }

      filtered_pg_vec.push_back(local_fgroup[i]);
    }
    local_fgroup.swap(filtered_pg_vec);
    std::vector<FeatureGroup>().swap(filtered_pg_vec);
    filtered_pg_vec.reserve(local_fgroup.size());

    for (Size i = 0; i < local_fgroup.size(); i++) {
      bool select = true;
      auto &pg = (local_fgroup)[i];

      bool pass = false;

      if (pass) {
        filtered_pg_vec.push_back(pg);
        continue;
      }

      if (pg.getMonoisotopicMass() <= 0) {
        continue;
      }
      double mass_tolerance = pg.getMonoisotopicMass() * tol * 2;

      int j = i + 1;
      for (int l = 0; l <= iso_length; l++) {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j < local_fgroup.size(); j++) {
          auto &pgo = (local_fgroup)[j];

          if (l != 0 && pgo.getMonoisotopicMass() - pg.getMonoisotopicMass() < off - mass_tolerance) {
            continue;
          }

          if (pgo.getMonoisotopicMass() - pg.getMonoisotopicMass() > off + mass_tolerance) {
            break;
          }
          select &= pg.getIsotopeCosine() >= pgo.getIsotopeCosine();
          if (!select) {
            break;
          }
        }
        if (!select) {
          break;
        }
      }

      if (!select) {
        continue;
      }

      j = i - 1;
      for (int l = 0; l <= iso_length; l++) {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j >= 0; j--) {
          auto &pgo = (local_fgroup)[j];

          if (l != 0 && pg.getMonoisotopicMass() - pgo.getMonoisotopicMass() < off - mass_tolerance) {
            continue;
          }

          if (pg.getMonoisotopicMass() - pgo.getMonoisotopicMass() > off + mass_tolerance) {
            break;
          }
          select &= pg.getIsotopeCosine() >= pgo.getIsotopeCosine();
          if (!select) {
            break;
          }
        }
      }
      if (!select) {
        continue;
      }
      filtered_pg_vec.push_back(pg);
    }
    local_fgroup.swap(filtered_pg_vec);
  }

  void FLASHDeconvQuant::getFeatureFromSpectrum_(std::vector<LogMassTrace *> &local_traces,
                                                 std::vector<FeatureGroup> &local_fgroup)
  {
    /// setting variables for binning (from FLASHDeconv)
    int current_charge_range = charge_upper_bound_ - charge_lower_bound_ + 1;
    int tmp_peak_cntr = current_charge_range - min_nr_mtraces_;

    mass_bin_min_value_ = log(std::max(1.0, min_mass_ - iso_model_.getAverageMassDelta(min_mass_)));
    tmp_peak_cntr = min_nr_mtraces_ - 1;

    /// deconvolution starts here
    mz_bin_min_value_ = local_traces[0]->getLogCentroidMz();
    double mz_bin_max_value = local_traces[local_traces.size()-1]->getLogCentroidMz();
    double mass_bin_max_value = std::min( mz_bin_max_value - filter_[tmp_peak_cntr],
                                          log(max_mass_ + iso_model_.getRightCountFromApex(max_mass_) + 1));

    // TODO : mz_bin_min_value_ remove from method
    Size mass_bin_number = getBinNumber_(mass_bin_max_value, mz_bin_min_value_, mz_bin_width_) + 1;

    // build bin offsets
    bin_offsets_.clear();
    bin_offsets_.reserve(current_charge_range);
    harmonic_bin_offset_matrix_.clear();
    for (int i = 0; i < current_charge_range; i++)
    {
      bin_offsets_.push_back((int) round((mz_bin_min_value_ - filter_[i] - mass_bin_min_value_) * mz_bin_width_));
    }

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
                                  mz_bin_width_));
      }
    }

    Size mz_bin_number = getBinNumber_(mz_bin_max_value, mz_bin_min_value_, mz_bin_width_) + 1;
    auto mz_bin_intensities = std::vector<float>(mz_bin_number, .0f);

    // From log mz to mz bins.
    updateMzBins_(local_traces, mz_bin_number, mz_bin_min_value_, mz_bin_intensities);

    // take the mass bins from previous overlapping spectra and put them in the candidate mass bins. // TODO: remove this? or??
    mass_bins_ = boost::dynamic_bitset<>(mass_bin_number);
    unionPrevMassBins_();
    auto per_mass_abs_charge_ranges = updateMassBins_(mz_bin_intensities);

    getCandidatePeakGroups_(local_traces, per_mass_abs_charge_ranges, local_fgroup);

    // filtering part
    scoreAndFilterPeakGroups_(local_fgroup);
    removeOverlappingPeakGroups_(local_fgroup, mz_tolerance_, 1);

    // TODO : revive rt binning?
//      while (!prev_rt_vector_.empty() &&
//          local_fgroup.getOriginalSpectrum().getRT() - prev_rt_vector_[0] > rt_window_)//
//      {
//        prev_rt_vector_.erase(prev_rt_vector_.begin());
//        prev_mass_bin_vector_.erase(prev_mass_bin_vector_.begin());
//      }

//      std::vector<Size> curr_mass_bin;
//      curr_mass_bin.reserve(local_fgroup.size());
//      for (auto &fg : local_fgroup)//filteredPeakGroups
//      {
//        fg.shrink_to_fit();
//
//        double mass_delta = iso_model_.getAverageMassDelta(fg.getMonoisotopicMass());
//        Size fg_bin = getBinNumber_(log(fg.getMonoisotopicMass() + mass_delta), 0, mz_bin_width_);
//        curr_mass_bin.push_back(fg_bin);
//      }

//      prev_rt_vector_.push_back(deconvoluted_spectrum_.getOriginalSpectrum().getRT());
//      prev_mass_bin_vector_.push_back(curr_mass_bin); //
//      prev_rt_vector_.shrink_to_fit();
//      prev_mass_bin_vector_.shrink_to_fit();
  }

  void FLASHDeconvQuant::addFeatureGroup_(std::vector<FeatureGroup> &features, std::vector<FeatureGroup>& new_fgroups) const
  {
    std::sort(new_fgroups.begin(), new_fgroups.end());

    Size fg_pointer = 0;
    for (auto new_fg_iter = new_fgroups.begin(); new_fg_iter != new_fgroups.end(); ++new_fg_iter)
    {
      // reached end of features. add the rest of new_fgroups to features
      if ( features.size() == fg_pointer )
      {
        features.insert(features.end(), new_fg_iter, new_fgroups.end());
        break;
      }

      // check if similar mass doesn't exist in the vector within tolerance
      if ( features[fg_pointer].getMonoisotopicMass() - new_fg_iter->getMonoisotopicMass() > mass_tolerance_ )
      {
        features.insert(features.begin()+fg_pointer, *new_fg_iter);
        continue;
      }

      // move fg_pointer to fit new_fg
      if ( new_fg_iter->getMonoisotopicMass() - features[fg_pointer].getMonoisotopicMass() > mass_tolerance_)
      {
        for (; std::abs(new_fg_iter->getMonoisotopicMass() - features[fg_pointer].getMonoisotopicMass()) > mass_tolerance_
              ; ++fg_pointer )
        {
          if (fg_pointer == features.size())
          {
//            ++fg_pointer;
            break;
          }
          auto tm = fg_pointer;
        }
//        --fg_pointer;

        // if reached at the end of features
        if (fg_pointer == features.size())
        {
          features.push_back(*new_fg_iter);
          continue;
        }
      }

      // find closest & within RT range
      double smallest_diff(std::abs(features[fg_pointer].getMonoisotopicMass()-new_fg_iter->getMonoisotopicMass()));
      Size selected_index = fg_pointer;
      for (Size p = fg_pointer-1; p < fg_pointer+2; ++p)
      {
        if (p < 0 || p >= features.size())
        {
          continue;
        }

        // if out of wanted RT range, ignore
        if( !doFWHMbordersOverlap(features[p].getFwhmRange(), new_fg_iter->getFwhmRange()))
        {
          continue;
        }

        double diff = std::fabs(features[p].getMonoisotopicMass() - new_fg_iter->getMonoisotopicMass());
        if (diff < smallest_diff)
        {
          selected_index = p;
          smallest_diff = diff;
        }
      }
      fg_pointer = selected_index;

      // if smallest_diff is not smaller than tolerance, add the new_fg to features
      if (smallest_diff > mass_tolerance_ ||
          !doFWHMbordersOverlap(features[fg_pointer].getFwhmRange(), new_fg_iter->getFwhmRange()) )
      {
        if (features[fg_pointer].getMonoisotopicMass() > new_fg_iter->getMonoisotopicMass())
        {
          --fg_pointer;
        }
        features.insert(features.begin()+fg_pointer, *new_fg_iter);
        continue;
      }

      // within mass_tolerance & rt overlaps -> combine
      // TODO : maybe only when overlaps?
      if ( features[fg_pointer].getTraceIndices() == new_fg_iter->getTraceIndices()) // if exactly same feature exist
      {
        continue;
      }

      auto trace_indices = features[fg_pointer].getTraceIndices();
      for (auto& fg : *new_fg_iter)
      {
        if ( std::find( trace_indices.begin(), trace_indices.end(), fg.getTraceIndex() )
              != trace_indices.end())
        {
          features[fg_pointer].push_back(fg);
        }
      }
      features[fg_pointer].setFwhmRange();
      features[fg_pointer].setTraceIndices();
      features[fg_pointer].updateMassesAndIntensity();
    }


    // sort features (for later)
    std::sort(features.begin(), features.end());
  }

  bool FLASHDeconvQuant::doFWHMbordersOverlap(const std::pair<double, double>& border1,
                                                  const std::pair<double, double>& border2) const
  {
    if ( (border1.first > border2.second) || (border2.first > border1.second))
      return false;

    const double overlap_length = std::min(border1.second, border2.second)-std::max(border1.first, border2.first);
    if ( (overlap_length/(border1.second-border1.first) < 0.7) &&
         (overlap_length/(border2.second-border2.first) < 0.7) ) return false;

    return true;
  }

  bool FLASHDeconvQuant::doMassTraceIndicesOverlap(const FeatureGroup& fg1, const FeatureGroup& fg2) const
  {
    // get overlapping charge states
    int min_overlapping_charge = std::max(get<0>(fg1.getChargeRange()), get<0>(fg2.getChargeRange()));
    int max_overlapping_charge = std::min(get<1>(fg1.getChargeRange()), get<1>(fg2.getChargeRange()));

    if(min_overlapping_charge > max_overlapping_charge) // no overlapping charge
    {
      return false;
    }

    // collect possible overlapping mass_traces based on charges
    std::vector<Size> mt_indices_1;
    std::vector<Size> mt_indices_2;
    mt_indices_1.reserve(fg1.size());
    mt_indices_2.reserve(fg2.size());
    for (Size fg1_idx = 0; fg1_idx < fg1.size(); ++fg1_idx)
    {
      if (fg1[fg1_idx].getCharge() >= min_overlapping_charge && fg1[fg1_idx].getCharge() <= max_overlapping_charge)
      {
        mt_indices_1.push_back(fg1[fg1_idx].getTraceIndex());
      }
    }
    for (Size fg2_idx = 0; fg2_idx < fg2.size(); ++fg2_idx)
    {
      if (fg2[fg2_idx].getCharge() >= min_overlapping_charge && fg2[fg2_idx].getCharge() <= max_overlapping_charge)
      {
        mt_indices_2.push_back(fg2[fg2_idx].getTraceIndex());
      }
    }
    std::sort(mt_indices_1.begin(), mt_indices_1.end());
    std::sort(mt_indices_2.begin(), mt_indices_2.end());

    Size min_vec_size = std::min(mt_indices_1.size(), mt_indices_2.size());
    std::vector<Size> inters_vec;
    inters_vec.reserve(min_vec_size);
//    std::vector<Size>::iterator  inters_it;
    std::set_intersection(mt_indices_1.begin(), mt_indices_1.end(), mt_indices_2.begin(), mt_indices_2.end(), std::back_inserter(inters_vec));

//    double overlap_p = inters_vec.size() / min_vec_size;
    double overlap_percentage = static_cast<double>(inters_vec.size()) / static_cast<double>(min_vec_size);
    // TODO : change this to overlapping only major cs?
    if(overlap_percentage < .2)
    {
      return false;
    }
    return true;
  }

  void FLASHDeconvQuant::updateFeatureGroupInformation_(FeatureGroup &fg) const
  {
    // update private members in FeatureGroup based on the changed LogMassTraces
    // algorithm based on MassFeatureTrace::findFeatures

    int min_feature_abs_charge = INT_MAX; // min feature charge
    int max_feature_abs_charge = INT_MIN; // max feature charge

    auto per_charge_intensity = std::vector<double>(charge_range_ + 1, 0);
//    auto per_charge_max_intensity = std::vector<double>(charge_range_ + 1, 0);
//    auto per_charge_mz = std::vector<double>(charge_range_ + 1, 0);
    auto per_isotope_intensity = std::vector<double>(iso_model_.getMaxIsotopeIndex(), 0);

    double max_intensity = 0;

    /// getting general information before scoring
    for (auto lmt : fg)
    {
      min_feature_abs_charge = min_feature_abs_charge < lmt.getCharge() ? min_feature_abs_charge : lmt.getCharge();
      max_feature_abs_charge = max_feature_abs_charge > lmt.getCharge() ? max_feature_abs_charge : lmt.getCharge();

      if (lmt.getIntensity() > max_intensity)
      {
        max_intensity = lmt.getIntensity();
      }

      per_charge_intensity[lmt.getCharge() - charge_lower_bound_] += lmt.getIntensity();
      per_isotope_intensity[lmt.getIsotopeIndex()] += lmt.getIntensity();
//      if (per_charge_max_intensity[lmt.getCharge() - charge_lower_bound_] > lmt.getIntensity())
//      {
//        continue;
//      }
//      per_charge_max_intensity[lmt.getCharge() - charge_lower_bound_] = lmt.getIntensity();
//      per_charge_mz[lmt.getCharge() - charge_lower_bound_] = lmt.getCentroidMz();
    }

    fg.setChargeRange(min_feature_abs_charge, max_feature_abs_charge);

    /// scoring starts here
    int offset = 0;

    double mass = fg.getMonoisotopicMass();
    double isotope_score = FLASHDeconvQuant::getIsotopeCosineAndDetermineIsotopeIndex(mass, per_isotope_intensity, offset, iso_model_);
    fg.setIsotopeCosine(isotope_score);

    // TODO : if less than min_isotope_cosine_, don't merge??
//    if (isotope_score < min_isotope_cosine_)
//    {
//      return;
//    }

    // TODO : remove? not used
    double c_score = getChargeFitScore_(per_charge_intensity);
    fg.setChargeScore(c_score);

    /// setting per_isotope_score
    auto iso_dist = iso_model_.get(mass);
    int iso_size = (int) iso_dist.size();
    float total_noise = .0;
    float total_signal = .0;

    auto current_charge_range = fg.getChargeRange();
    // re initialize private vectors
    fg.initializePerChargeVectors();

    for (int abs_charge = std::get<0>(current_charge_range);
         abs_charge <= std::get<1>(current_charge_range);
         ++abs_charge) {
      int j = abs_charge - charge_lower_bound_;//current_min_charge_;
      if (per_charge_intensity[j] <= 0) {
        continue;
      }
      auto current_per_isotope_intensities = std::vector<double>(iso_model_.getMaxIsotopeIndex(), 0);

      int min_isotope_index = iso_model_.getMaxIsotopeIndex();
      int max_isotope_index = 0;

      double max_intensity = .0;

      for (auto &peak: fg) {
        if (peak.getCharge() != abs_charge) {
          continue;
        }

        if (peak.getIsotopeIndex() > iso_size) {
          continue;
        }

        current_per_isotope_intensities[peak.getIsotopeIndex()] += peak.getIntensity();
        //sumIntensity += p.intensity;
        min_isotope_index = min_isotope_index < peak.getIsotopeIndex() ? min_isotope_index : peak.getIsotopeIndex();
        max_isotope_index = max_isotope_index < peak.getIsotopeIndex() ? peak.getIsotopeIndex() : max_isotope_index;

        if (max_intensity < peak.getIntensity()) {
          max_intensity = peak.getIntensity();
        }
        // sp += p.intensity * p.intensity;
      }
      if (max_intensity <= 0) {
        continue;
      }

      for (int k = min_isotope_index; k <= max_isotope_index; ++k) {
        if (k > iso_size) {
          break;
        }
      }

      double cos_score = getCosine_(current_per_isotope_intensities,
                                    min_isotope_index,
                                    max_isotope_index,
                                    iso_dist,
                                    iso_size,
                                    0);

      fg.setChargeIsotopeCosine(abs_charge, cos_score);
      fg.setChargeIntensity(abs_charge, per_charge_intensity[j]);
    }

    fg.setAvgPPMError(getAvgPPMError_(fg));
  }

  void FLASHDeconvQuant::refineFeatureGroups_(std::vector<FeatureGroup>& in_features)
  {
    // change min, max charges based on built FeatureGroups (for later use in scoring)
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;

    // output features
    std::vector<FeatureGroup> out_feature;
    out_feature.reserve(in_features.size());

    // sort features by masses
    std::sort(in_features.begin(), in_features.end());

    for (auto & f : in_features)
    {
      min_abs_charge = min_abs_charge < std::get<0>(f.getChargeRange()) ? min_abs_charge : std::get<0>(f.getChargeRange());
      max_abs_charge = max_abs_charge > std::get<1>(f.getChargeRange()) ? max_abs_charge : std::get<1>(f.getChargeRange());
    }
    charge_lower_bound_ = min_abs_charge;
    charge_upper_bound_ = max_abs_charge;
    charge_range_ = charge_upper_bound_ - charge_lower_bound_ + 1;

    // insert FeatureGroup with highest score to out_features, merge if other FeatureGroup exist within mass_tol
    while(!in_features.empty())
    {
      // get a feature with the highest IsotopeCosineScore
      auto candidate_fg = std::max_element(in_features.begin(), in_features.end(), CmpFeatureGroupByScore());

      // get all features within mass_tol from candidate FeatureGroup
      std::vector<FeatureGroup>::iterator low_it, up_it;
      FeatureGroup lower_fg(candidate_fg->getMonoisotopicMass() - mass_tolerance_);
      FeatureGroup upper_fg(candidate_fg->getMonoisotopicMass() + mass_tolerance_);

      low_it = std::lower_bound(in_features.begin(), in_features.end(), lower_fg);
      up_it = std::upper_bound(in_features.begin(), in_features.end(), upper_fg);

      // no matching in features (found only candidate itself)
      if (up_it-low_it == 1)
      {
        // save it to out_features
        out_feature.push_back(*candidate_fg);

        // remove candidate from features
        in_features.erase(candidate_fg);
        continue;
      }

      // check if found features are overlapping with the candidate feature
      std::vector<int> v_indices_to_remove;
      v_indices_to_remove.reserve(up_it-low_it);
      std::vector<Size> mt_indices_to_add;
      mt_indices_to_add.reserve( (up_it-low_it) * candidate_fg->size() );
      std::vector<LogMassTrace *> mts_to_add;
      mts_to_add.reserve( (up_it-low_it) * candidate_fg->size() );
      std::vector<FeatureGroup *> mts_origin_fg;
      mts_origin_fg.reserve( (up_it-low_it) * candidate_fg->size() );

      for (; low_it != up_it; ++low_it)
      {
        // if low_it is candidate feature, ignore
        if (candidate_fg == low_it)
        {
          v_indices_to_remove.push_back( low_it - in_features.begin() );
          continue;
        }

        // check if fwhm overlaps
        if (!doFWHMbordersOverlap( low_it->getFwhmRange(), candidate_fg->getFwhmRange() ))
        {
          continue;
        }

        // check if masstrace overlaps
        if (!doMassTraceIndicesOverlap( *low_it, *candidate_fg ))
        {
          continue;
        }

        // merge found feature to candidate feature
        auto trace_indices = candidate_fg->getTraceIndices();
        for (auto& new_mt : *low_it)
        {
          if ( std::find( trace_indices.begin(), trace_indices.end(), new_mt.getTraceIndex() ) == trace_indices.end())
          {
            auto mt_already_added_idx = std::find(mt_indices_to_add.begin(), mt_indices_to_add.end(), new_mt.getTraceIndex());

            if (mt_already_added_idx != mt_indices_to_add.end()) // found already added logmasstrace
            {
              int found_index = mt_already_added_idx-mt_indices_to_add.begin();
              auto& mt_already_added = mts_to_add[found_index];

              // check if charge and isotope_index_ is same
              if (mt_already_added->getCharge() == new_mt.getCharge() && mt_already_added->getIsotopeIndex() == new_mt.getIsotopeIndex())
              {
                continue;
              }

              // check if cs & iso position is already taken in candidate_feature
              if (candidate_fg->doesThisIsotopeInChargeExist(new_mt.getCharge(), new_mt.getIsotopeIndex()))
              {
                continue;
              }

              // TODO: how should i address this? this doesn't make sense...
              auto &origin_fg = mts_origin_fg[found_index];
              if (mt_already_added->getCharge() == new_mt.getCharge() && mt_already_added->getIsotopeIndex() != new_mt.getIsotopeIndex())
              {
                if (origin_fg->getIsotopeCosine() < low_it->getIsotopeCosine())
                {
                  mts_to_add[found_index] = &new_mt;
                  mts_origin_fg[found_index] = &(in_features[ low_it - in_features.begin() ]);
                }
                continue;
              }

              auto m = mt_already_added;
            }

            // check if cs & iso position is already taken in candidate_feature
            if (candidate_fg->doesThisIsotopeInChargeExist(new_mt.getCharge(), new_mt.getIsotopeIndex()))
            {
              continue;
            }

            mt_indices_to_add.push_back(new_mt.getTraceIndex());
            mts_to_add.push_back(&new_mt);
            mts_origin_fg.push_back(&(in_features[low_it-in_features.begin()]));
          }
        }
        // add index of found feature to "to_be_removed_vector"
        v_indices_to_remove.push_back(low_it - in_features.begin());
      }
      // add extra masstraces to candidate_feature
      FeatureGroup final_candidate_fg = *candidate_fg; // copy of candidate_feature
      for(auto &new_mt : mts_to_add)
      {
        final_candidate_fg.push_back(*new_mt);
      }
      final_candidate_fg.setFwhmRange();
      final_candidate_fg.setTraceIndices();
      final_candidate_fg.updateMassesAndIntensity();

      // TODO : re-score based on current feature info
      updateFeatureGroupInformation_(final_candidate_fg);

      // save it to out_features
      out_feature.push_back(final_candidate_fg);

      // remove candidate from features
      std::sort(v_indices_to_remove.begin(), v_indices_to_remove.end());

      std::vector<FeatureGroup> tmp_out_fgs;
      tmp_out_fgs.reserve( in_features.size() - v_indices_to_remove.size() );
      for (Size i = 0; i<in_features.size(); ++i )
      {
        if (std::find(v_indices_to_remove.begin(), v_indices_to_remove.end(), i)==v_indices_to_remove.end())
        {
          tmp_out_fgs.push_back(in_features[i]);
        }
      }
      in_features.swap(tmp_out_fgs);

//      for(auto& r_idx : v_indices_to_remove)
//      {
//        auto test = r_idx;
//        features.erase( features.begin() + r_idx );
//      }
    }

    in_features.swap(out_feature);
  }

  void FLASHDeconvQuant::buildMassTraceGroups_(std::vector<LogMassTrace> &log_mtraces, std::vector<FeatureGroup>& features)
  {
    /// group mass traces to spectrum
    Size last_starting_point = 0;
    Size last_index = log_mtraces.size();

    std::vector<std::pair<double, LogMassTrace*>> mt_rt_starts;
    std::vector<std::pair<double, LogMassTrace*>> mt_rt_ends;
    mt_rt_starts.reserve(log_mtraces.size());
    mt_rt_ends.reserve(log_mtraces.size());
    int counter = 0;

    // collect rt information from mtraces to generate spectrum
    for (auto &lmt : log_mtraces)
    {
      mt_rt_starts.push_back(std::make_pair(lmt.getFwhmStart(), &lmt));
      mt_rt_ends.push_back(std::make_pair(lmt.getFwhmEnd(), &lmt));
    }

    // sorting mass traces in rt
    std::sort(mt_rt_starts.begin(), mt_rt_starts.end());
    std::sort(mt_rt_ends.begin(), mt_rt_ends.end());

    std::vector<std::pair<double, LogMassTrace*>>::const_iterator rt_s_iter = mt_rt_starts.begin();
    std::vector<std::pair<double, LogMassTrace*>>::const_iterator rt_e_iter = mt_rt_ends.begin();
    auto end_of_iter = mt_rt_starts.end();
    double current_rt = mt_rt_starts[0].first;
    std::vector<LogMassTrace*> local_traces;
    local_traces.reserve(log_mtraces.size());

    while(rt_s_iter != end_of_iter)
    {
      // initial rt binning is 1 sec (for generating spectrum)
      current_rt += 1;

      // add mass traces within rt range
      bool is_new_mt_added = false;
      for(;rt_s_iter != end_of_iter && rt_s_iter->first <= current_rt; ++rt_s_iter)
      {
        local_traces.push_back(rt_s_iter->second);
        is_new_mt_added = true;
      }

      // if nothing is added, increase current_rt
      if (!is_new_mt_added)
      {
        continue;
      }

      // remove mass traces out of rt range
      for(;rt_e_iter != mt_rt_ends.end() && rt_e_iter->first < current_rt ; ++rt_e_iter)
      {
        local_traces.erase(std::remove_if(local_traces.begin(), local_traces.end(),
                                          [&rt_e_iter](auto const& p){return rt_e_iter->second==p; }));
      }

      // sort local traces in mz
      sort(local_traces.begin(), local_traces.end(), CmpLogMassTraceByMZ());

      std::vector<FeatureGroup> local_fgroup;
      getFeatureFromSpectrum_(local_traces, local_fgroup);
      ++counter; // to track the number of generated spectra
      if (local_fgroup.size() == 0)
      {
        continue;
      }

      for (auto &pg : local_fgroup) {
        sort(pg.begin(), pg.end());
        pg.setFwhmRange();
        pg.setTraceIndices();
      }

      /// TODO: check if local features are duplicate - if not, add it to features
      // check if feature already exist ( 80% of overlap), then add to the output
//      addFeatureGroup_(features, local_fgroup);

      features.insert(features.end(), local_fgroup.begin(), local_fgroup.end());
//      std::sort(features.begin(), features.end());
    }

    OPENMS_LOG_INFO << "# generated spec from mass traces : " << counter << endl;
    OPENMS_LOG_INFO << "# generated feature groups from mass traces : " << features.size() << endl;
  }

  void FLASHDeconvQuant::clusterFeatureGroups_(std::vector<FeatureGroup>& fgroups, std::vector<std::vector<Size>>& shared_m_traces) const
  {
    // test writing
    String out_path = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "groups.tsv";
    ofstream out;
    out.open(out_path, ios::out);

    // header
//    out << "feature_label\tcs\tscore\tquant\tbounding_box_pos\tbounding_box_width\tbounding_box_height\t"
//           "iso_position\tmasstrace_centroid_rts\tmasstrace_centroid_mzs\tclustername\n";
      out << "mono_mass\tclustername\tshared\tcentroid_mzs\trts\tmzs\tintys\n";


    // generate array with indices of shared mass traces
    for (Size fg_index = 0; fg_index<fgroups.size(); ++fg_index)
    {
      for (auto& mt_i : fgroups[fg_index].getTraceIndices())
      {
        shared_m_traces[mt_i].push_back(fg_index);
      }
    }

    // *********************************************************** //
    // Step 1 constructing hypergraph from featurehypotheses
    //        node = mass traces
    //        hyperedge = hypotheses
    // *********************************************************** //
    Size num_nodes = shared_m_traces.size();
    std::vector<bool> bfs_visited;
    bfs_visited.resize(num_nodes, false);
    std::queue<Size> bfs_queue;
    Size search_pos = 0; // keeping track of mass trace index to look for seed

    std::vector<FeatureGroup> out_features;
    out_features.reserve(fgroups.size());
    std::vector<Size> out_feature_idx;
    out_feature_idx.reserve(fgroups.size());

    // BFS
    // TODO : progress logger
    this->startProgress(0, shared_m_traces.size(), "clustering features based on the shared mass traces");
    Size progress = 0;
    Size cluster_counter = 0;
    while (true)
    {
      this->setProgress(progress);
      ++progress;
      // finding a seed 'shared_mass_trace' to start with (for constructing a cluster)
      bool finished = true;
      for (Size i = search_pos; i < num_nodes; ++i)
      {
        if (!bfs_visited[i])
        {
          bfs_queue.push(i);
          bfs_visited[i] = true;
          finished = false;
          search_pos = i + 1;
          break;
        }
      }
      if (finished) // if no possible seed is left
        break;

      set<Size> hypo_indices_in_current_cluster;

      while (!bfs_queue.empty())
      {
        Size i = bfs_queue.front(); // index of seed
        bfs_queue.pop();

        // get feature indices sharing this seed
        for (vector<Size>::const_iterator it = shared_m_traces[i].begin();
             it != shared_m_traces[i].end();
             ++it)
        {
          hypo_indices_in_current_cluster.insert(*it);
          FeatureGroup &current_fg = fgroups[*it];

          for (const auto &mt_index : current_fg.getTraceIndices())
          {
            if (!bfs_visited[mt_index])
            {
              bfs_queue.push(mt_index);
              bfs_visited[mt_index] = true;
            }
          }
        }
      }
      // current cluster
      if (hypo_indices_in_current_cluster.empty()) continue; // no cluster out of current mt
      cluster_counter++;
      if (hypo_indices_in_current_cluster.size() == 1){
        //         no conflict, but cannot happen.
        out_features.push_back(fgroups[*(hypo_indices_in_current_cluster.begin())]);
        out_feature_idx.push_back(*(hypo_indices_in_current_cluster.begin()));
//        continue;
      }
      String cluster_name = "cluster" + to_string(progress);
      resolveConflictInCluster_(fgroups, shared_m_traces, hypo_indices_in_current_cluster, out_features, out_feature_idx,
                                out, cluster_name);
    }
    this->endProgress();

    out_features.shrink_to_fit();
    out_feature_idx.shrink_to_fit();
    out.close();
    OPENMS_LOG_INFO << "#cluster :" << cluster_counter << endl;

    // calculating quants with returned out_features
//    for( auto curr_it = deconv_masses.begin(); curr_it != deconv_masses.end();  )
//    {
//      bool does_update_needed = false;
//      std::set<int> new_cs_set;
//      auto& f_idx_vec = curr_it->second.feature_idx;
//      for (auto f_idx_itr = f_idx_vec.begin(); f_idx_itr != f_idx_vec.end(); )
//      {
//        auto iter = std::find(out_feature_idx.begin(), out_feature_idx.end(), *f_idx_itr);
//        if (iter != out_feature_idx.end())
//        {
//          FeatureHypothesis& fh = out_features[iter-out_feature_idx.begin()];
//          curr_it->second.quant_values += fh.computeQuant();
//          new_cs_set.insert(fh.getCharge());
//          ++f_idx_itr;
//        }
//        else
//        {
//          curr_it->second.removeFeatureHypothesis(hypotheses[*f_idx_itr].getFeatureMass(), hypotheses[*f_idx_itr].getScore());
//          f_idx_itr = f_idx_vec.erase(f_idx_itr);
//          does_update_needed = true;
//        }
//      }
//      if (does_update_needed)
//      {
//        curr_it->second.charges = new_cs_set;
//      }
//      filterDeconvMassStruct(deconv_masses, hypotheses, curr_it, does_update_needed);
//    }
//
//    OPENMS_LOG_INFO << "#final masses :" << deconv_masses.size() << endl;
  }

  void FLASHDeconvQuant::resolveConflictInCluster_(const std::vector<FeatureGroup>& in_features,
                                                       const std::vector<std::vector<Size> >& shared_m_traces_indices,
                                                       const std::set<Size>& hypo_indices,
                                                       std::vector<FeatureGroup>& out_features,
                                                       std::vector<Size>& out_feature_idx,
                                                       ofstream& out,
                                                       String& cluster_name) const
  {
    OPENMS_LOG_INFO << "-------- " << cluster_name <<  " -------- \n";
    for (const auto& feat_idx : hypo_indices)
    {
      OPENMS_LOG_INFO << feat_idx << "\t" << in_features[feat_idx].getMonoisotopicMass() << "\n";
    }

    for (const auto& f_idx : hypo_indices)
    {
//      double mz_upper_limit = 0.0;
//      double rt_upper_limit  = 0.0;
//      double mz_lower_limit = std::numeric_limits<double>::max();
//      double rt_lower_limit = std::numeric_limits<double>::max();

      auto feat = in_features[f_idx];

      double quant(0.0);

      auto mt_idxs = feat.getTraceIndices();

      for (const auto& mt : feat)
      {
//        const auto& boundingbox = mt.getMassTrace()->getConvexhull().getBoundingBox();
//
//        if (boundingbox.minX() < rt_lower_limit)
//          rt_lower_limit = boundingbox.minX();
//        if (boundingbox.maxX() > rt_upper_limit)
//          rt_upper_limit = boundingbox.maxX();
//        if (boundingbox.minY() < mz_lower_limit)
//          mz_lower_limit = boundingbox.minY();
//        if (boundingbox.maxY() > mz_upper_limit)
//          mz_upper_limit = boundingbox.maxY();

        // check if current mt is shared with other features
        bool isShared=false;
        if (shared_m_traces_indices[mt.getTraceIndex()].size() > 1)
        {
          isShared = true;
        }

        stringstream rts;
        stringstream mzs;
        stringstream intys;

        for (const auto& peak : *(mt.getMassTrace()))
        {
          mzs << peak.getMZ() << ",";
          rts << peak.getRT() << ",";
          intys << peak.getIntensity() << ",";
        }
        std::string peaks = rts.str();
        peaks.pop_back();
        peaks = peaks + "\t" + mzs.str();
        peaks.pop_back();
        peaks = peaks + "\t" + intys.str();
        peaks.pop_back();

        quant += mt.getMassTrace()->computePeakArea();

        out << std::to_string(feat.getMonoisotopicMass()) << "\t" << cluster_name << "\t"
            << isShared << "\t" << std::to_string(mt.getCentroidMz()) << "\t"
            << peaks + "\n";
      }

//      stringstream isos;
//      for (const auto& pair : feat.getIndicesOfMassTraces())
//      {
//        isos << pair.first << ",";
//      }
//      std::string iso_str = isos.str();
//      iso_str.pop_back();

//      String label = std::to_string(feat.get) + "\t" + std::to_string(feat.getCharge());
//      out << label << "\t" << to_string(feat.getScore()) << "\t" << std::to_string(quant) << "\t"
//          << to_string(rt_lower_limit) << "," << to_string(mz_lower_limit) << "\t"
//          << to_string((rt_upper_limit-rt_lower_limit)) << "\t" << to_string((mz_upper_limit-mz_lower_limit)) << "\t"
//          << iso_str << "\t" << centroids
//          << "\t" << cluster_name << "\n";
    }
  }

}