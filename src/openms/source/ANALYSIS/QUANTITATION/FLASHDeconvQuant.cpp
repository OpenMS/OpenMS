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
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>

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

  void FLASHDeconvQuant::writeFeatureGroupsInFile(std::vector<FeatureGroup>& fgroups, std::vector<MassTrace>& input_mtraces)
  {
    // to get "shared" information
    std::vector<std::vector<Size>> shared_m_traces(input_mtraces.size(), std::vector<Size>());
    for (Size fg_index = 0; fg_index<fgroups.size(); ++fg_index)
    {
      for (auto& mt_i : fgroups[fg_index].getTraceIndices())
      {
        shared_m_traces[mt_i].push_back(fg_index);
      }
    }

    std::vector<bool> fg_with_shared(fgroups.size(), false);
    for (auto &sharing_fgs : shared_m_traces)
    {
      if(sharing_fgs.size() < 2)
        continue;

      for (auto &fg_idx : sharing_fgs)
      {
        fg_with_shared[fg_idx] = true;
      }
    }

    ofstream out;
    OPENMS_LOG_INFO << "writing output..." << outfile_path << endl;
    out.open(outfile_path, ios::out);

    // header
    out << "mono_mass\tmin_charge\tmax_charge\tquant_value\tsummed_inty\tmost_abundant_charge\t"
           "quant_of_most_abundant_charge\tsummed_inty_of_most_abundant_cs\tarea_of_most_abundant_cs\t"
           "quant_of_apices_of_cs\tsummed_inty_of_apices_of_cs\t"
           "quant_of_top_three_cs\tsummed_inty_of_top_three_cs\tarea_of_top_three_cs\t"
           "all_area_under_the_curve\t"
           "charge_score\tiso_cosine\tweighted_iso_score\t"
           "fwhm_start\tfwhm_end\tfwhm_avg_length\trt_start\trt_end\tcentroid_rt_of_apices\trt_of_apex\t"
           "iso_length\tperIsotopeIntensity\tis_shared\n"; //

    int counter = -1;
    for (auto& fg : fgroups)
    {
      counter++;
      double quant_value = fg.getIntensity();

      // intys
      std::vector<double> per_iso_area = std::vector<double>(iso_model_.getMaxIsotopeIndex(), .0);
      std::vector<double> per_cs_intys = std::vector<double>(std::get<1>(fg.getChargeRange())+1 , .0);
      std::vector<double> per_cs_area = std::vector<double>(std::get<1>(fg.getChargeRange())+1 , .0);
      std::vector<double> per_cs_area_all = std::vector<double>(std::get<1>(fg.getChargeRange())+1 , .0);

      // rt range
      double min_rt = std::numeric_limits<double>::max();
      double max_rt = 0;

      // getting information while looping through mass traces in featuregroup
      for (auto& lmt : fg)
      {
        auto lmt_ptr = lmt.getMassTrace();
        double tmp_min = lmt_ptr->getConvexhull().getBoundingBox().minX();
        double tmp_max = lmt_ptr->getConvexhull().getBoundingBox().maxX();

        double int_before = (*lmt_ptr)[0].getIntensity();
        double rt_before = (*lmt_ptr)[0].getRT();

        double tmp_area = 0;
        for (auto& peaks : *lmt_ptr)
        {
          per_cs_intys[lmt.getCharge()]  += peaks.getIntensity();
          tmp_area += (int_before + peaks.getIntensity())/2 * (peaks.getRT() - rt_before);
          int_before = peaks.getIntensity();
          rt_before = peaks.getRT();
        }
        per_cs_area_all[lmt.getCharge()] += tmp_area;
        per_iso_area[lmt.getIsotopeIndex()] += lmt.getIntensity();
        per_cs_area[lmt.getCharge()] += lmt.getIntensity();

        if (tmp_min < min_rt)
          min_rt = tmp_min;
        if (tmp_max > max_rt)
          max_rt = tmp_max;
      }

      // get centroid rt of Apices
      fg.setCentroidRtOfApices();

      // get adundances of CS apices
      auto abundances_per_cs = fg.getSummedIntensityOfMostAbundantMTperCS();

      // get most abundant charge
      int most_abundant_cs = std::distance(per_cs_intys.begin(), std::max_element(per_cs_intys.begin(), per_cs_intys.end()));

      // intensities
      auto summedIntensities = std::accumulate(per_cs_intys.begin(), per_cs_intys.end(), .0);
      double summedInty_most_abundant_cs = *std::max_element(per_cs_intys.begin(), per_cs_intys.end());
      double area_most_abundant_cs = *std::max_element(per_cs_area.begin(), per_cs_area.end());
      double area_under_the_curve = std::accumulate(per_cs_area_all.begin(), per_cs_area_all.end(), .0);
      double all_area_most_abundant_cs = *std::max_element(per_cs_area_all.begin(), per_cs_area_all.end());

      // top3
      std::sort(per_cs_area.begin(), per_cs_area.end(), std::greater<double>());
      std::sort(per_cs_intys.begin(), per_cs_intys.end(), std::greater<double>());
      std::sort(per_cs_area_all.begin(), per_cs_area_all.end(), std::greater<double>());
      double top3_area = per_cs_area[0] + per_cs_area[1] + per_cs_area[2];
      double top3_intys = per_cs_intys[0] + per_cs_intys[1] + per_cs_intys[2];
      double top3_area_all = per_cs_area_all[0] + per_cs_area_all[1] + per_cs_area_all[2];

      auto iso_start = std::distance(per_iso_area.begin(), find_if( per_iso_area.begin(), per_iso_area.end(), [](double x) { return x != 0; }));
      auto iso_end = std::distance(per_iso_area.rbegin(), find_if( per_iso_area.rbegin(), per_iso_area.rend(), [](double x) { return x != 0; }));
      iso_end = per_iso_area.size() - iso_end;

      stringstream iso_ss;
      for(int i = iso_start; i < iso_end; ++i)
      {
        iso_ss << std::to_string(per_iso_area[i]) << ";";
      }
      std::string iso_ss_str = iso_ss.str();
      iso_ss_str.pop_back();

      out << std::to_string(fg.getMonoisotopicMass()) << "\t"
          << std::get<0>(fg.getChargeRange()) << "\t"
          << std::get<1>(fg.getChargeRange()) << "\t"
          << std::to_string(quant_value) << "\t"
          << std::to_string(summedIntensities) << "\t"
          << most_abundant_cs << "\t"
          << std::to_string(area_most_abundant_cs) << "\t"
          << std::to_string(summedInty_most_abundant_cs) << "\t"
          << std::to_string(all_area_most_abundant_cs) << "\t"
          << std::to_string(abundances_per_cs.first) << "\t"
          << std::to_string(abundances_per_cs.second) << "\t"
          << std::to_string(top3_area) << "\t"
          << std::to_string(top3_intys) << "\t"
          << std::to_string(top3_area_all) << "\t"
          << std::to_string(area_under_the_curve) << "\t"
          << std::to_string(fg.getChargeScore()) << "\t"
          << std::to_string(fg.getIsotopeCosine()) << "\t"
          << std::to_string(fg.getFeatureGroupScore()) << "\t"
          << std::to_string(fg.getFwhmRange().first) << "\t"
          << std::to_string(fg.getFwhmRange().second) << "\t"
          << std::to_string(fg.getAvgFwhmLength()) << "\t"
          << std::to_string(min_rt) << "\t"
          << std::to_string(max_rt) << "\t"
          << std::to_string(fg.getCentroidRtOfApices()) << "\t"
          << std::to_string(fg.getRtOfApex()) << "\t"
          << (iso_end - iso_start) << "\t"
          << iso_ss_str << "\t"
          << fg_with_shared[counter] << "\n";

//      if (counter == 42)
//      {
//        auto k = "1";
//      }
    }

    out.close();

    /// only for writing purpose belowe here :
    // writing mass traces (when no resolution is done)
    String out_path = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "groups.tsv";
    ofstream groups_out;
    groups_out.open(out_path, ios::out);
    groups_out << "mono_mass\tcharge\tisotope_index\tquant_value\tshared\tcentroid_mzs\trts\tmzs\tintys\n"; // header
    writeMassTracesOfFeatureGroup(fgroups, shared_m_traces, groups_out);
    groups_out.close();
  }

  void printSharedMassTraces(std::vector<std::vector<Size>>& shared_m_traces,
                             String outputPath,
                             std::vector<FeatureGroup>& fgroups,
                             std::vector<MassTrace>& input_mtraces)
  {
    ofstream out;
    out.open(outputPath, ios::out);
    out << "mass_trace_idx\tmass_trace_centroid_mz\tmass_trace_centroid_rt\tfeature_indices\tfeature_masses\tfeature_cs_iso\trts_or_apices\n" ;

    for (Size mt_idx = 0; mt_idx < shared_m_traces.size(); ++mt_idx)
    {
      if (shared_m_traces[mt_idx].size() < 2)
        continue;

      out << mt_idx << "\t"
          << std::to_string(input_mtraces[mt_idx].getCentroidMZ()) << "\t"
          << std::to_string(input_mtraces[mt_idx].getCentroidRT()) << "\t";

      stringstream f_indices;
      stringstream f_masses;
      stringstream f_infos;
      stringstream f_rts;
      for (auto & f_idx : shared_m_traces[mt_idx])
      {
        auto& feature = fgroups[f_idx];

        feature.setCentroidRtOfApices();
        f_indices << f_idx << ", ";
        f_masses << std::to_string(feature.getMonoisotopicMass()) << ", ";
        f_rts << std::to_string(feature.getCentroidRtOfApices()) << ", ";

        for (auto& mt : feature)
        {
          if (mt.getTraceIndex() == mt_idx)
          {
            f_infos << mt.getCharge() << "(" << mt.getIsotopeIndex() << ")" << ", ";
            break;
          }
        }
      }
      std::string f_str = f_indices.str();
      f_str.pop_back(); // comma
      f_str.pop_back(); // white space
      f_str = f_str + "\t" + f_masses.str();
      f_str.pop_back();
      f_str.pop_back();
      f_str = f_str + "\t" + f_infos.str();
      f_str.pop_back();
      f_str.pop_back();
      f_str = f_str + "\t" + f_rts.str();
      f_str.pop_back();
      f_str.pop_back();
      out << f_str << "\n";
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

    // filter out mass traces in each features and re-score them
//    std::vector<FeatureGroup> tmpFGroups;
//    tmpFGroups.swap(features);
//    features.reserve(tmpFGroups.size());
//    for (auto& f : tmpFGroups)
//    {
//      f.filterIsoMassTracesWithLowIntensities();
//      if(!rescoreFeatureGroup_(f))
//      {
//        continue;
//      }
//      f.setCentroidRtOfApices();
//      features.push_back(f);
//    }
    OPENMS_LOG_INFO << "total #feature groups : " << features.size() << endl;

    // *********************************************************** //
    // Step 3 clustering features
    // *********************************************************** //
    // TODO : move this to cluster method later


//    resolveSharedMassTraces_simple(features, shared_m_traces, input_mtraces);
//    resolveSharedMassTraces(features, shared_m_traces, input_mtraces);
    clusterFeatureGroups_(features, input_mtraces);
    writeFeatureGroupsInFile(features, input_mtraces);
  }

  void FLASHDeconvQuant::logTransformMassTraces_(std::vector<MassTrace> &input_mtraces, std::vector<LogMassTrace> &log_mtraces)
  {
    // shortest mass trace length? -> update local_rt_range_;
    double shortest_mt_length = local_rt_range_;
    double min_mz = std::numeric_limits<double>::max();
    double max_mz = 0;
    Size index = 0;
    double sum_intensity = .0;

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
      sum_intensity += tmp_lmt.getIntensity();
    }
    local_rt_range_ = shortest_mt_length;
    upper_bound_mz_ = max_mz;
    lower_bound_mz_ = min_mz;
    total_intensity_ = sum_intensity;

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
        harmonic_filter_matrix_.setValue(k, i, log(1.0 / (-1.0 * n / hc + (i + charge_lower_bound_))));
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

  double FLASHDeconvQuant::getBinValue_(const Size &bin, const double &min_value) const
  {
    return min_value + bin / mz_bin_width_;
  }

  Size FLASHDeconvQuant::getBinNumber_(const double &value, const double &min_value) const
  {
    if (value < min_value)
    {
      return 0;
    }
    return (Size) (((value - min_value) * mz_bin_width_) + .5);
  }

  // From log mz to mz bins. To reduce edge effect, its adjacent bin is set when a bin is set.
  void FLASHDeconvQuant::updateMzBins_(std::vector<LogMassTrace*> &local_traces, const Size &bin_number,
                                       const double& mz_bin_min, std::vector<float> &mz_bin_intensities)
  {
    mz_bins_for_edge_effect_ = boost::dynamic_bitset<>(bin_number);
    mz_bins_ = boost::dynamic_bitset<>(bin_number);
    for (auto &trace : local_traces)
    {
      Size bi = getBinNumber_(trace->getLogCentroidMz(), mz_bin_min);
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
      Size bi = getBinNumber_(trace->getLogCentroidMz(), mz_bin_min);
      double delta = (trace->getLogCentroidMz() - getBinValue_(bi, mz_bin_min));

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
        if ((Size) j >= mass_bins_.size()-1)
        {
          break;
        }
        mass_bins_[j-1] = true;
        mass_bins_[j] = true;
        mass_bins_[j+1] = true;
      }
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
    Size mass_bin_size = mass_bins_.size();
    int log_mz_peak_size = (int) log_mtraces.size();
    // this stores which peak is now being considered per charge. Per charge, peak is considered from left (lowest m/z) to right (highest m/z).
    auto current_peak_index = std::vector<int>(charge_range_, 0);
    fgroup.reserve(mass_bins_.count());
    Size mass_bin_index = mass_bins_.find_first();
    auto peak_bin_numbers = std::vector<Size>(log_mz_peak_size);

    // per peak, store the m/z bin number for fast processing
    for (int i = 0; i < log_mz_peak_size; ++i)
    {
      peak_bin_numbers[i] = getBinNumber_(log_mtraces[i]->getLogCentroidMz(), mz_bin_min_value_);
    }

    // main iteration. per_mass_abs_charge_ranges gives the range of charges for each mass bin
    while (mass_bin_index != mass_bins_.npos)
    {
      double log_m = getBinValue_(mass_bin_index, mass_bin_min_value_);
      double mass = exp(log_m);
      FeatureGroup fg(charge_lower_bound_,
                   per_mass_abs_charge_ranges.getValue(1, mass_bin_index) + charge_lower_bound_);
      auto test = per_mass_abs_charge_ranges.getValue(1, mass_bin_index);

      fg.reserve(charge_range_ * max_nr_traces_);
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

        while (cpi < log_mz_peak_size - 1) // scan through peaks from cpi
        {
          if (peak_bin_numbers[cpi] ==
              b_index) // if the peak of consideration matches to this mass with charge abs_charge
          {
            double intensity = log_mtraces[cpi]->getIntensity();
            if (intensity >
                max_intensity) // compare with other matching peaks and select the most intense peak (in max_peak_index)
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

        // now we have a mathcing peak for this mass of charge  abs_charge. From here, istope peaks are collected
        const double mz = log_mtraces[max_peak_index]->getCentroidMz();
        const double iso_delta = Constants::ISOTOPE_MASSDIFF_55K_U / (abs_charge);
        double mz_delta = mz_tolerance_ * mz; //

        // double peak_pwr = .0;

        double max_mz = mz;
        double max_peak_intensity = log_mtraces[max_peak_index]->getIntensity();
        for (int peak_index = max_peak_index; peak_index < log_mz_peak_size; peak_index++)
        {
          const double observed_mz = log_mtraces[peak_index]->getCentroidMz();

          const double intensity = log_mtraces[peak_index]->getIntensity();
          double mz_diff = observed_mz - mz;

          int tmp_i = (int) (.5 + mz_diff / iso_delta);

          if (max_peak_intensity < intensity)
          {
            max_peak_intensity = intensity;
            max_mz = observed_mz;
          }

          int tmp_i_for_stop = (int) (.5 + (observed_mz - max_mz) / iso_delta);;
          if (tmp_i_for_stop > (int) right_index)
          {
            break;
          }

          // peak_pwr += intensity * intensity;
          if (abs(mz_diff - tmp_i * iso_delta) < mz_delta) // noise   max_intensity  vs   intensity
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size)
            {
              LogMassTrace p(*log_mtraces[peak_index]);
              p.setCharge(abs_charge);
              p.setIsotopeIndex(tmp_i);
              fg.push_back(p);
            }
          }
        }

        for (int peak_index = max_peak_index - 1; peak_index >= 0; peak_index--)
        {
          const double observed_mz = log_mtraces[peak_index]->getCentroidMz();
          const double intensity = log_mtraces[peak_index]->getIntensity();

          //observedMz = mz + isof * i * d - d * mzDelta;
          double mz_diff = mz - observed_mz;
          int tmp_i = (int) (.5 + mz_diff / iso_delta);

          if (max_peak_intensity < intensity)
          {
            max_peak_intensity = intensity;
            max_mz = observed_mz;
          }
          int tmp_i_for_stop = (int) (.5 + (max_mz - observed_mz) / iso_delta);;
          if (tmp_i_for_stop > (int) left_index)
          {
            break;
          }

          // peak_pwr += intensity * intensity;
          if (abs(mz_diff - tmp_i * iso_delta) >= mz_delta)
          {
            const Size bin = peak_bin_numbers[peak_index] + bin_offset;
            if (bin < mass_bin_size)
            {
              LogMassTrace p(*log_mtraces[peak_index]);
              p.setCharge(abs_charge);
              p.setIsotopeIndex(tmp_i);
              fg.push_back(p);
            }
          }
        }
        // TODO : revive charge snr?
        // pg.setChargePower(abs_charge, peak_pwr);
      }

      if (!fg.empty())
      {

        double max_intensity = -1.0;
        //double sum_intensity = .0;
        double t_mass = .0;
        auto new_peaks = std::vector<LogMassTrace>();
        new_peaks.reserve(fg.size());

        for (auto &p: fg)
        {
          if (max_intensity < p.getIntensity())
          {
            max_intensity = p.getIntensity();
            t_mass = p.getUnchargedMass();
          }
        }

        double iso_tolerance = mz_tolerance_ * t_mass;
        int min_off = 10000;
        int max_off = -1;
        int max_charge = -1;
        //        std::vector<double> noise_power(charge_upper_bound_ + 1, .0);
        for (auto &p: fg)
        {
          p.setIsotopeIndex(round((p.getUnchargedMass() - t_mass) / Constants::ISOTOPE_MASSDIFF_55K_U));
          if (abs(t_mass - p.getUnchargedMass() + Constants::ISOTOPE_MASSDIFF_55K_U * p.getIsotopeIndex()) >
              iso_tolerance)
          {
            continue;
          }
          new_peaks.push_back(p);
          min_off = min_off > p.getIsotopeIndex() ? p.getIsotopeIndex() : min_off;
          max_off = max_off < p.getIsotopeIndex() ? p.getIsotopeIndex() : max_off;
          max_charge = max_charge < p.getCharge() ? p.getCharge() : max_charge;
          //          signal_power[p.abs_charge] += p.intensity * p.intensity;
        }
        if (max_charge < low_charge_ && min_off == max_off)
        {
          continue;
        }

        fg.swap(new_peaks);

        for (auto &p: fg)
        {
          p.setIsotopeIndex(p.getIsotopeIndex() - min_off);
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
    int h_charge_size = (int) harmonic_charges_.size();
    int min_peak_cntr = min_nr_mtraces_;
    long bin_end = (long) mass_bins_.size();
    // how many peaks of continuous charges per mass
    auto support_peak_count = std::vector<int>(mass_bins_.size(), 0);

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
    auto prev_charges = std::vector<int>(mass_bins_.size(), charge_range_ + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prev_intensities = std::vector<float>(mass_bins_.size(), 1.0f);

    // intensity change ratio should not exceed the factor.
    const float factor = 5.0;
    // intensity ratio between consecutive charges for possible harmonic should be within this factor
    const float hfactor = 2.0;
    for (int i = mz_bin_index_reverse.size() - 1; i >= 0; i--)
    {
      mz_bin_index = mz_bin_index_reverse[i];
      float intensity = mz_intensities[mz_bin_index];
      double mz = -1.0, log_mz = 0;
      log_mz = getBinValue_(mz_bin_index, mz_bin_min_value_);
      mz = exp(log_mz);

      // scan through charges
      for (int j = 0; j < charge_range_; j++) //  increasing.
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
        bool charge_not_continous = prev_charge - j != -1 && (prev_charge <= charge_range_);
        bool pass_first_check = false;

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
        if (!pass_first_check && abs_charge <= low_charge_)
        { // for low charges
          double diff = Constants::ISOTOPE_MASSDIFF_55K_U / abs_charge / mz;
          Size next_iso_bin = getBinNumber_(log_mz + diff, mz_bin_min_value_);

          if (next_iso_bin < mz_bins_.size() && mz_bins_[next_iso_bin])
          {
            pass_first_check = true;
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
              max_intensity = intensity; //
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
              for (int off = 0; off <= 0; off++)
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
                    break;
                  }
                }
              }
              if (is_harmonic)
              {
                break;
              }
            }
            if (!is_harmonic) // if it is not harmonic
            {
              mass_intensitites[mass_bin_index] += intensity + support_peak_intensity;
              spc++;
              if (spc >= min_peak_cntr || spc >= abs_charge / 2) //
              {
                mass_bins_[mass_bin_index] = true;
              }
            }
            else // if harmonic
            {
              mass_intensitites[mass_bin_index] -= intensity;
            }
          }
          else // for low charge, include the mass if isotope is present
          {
            mass_intensitites[mass_bin_index] += intensity;
            spc++;
            //if (spc >= min_peak_cntr || spc >= abs_charge * .75) //
            {
              mass_bins_[mass_bin_index] = true;
            }
          }
        }
        prev_intensity = intensity;
        prev_charge = j;
      }
    }
  }

  // Subfunction of updateMassBins_. If a peak corresponds to multiple masses, only one mass is selected for the peak based on intensities..
  // mass level harmonic check is also performed in this function
  // it also outputs the charge range of each mass bin
  Matrix<int> FLASHDeconvQuant::filterMassBins_(const std::vector<float> &mass_intensities)
  {
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
      int max_intensity_abs_charge = 0;

      for (int j = 0; j < charge_range_; j++)
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

          int abs_charge = (j + charge_lower_bound_);

          double original_log_mass = getBinValue_(mass_bin_index, mass_bin_min_value_);
          double mass = exp(original_log_mass);
          for (int iso_off = -2; iso_off <= 2 && !artifact; ++iso_off)
          {
            double log_mass = log(mass + iso_off * Constants::ISOTOPE_MASSDIFF_55K_U);   //original_log_mass + diff * iso_off;
            if (log_mass < 1)
            {
              continue;
            }
            for (int h = 2; h <= 3 && !artifact; h++)
            {
              for (int f = -1; f <= 1 && !artifact; f += 2) //
              {
                double hmass = log_mass - log(h) * f;
                Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_);
                if (hmass_index > 0 && hmass_index < mass_bins_.size() - 1)
                {
                  if (mass_intensities[hmass_index] >= t)
                  {
                    artifact = true; //
                    break;
                  }
                }
              }
            }
            //   max_intensity_abs_charge off by one here
            if (!artifact)
            {
              for (int coff = 1; coff <= 2 && !artifact; coff++)
              {
                for (int f = -1; f <= 1 && !artifact; f += 2)
                {
                  if (abs_charge + f * coff <= 0)
                  {
                    continue;
                  }
                  double hmass = log_mass - log(abs_charge) + log(abs_charge + f * coff);
                  Size hmass_index = getBinNumber_(hmass, mass_bin_min_value_);
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

  bool FLASHDeconvQuant::checkChargeDistribution_(const std::vector<double> &per_charge_intensity) const {

//    double max_per_charge_intensity = *std::max_element(per_charge_intensity.begin(), per_charge_intensity.end());
//    double sum_charge_intensity = .0;
    int cntr = 0;
    int non_zero_start = -1, non_zero_end = 0;
    for (int i = 0; i < charge_range_; i++) {
      if (per_charge_intensity[i] > 0) {
//        sum_charge_intensity += per_charge_intensity[i];
        cntr++;
//        max_per_charge_intensity = std::max(max_per_charge_intensity, per_charge_intensity[i]);
        if (non_zero_start < 0) {
          non_zero_start = i;
        }
        non_zero_end = i;
      }
    }
    if (cntr < min_nr_mtraces_) { // less than 3 charges in this FeatureGroup
      return false;
    }

    int prev_charge = non_zero_start;
    int n_r = 0;
//    float factor = 1;
//    double int_threshold = (sum_charge_intensity / cntr) * factor;
    for (int k = prev_charge + 1; k <= non_zero_end; k++) {
      if (per_charge_intensity[k] <= 0) {
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

  void FLASHDeconvQuant::calculatePerChargeIsotopeIntensity_(std::vector<double> &per_isotope_intensity,
                                                              std::vector<double> &per_charge_intensity,
                                                              const int max_isotope_count,
                                                              FeatureGroup &fg) const {
    int min_pg_charge = INT_MAX;
    int max_pg_charge = INT_MIN;

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

//    return std::vector<int>{max_intensity_charge_index, max_intensity_iso_index};
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
//      double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i + 1]);
      if (diff <= 0) {
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
    for (int j = a_start; j <= a_end; j++)
    {
      int i = j - offset;
      a_norm += a[j] * a[j];

      if (i < 0 || i >= b_size || b[i].getIntensity() <= 0)
      {
        continue;
      }
      n += a[j] * b[i].getIntensity(); //
    }
    if (a_norm <= 0)
    {
      return 0;
    }
    return n / sqrt(a_norm);
  }

  double FLASHDeconvQuant::getIsotopeCosineAndDetermineIsotopeIndex(const double mono_mass,
                                                                    const std::vector<double> &per_isotope_intensities,
                                                                    int &offset) const
  {
    auto iso = iso_model_.get(mono_mass);

    int iso_size = (int) iso.size();
    int apex_index = iso_model_.getApexIndex(mono_mass);
    int iso_range = iso_model_.getRightCountFromApex(mono_mass) + iso_model_.getLeftCountFromApex(mono_mass);

    offset = 0;
    double min_diff = -1000;
    int max_isotope_index = 0, min_isotope_index = -1;

    for (int i = 0; i < iso_model_.getMaxIsotopeIndex(); i++)
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

    for (int tmp_offset = -apex_index - 1; tmp_offset <= -apex_index + max_isotope_index + 1; tmp_offset++)
    {
      double tmp_cos = getCosine_(per_isotope_intensities,
                                  min_isotope_index,
                                  max_isotope_index,
                                  iso,
                                  iso_size, //apex_index,
                                  tmp_offset);

      if (min_diff < tmp_cos) {
        min_diff = tmp_cos;
        offset = tmp_offset;
      }
    }

    return getCosine_(per_isotope_intensities,
                      min_isotope_index,
                      max_isotope_index,
                      iso,
                      iso_size,
                      offset);
  }

  float FLASHDeconvQuant::getAvgPPMError_(FeatureGroup &pg) const{
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
      if (scoreFeatureGroup_(feature_group))
      {
        filtered_peak_groups.push_back(feature_group);
      }
    }
    local_fgroup.swap(filtered_peak_groups);

    // TODO : revive here?
//    filterPeakGroupsByIsotopeCosine_(max_c);
  }

  void FLASHDeconvQuant::removeOverlappingPeakGroups_(std::vector<FeatureGroup> &local_fgroup) const{
    int iso_length = 1;

    std::vector<FeatureGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(local_fgroup.size());
    sort(local_fgroup.begin(), local_fgroup.end());

    filtered_pg_vec.push_back(local_fgroup[0]); // when i = 0
    for (Size i = 1; i < local_fgroup.size(); ++i) {
      if (abs(local_fgroup[i - 1].getMonoisotopicMass() - local_fgroup[i].getMonoisotopicMass()) < 1e-3
          &&
          local_fgroup[i - 1].getIntensity() >= local_fgroup[i].getIntensity()) {
        continue;
      }

      filtered_pg_vec.push_back(local_fgroup[i]);
    }
    local_fgroup.swap(filtered_pg_vec);
    std::vector<FeatureGroup>().swap(filtered_pg_vec);
    filtered_pg_vec.reserve(local_fgroup.size());

    for (Size i = 0; i < local_fgroup.size(); i++) {
      bool select = true;
      auto &pg = (local_fgroup)[i];

      if (pg.getMonoisotopicMass() <= 0) {
        continue;
      }
      double mass_tolerance = pg.getMonoisotopicMass() * mz_tolerance_ * 2;

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
        for (; j >= 0; j--)
        {
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
    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;

    double mz_bin_max_value = local_traces[local_traces.size()-1]->getLogCentroidMz();
    double mass_bin_max_value = std::min( mz_bin_max_value - filter_[tmp_peak_cntr],
                                          log(max_mass_ + iso_model_.getRightCountFromApex(max_mass_) + 1));

    tmp_peak_cntr = min_nr_mtraces_ - 1;
    tmp_peak_cntr = tmp_peak_cntr < 0 ? 0 : tmp_peak_cntr;
    mass_bin_min_value_ = log(std::max(1.0, min_mass_ - iso_model_.getAverageMassDelta(min_mass_)));
    mz_bin_min_value_ = local_traces[0]->getLogCentroidMz();
    Size mass_bin_number = getBinNumber_(mass_bin_max_value, mz_bin_min_value_) + 1;

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

    Size mz_bin_number = getBinNumber_(mz_bin_max_value, mz_bin_min_value_) + 1;
    auto mz_bin_intensities = std::vector<float>(mz_bin_number, .0f);

    // From log mz to mz bins.
    updateMzBins_(local_traces, mz_bin_number, mz_bin_min_value_, mz_bin_intensities);

    // take the mass bins from previous overlapping spectra and put them in the candidate mass bins. // TODO: remove this? or??
    mass_bins_ = boost::dynamic_bitset<>(mass_bin_number);
    unionPrevMassBins_();
    auto per_mass_abs_charge_ranges = updateMassBins_(mz_bin_intensities);

    getCandidatePeakGroups_(local_traces, per_mass_abs_charge_ranges, local_fgroup);

    if (local_fgroup.size() > 0)
    {
      scoreAndFilterPeakGroups_(local_fgroup);
    }
    if (local_fgroup.size() > 1)
    {
      removeOverlappingPeakGroups_(local_fgroup);
    }

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

  bool FLASHDeconvQuant::doFWHMbordersOverlap(const std::pair<double, double>& border1,
                                                  const std::pair<double, double>& border2) const
  {
    if ( (border1.first > border2.second) || (border2.first > border1.second))
      return false;

    const double overlap_length = std::min(border1.second, border2.second)-std::max(border1.first, border2.first);
    if ( (overlap_length/(border1.second-border1.first) < 0.5) &&
         (overlap_length/(border2.second-border2.first) < 0.5) ) return false;

    return true;
  }

  bool FLASHDeconvQuant::doMassTraceIndicesOverlap(const FeatureGroup& fg1, const FeatureGroup& fg2) const
  {
    // get overlapping charge states
    int min_overlapping_charge = std::max(get<0>(fg1.getChargeRange()), get<0>(fg2.getChargeRange()));
    int max_overlapping_charge = std::min(get<1>(fg1.getChargeRange()), get<1>(fg2.getChargeRange()));

    if(min_overlapping_charge > max_overlapping_charge) // no overlapping charge
    {
      return true;
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
//    double overlap_percentage = static_cast<double>(inters_vec.size()) / static_cast<double>(min_vec_size);
    // TODO : change this to overlapping only major cs?
    if(inters_vec.size() < 1)
    {
      return false;
    }
    return true;
  }

  bool FLASHDeconvQuant::rescoreFeatureGroup_(FeatureGroup &fg) const
  {
    // update private members in FeatureGroup based on the changed LogMassTraces

    if( !scoreFeatureGroup_(fg) )
    {
      return false;
    }

    fg.setFwhmRange();
    fg.setTraceIndices();
    return true;
  }

  bool FLASHDeconvQuant::scoreFeatureGroup_(FeatureGroup &fg, double given_iso_score) const
  {
    // return false when scoring is not done (filtered out)

    if (given_iso_score == 0)
    {
      given_iso_score = min_isotope_cosine_;
    }

    auto per_isotope_intensities = std::vector<double>(iso_model_.getMaxIsotopeIndex(), 0);
    auto per_charge_intensities = std::vector<double>(charge_range_, 0);

    calculatePerChargeIsotopeIntensity_(per_isotope_intensities, per_charge_intensities,
                                        iso_model_.getMaxIsotopeIndex(), fg);

    // if the number of charges are not enough
    if (std::get<1>(fg.getChargeRange()) - std::get<0>(fg.getChargeRange()) < min_nr_mtraces_ )
    {
      return false;
    }

    /// isotope cosine calculation
    int offset = 0;
    double isotope_score = getIsotopeCosineAndDetermineIsotopeIndex(fg.getMonoisotopicMass(),
                                                                    per_isotope_intensities,
                                                                    offset);
    fg.setIsotopeCosine(isotope_score);
    if (fg.empty() || isotope_score < given_iso_score) // min_isotope_cosine_
    {
      return false;
    }

    // TODO : remove? does charge distribution makes big difference?
    double cs = getChargeFitScore_(per_charge_intensities);
    fg.setChargeScore(cs);

    // NOTE : if not using cs dist, too many harmonics
    bool is_charge_well_distributed = checkChargeDistribution_(per_charge_intensities);
//    //double tmp = getChargeFitScore_(per_abs_charge_intensities, charge_range);
    if (!is_charge_well_distributed) {
        return false;
    }
    /// TODO_ends_here

    /// update monoisotopic mass of this FeatureGroup
    fg.updateMassesAndIntensity(offset, iso_model_.getMaxIsotopeIndex());
    if (fg.getMonoisotopicMass() < min_mass_ || fg.getMonoisotopicMass() > max_mass_) {
      return false;
    }

    /// setting per_isotope_score
    auto iso_dist = iso_model_.get(fg.getMonoisotopicMass());
    int iso_size = (int) iso_dist.size();

    // initialize private vectors
    fg.initializePerChargeVectors();

    double weighted_iso_score = fg.getIntensity() / total_intensity_;

    auto current_charge_range = fg.getChargeRange();
    for (int abs_charge = std::get<0>(current_charge_range);
         abs_charge <= std::get<1>(current_charge_range);
         ++abs_charge) {
      int j = abs_charge - charge_lower_bound_;//current_min_charge_;
      if (per_charge_intensities[j] <= 0) {
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

        if (max_intensity < peak.getIntensity())
        {
          max_intensity = peak.getIntensity();
        }
      }
      if (max_intensity <= 0) {
        continue;
      }

      double cos_score = getCosine_(current_per_isotope_intensities,
                                    min_isotope_index,
                                    max_isotope_index,
                                    iso_dist,
                                    iso_size,
                                    0);

      weighted_iso_score +=  (per_charge_intensities[j] / total_intensity_) * cos_score;

      fg.setChargeIsotopeCosine(abs_charge, cos_score);
      fg.setChargeIntensity(abs_charge, per_charge_intensities[j]);
    }

    fg.setFeatureGroupScore(weighted_iso_score);
    fg.setAvgPPMError(getAvgPPMError_(fg));
    return true;
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

      // calculate most intensive mass traces (per charge, per iso) in each feature
    }
    charge_lower_bound_ = min_abs_charge;
    charge_upper_bound_ = max_abs_charge;
    charge_range_ = charge_upper_bound_ - charge_lower_bound_ + 1;

    Size initial_size = in_features.size();
    this->startProgress(0, initial_size, "refining feature groups");
    // insert FeatureGroup with highest score to out_features, merge if other FeatureGroup exist within mass_tol
    while(!in_features.empty())
    {
      this->setProgress(initial_size-in_features.size());

      // get a feature with the highest IsotopeCosineScore
      auto candidate_fg = std::max_element(in_features.begin(), in_features.end(), CmpFeatureGroupByScore());

      // get all features within mass_tol from candidate FeatureGroup
      std::vector<FeatureGroup>::iterator low_it, up_it;
//      double mass_tolerance = .0;
//      if (candidate_fg->getMonoisotopicMass() > 10000)
//      {
//        mass_tolerance = mass_tolerance_da_;
//      }
//      else{
//        mass_tolerance = candidate_fg->getMonoisotopicMass() * mass_tolerance_ppm_ * 1e-6;
//      }

//      if (candidate_fg->getMonoisotopicMass() < 13475  && candidate_fg->getMonoisotopicMass() > 13470)
//        auto& here = candidate_fg;

      FeatureGroup lower_fg(candidate_fg->getMonoisotopicMass() - mass_tolerance_da_);
      FeatureGroup upper_fg(candidate_fg->getMonoisotopicMass() + mass_tolerance_da_);

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
      std::set<Size> mt_indices_to_add;
//      mt_indices_to_add.reserve( (up_it-low_it) * candidate_fg->size() );
      std::vector<LogMassTrace *> mts_to_add;
      mts_to_add.reserve( (up_it-low_it) * candidate_fg->size() );
//      std::vector<FeatureGroup *> mts_origin_fg;
//      mts_origin_fg.reserve( (up_it-low_it) * candidate_fg->size() );

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
          // if this mass trace is not used in candidate_fg
          if ( std::find( trace_indices.begin(), trace_indices.end(), new_mt.getTraceIndex() ) == trace_indices.end())
          {
//            auto mt_already_added_idx = std::find(mt_indices_to_add.begin(), mt_indices_to_add.end(), new_mt.getTraceIndex());
//
//            // if this mass trace is already added to candidate_fg by other merged fg
//            if (mt_already_added_idx != mt_indices_to_add.end())
//            {
//              int found_index = mt_already_added_idx-mt_indices_to_add.begin();
//              auto& mt_already_added = mts_to_add[found_index];
//
//              // check if charge and isotope_index are same
//              if (mt_already_added->getCharge() == new_mt.getCharge() && mt_already_added->getIsotopeIndex() == new_mt.getIsotopeIndex())
//              {
//                continue;
//              }
//
//              // check if cs & iso position is already taken in candidate_feature
//              if (candidate_fg->doesThisIsotopeInChargeExist(new_mt.getCharge(), new_mt.getIsotopeIndex())!= nullptr)
//              {
//                continue;
//              }
//
//              // TODO: how should i address this? this doesn't make sense...
//              auto &origin_fg = mts_origin_fg[found_index];
//              if (mt_already_added->getCharge() == new_mt.getCharge() && mt_already_added->getIsotopeIndex() != new_mt.getIsotopeIndex())
//              {
//                if (origin_fg->getIsotopeCosine() < low_it->getIsotopeCosine())
//                {
//                  mts_to_add[found_index] = &new_mt;
//                  mts_origin_fg[found_index] = &(in_features[ low_it - in_features.begin() ]);
//                }
//                continue;
//              }
//            }
//
//            // check if cs & iso position is already taken in candidate_feature
//            const LogMassTrace* existing_mt = candidate_fg->doesThisIsotopeInChargeExist(new_mt.getCharge(), new_mt.getIsotopeIndex());
//            if (existing_mt != nullptr)
//            {
//              // checking if this mt could replace the cs & iso position
//
//              // does new mt have higher intensity?
//              if (existing_mt->getIntensity() > new_mt.getIntensity())
//                continue;
//
//              // does fwhm overlap with candidate_fg?
//              auto border = candidate_fg->getFwhmRange();
//              if ( (border.first > new_mt.getFwhmEnd()) || (border.second < new_mt.getFwhmStart()))
//                continue;
//
//              // replace cs & iso position with current one (remove existing_mt)
//              for (auto lmt_iter = candidate_fg->begin(); lmt_iter != candidate_fg->end(); ++lmt_iter)
//              {
//                if (existing_mt->getTraceIndex() != lmt_iter->getTraceIndex()) continue;
//                candidate_fg->erase(lmt_iter);
//                break;
//              }
//            }

            mt_indices_to_add.insert(new_mt.getTraceIndex());
            mts_to_add.push_back(&new_mt);
//            mts_origin_fg.push_back(&(in_features[low_it-in_features.begin()]));
          }
        }
        // add index of found feature to "to_be_removed_vector"
        v_indices_to_remove.push_back(low_it - in_features.begin());
      }

      // sort mts_to_add by abundance
      std::sort(mts_to_add.begin(), mts_to_add.end(), CmpLogMassTraceByIntensity());

      // add extra masstraces to candidate_feature
      FeatureGroup final_candidate_fg = *candidate_fg; // copy of candidate_feature
      for(auto &new_mt : mts_to_add)
      {
        // check if this mt of this lmt was already added to candidate_fg
        if (mt_indices_to_add.find(new_mt->getTraceIndex()) == mt_indices_to_add.end())
        {
          continue;
        }
        mt_indices_to_add.erase(new_mt->getTraceIndex());

        LogMassTrace* apex_lmt_in_this_cs = final_candidate_fg.getApexLMTofCharge(new_mt->getCharge());
        // if this mt is introducing new charge
        if (apex_lmt_in_this_cs == nullptr)
        {
          final_candidate_fg.push_back(*new_mt);
          continue;
        }

        /// re-calculate isotope index
        double apex_lmt_mass = apex_lmt_in_this_cs->getUnchargedMass();
        double new_lmt_mass = new_mt->getUnchargedMass();

        int tmp_iso_idx = round((new_lmt_mass - apex_lmt_mass) / Constants::ISOTOPE_MASSDIFF_55K_U);
        if (abs(apex_lmt_mass - new_lmt_mass + Constants::ISOTOPE_MASSDIFF_55K_U * tmp_iso_idx) >
            apex_lmt_mass * mz_tolerance_)
        {
          continue;
        }
        if (apex_lmt_in_this_cs->getIsotopeIndex() + tmp_iso_idx < 0)
        {
          continue;
        }
        new_mt->setIsotopeIndex(apex_lmt_in_this_cs->getIsotopeIndex() + tmp_iso_idx);

        final_candidate_fg.push_back(*new_mt);
      }

      if (!rescoreFeatureGroup_(final_candidate_fg)) // don't merge when it failed to exceed filtering threshold
      {
        out_feature.push_back(*candidate_fg);
        // remove it from features
        in_features.erase(candidate_fg);
        continue;
      }

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
    }
    this->endProgress();

    in_features.swap(out_feature);
  }

  void FLASHDeconvQuant::buildMassTraceGroups_(std::vector<LogMassTrace> &log_mtraces, std::vector<FeatureGroup>& features)
  {
    /// group mass traces to spectrum
    std::vector<std::pair<double, LogMassTrace*>> mt_rt_starts;
    std::vector<std::pair<double, LogMassTrace*>> mt_rt_ends;
    mt_rt_starts.reserve(log_mtraces.size());
    mt_rt_ends.reserve(log_mtraces.size());
    int counter = 0;

    // collect rt information from mtraces to generate spectrum
    double min_fwhm_length = numeric_limits<double>::max();
    for (auto &lmt : log_mtraces)
    {
      mt_rt_starts.push_back(std::make_pair(lmt.getFwhmStart(), &lmt));
      mt_rt_ends.push_back(std::make_pair(lmt.getFwhmEnd(), &lmt));
      if( lmt.getMassTrace()->getFWHM() < min_fwhm_length)
      {
        min_fwhm_length = lmt.getMassTrace()->getFWHM();
      }
    }

    if (min_fwhm_length > rt_window_)
    {
      rt_window_ = min_fwhm_length;
    }
    OPENMS_LOG_INFO << "spectrum rt window : " << std::to_string(rt_window_) << endl;

    // sorting mass traces in rt
    std::sort(mt_rt_starts.begin(), mt_rt_starts.end());
    std::sort(mt_rt_ends.begin(), mt_rt_ends.end());

    std::vector<std::pair<double, LogMassTrace*>>::const_iterator rt_s_iter = mt_rt_starts.begin();
    std::vector<std::pair<double, LogMassTrace*>>::const_iterator rt_e_iter = mt_rt_ends.begin();
    auto end_of_iter = mt_rt_starts.end();
    double end_of_current_rt_window = mt_rt_starts[0].first;
    double last_rt = mt_rt_ends[mt_rt_ends.size()-1].first;
    std::vector<LogMassTrace*> local_traces;
    local_traces.reserve(log_mtraces.size());

    int possible_spec_size = int((mt_rt_starts[mt_rt_starts.size()-1].first - end_of_current_rt_window)/rt_window_);
    this->startProgress(0, possible_spec_size, "assembling mass traces to features");

    while(rt_s_iter != end_of_iter && end_of_current_rt_window < last_rt)
    {
      this->setProgress(counter);

      // initial rt binning is 1 sec (for generating spectrum)
      end_of_current_rt_window += rt_window_;

      // add mass traces within rt range
      bool is_new_mt_added = false;
      for(;rt_s_iter != end_of_iter && rt_s_iter->first <= end_of_current_rt_window; ++rt_s_iter)
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
      for(;rt_e_iter != mt_rt_ends.end() && rt_e_iter->first < end_of_current_rt_window ; ++rt_e_iter)
      {
        local_traces.erase(std::remove_if(local_traces.begin(), local_traces.end(),
                                          [&rt_e_iter](auto const& p){return rt_e_iter->second==p; }));
      }
      if(local_traces.empty())
      {
        continue;
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

      features.insert(features.end(), local_fgroup.begin(), local_fgroup.end());
    }

    this->endProgress();
    OPENMS_LOG_INFO << "# generated spec from mass traces : " << counter << endl;
    OPENMS_LOG_INFO << "# generated feature groups from mass traces : " << features.size() << endl;
  }

  void printQuantValuesPerFeature(FeatureGroup fg, String outPath)
  {
    ofstream out;
    out.open(outPath, ios::out);

    out << "centroid_mz\tcs\tiso_index\tintensity\tpeak_area\tsmoothed_peak_area\n";

    for (auto& mt : fg)
    {
      out << std::to_string(mt.getCentroidMz()) << "\t"
          << mt.getCharge() << "\t"
          << mt.getIsotopeIndex() << "\t"
          << mt.getIntensity() << "\t"
          << mt.getMassTrace()->computePeakArea() << "\t"
          << mt.getMassTrace()->computeSmoothedPeakArea() << "\n";
    }

    out.close();
  }

  void FLASHDeconvQuant::resolveSharedMassTraces_simple(std::vector<FeatureGroup> &fgroups,
                                                        std::vector<std::vector<Size>> &shared_m_traces,
                                                        std::vector<MassTrace> &input_mtraces) const
  {
    /// test writing for shared mass traces
    String s_out = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "_sharedMTs.tsv";
    printSharedMassTraces(shared_m_traces, s_out, fgroups, input_mtraces);

    /// assign intensity of shared mass traces to the "feature"(in FeatureGroup) with weights (w: based on apex inty)
    for (Size mt_idx = 0; mt_idx < shared_m_traces.size(); ++mt_idx)
    {
      if (shared_m_traces[mt_idx].size() < 2) // this mass trace is not shared
      {
        continue;
      }

      // find two features to assign shared mass trace
      std::vector<std::pair<double, std::pair<Size, Size>>> feat_quant_pair; // first: apex feature(w/ charge) intensity, second : (f_idx, lmt_idx)
      for (auto& f_idx : shared_m_traces[mt_idx])
      {
        const auto& feature = fgroups[f_idx];

        // find charge of shared mass traces
        int feat_cs;
        Size lmt_idx_in_feat;
        for (Size lmt_idx = 0; lmt_idx < feature.size(); ++lmt_idx)
        {
          auto &lmt = feature[lmt_idx];
          if (mt_idx != lmt.getTraceIndex()) continue;
          feat_cs = lmt.getCharge();
          lmt_idx_in_feat = lmt_idx;
          break;
        }

        // find max quant of feature where shared mass trace belongs (shared cannot be the max)
        double max_quant = .0;
        for (Size lmt_idx = 0; lmt_idx < feature.size(); ++lmt_idx)
        {
          auto &lmt = feature[lmt_idx];
          if (mt_idx == lmt.getTraceIndex()) continue; // if this lmt is the shared one
          if (feat_cs != lmt.getCharge()) continue;

          if (lmt.getIntensity() > max_quant)
          {
            max_quant = lmt.getIntensity();
          }
        }
        feat_quant_pair.push_back(std::make_pair(max_quant, std::make_pair(f_idx, lmt_idx_in_feat)));
      }
//      std::sort(feat_quant_pair.rbegin(), feat_quant_pair.rend());

      // use top 2 featuregroups to assign shared mass traces
//      double weight_1 = feat_quant_pair[0].first/(feat_quant_pair[0].first + feat_quant_pair[1].first);
//      double weight_2 = feat_quant_pair[1].first/(feat_quant_pair[0].first + feat_quant_pair[1].first);
//
//      auto &fg1_lmt = fgroups[feat_quant_pair[0].second.first][feat_quant_pair[0].second.second];
//      auto &fg2_lmt = fgroups[feat_quant_pair[1].second.first][feat_quant_pair[1].second.second];
//
//      fg1_lmt.setIntensity( fg1_lmt.getIntensity() * weight_1 );
//      fg2_lmt.setIntensity( fg2_lmt.getIntensity() * weight_2 );

      double sum_of_quant = std::accumulate(feat_quant_pair.begin(), feat_quant_pair.end(), 0.0,
                                            [](auto &a, auto &b){return a + b.first;});

      for (auto& q_pair : feat_quant_pair)
      {
        double weight = q_pair.first/sum_of_quant;
        auto &fg_lmt = fgroups[q_pair.second.first][q_pair.second.second];
        fg_lmt.setIntensity(fg_lmt.getIntensity() * weight);
      }
//
//      // remove shared mass traces from "the other" features
//      for (Size i = 2; i < feat_quant_pair.size(); ++i)
//      {
//        Size f_idx = feat_quant_pair[i].second.first;
//        auto& feature = fgroups[f_idx];
//        for (auto lmt_iter = feature.begin(); lmt_iter != feature.end(); ++lmt_iter)
//        {
//          if (mt_idx != lmt_iter->getTraceIndex()) continue;
//          feature.erase(lmt_iter);
//          break;
//        }
//      }
    }

    /// re-calcualte quants
    for (auto& f : fgroups)
    {
      f.updateIntensity();
    }

    /// test writing for all mass traces in FeatureGroup
    String out_path = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "groups.tsv";
    ofstream out;
    out.open(out_path, ios::out);
    out << "mono_mass\tcharge\tisotope_index\tquant_value\tshared\tcentroid_mzs\trts\tmzs\tintys\n"; // header
    writeMassTracesOfFeatureGroup(fgroups, shared_m_traces, out);
    out.close();
  }


  void FLASHDeconvQuant::resolveSharedMassTraces(std::vector<FeatureGroup> &fgroups,
                                                 std::vector<std::vector<Size>> &shared_m_traces,
                                                 std::vector<MassTrace> &input_mtraces) const
  {
    /// generate array with indices of shared mass traces
    for (Size fg_index = 0; fg_index<fgroups.size(); ++fg_index)
    {
      for (auto& mt_i : fgroups[fg_index].getTraceIndices())
      {
        shared_m_traces[mt_i].push_back(fg_index);
      }
    }

    /// test writing for shared mass traces
    String s_out = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "_sharedMTs.tsv";
    printSharedMassTraces(shared_m_traces, s_out, fgroups, input_mtraces);


    /// assign shared mass traces to the "feature"(in FeatureGroup) with higher isotopeCosineScorePerCharge
    for (Size mt_idx = 0; mt_idx < shared_m_traces.size(); ++mt_idx)
    {
      if (shared_m_traces[mt_idx].size() < 2) // this mass trace is not shared
      {
        continue;
      }

      double max_score = 0;
      Size feat_idx_with_max_score = 0;

      // find the feature to assign shared mass trace
      for (auto& f_idx : shared_m_traces[mt_idx])
      {
        const auto& feature = fgroups[f_idx];
        if(feature.getIsotopeCosine() == 0) // if this featureGroup is rejected in other competition
        {
          continue;
        }
        for (Size lmt_idx = 0; lmt_idx < feature.size(); ++lmt_idx)
        {
          auto &lmt = feature[lmt_idx];
          if (mt_idx != lmt.getTraceIndex()) continue;

          double tmp_score = feature.getIsotopeCosineOfCharge(lmt.getCharge());
          if (tmp_score > max_score)
          {
            max_score = tmp_score;
            feat_idx_with_max_score = f_idx;
          }
          break;
        }
      }

      // remove shared mass traces from "the other" features
      for (auto& f_idx : shared_m_traces[mt_idx])
      {
        if (f_idx != feat_idx_with_max_score)
        {
          auto& feature = fgroups[f_idx];
          for (auto lmt_iter = feature.begin(); lmt_iter != feature.end(); ++lmt_iter)
          {
            if (mt_idx != lmt_iter->getTraceIndex()) continue;
            feature.erase(lmt_iter);
            break;
          }
        }
        if(!rescoreFeatureGroup_(fgroups[f_idx]))
        {
          fgroups[f_idx].setIsotopeCosine(0.0);
        }
      }
    }

    // rescore
    std::vector<FeatureGroup> tmpFGroups;
    tmpFGroups.swap(fgroups);
    fgroups.reserve(tmpFGroups.size());
    for (auto& f : tmpFGroups)
    {
      if(!rescoreFeatureGroup_(f))
      {
        continue;
      }
      fgroups.push_back(f);
    }

    /// test writing for all mass traces in FeatureGroup
    String out_path = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "groups.tsv";
    ofstream out;
    out.open(out_path, ios::out);
    out << "mono_mass\tcharge\tisotope_index\tquant_value\tshared\tcentroid_mzs\trts\tmzs\tintys\n"; // header
    writeMassTracesOfFeatureGroup(fgroups, shared_m_traces, out);
    out.close();

    // print per feature info
    /// test writing
//    for (Size i =0; i< fgroups.size(); ++i)
//    {
//      auto& f = fgroups[i];
//      String mass = int(f.getMonoisotopicMass());
//      String path_for_f = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "_" + mass + "_" + i + "_mts.tsv";
//      printQuantValuesPerFeature(f, path_for_f);
//    }

  }

  /// cluster FeatureGroups with shared mass traces. If not, report as output
  void FLASHDeconvQuant::clusterFeatureGroups_(std::vector<FeatureGroup>& fgroups, std::vector<MassTrace>& input_mtraces) const
  {
    // *********************************************************** //
    // Step 1 preparation for hypergraph : collect feature idx with shared mass traces
    // *********************************************************** //
    std::vector<std::vector<Size>> shared_m_traces(input_mtraces.size(), std::vector<Size>());
    for (Size fg_index = 0; fg_index<fgroups.size(); ++fg_index)
    {
      for (auto& mt_i : fgroups[fg_index].getTraceIndices())
      {
        shared_m_traces[mt_i].push_back(fg_index);
      }
    }

    // print per feature info
//    /// test writing
//    for (Size i =0; i< fgroups.size(); ++i)
//    {
//      auto& f = fgroups[i];
//      String mass = int(f.getMonoisotopicMass());
//      String path_for_f = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "_" + mass + "_" + i + "_mts.tsv";
//      printQuantValuesPerFeature(f, path_for_f);
//    }
//
    //// test writing
//    String out_path = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "groups.tsv";
    ofstream out;
//    out.open(out_path, ios::out);
//
//    // header
////    out << "feature_label\tcs\tscore\tquant\tbounding_box_pos\tbounding_box_width\tbounding_box_height\t"
////           "iso_position\tmasstrace_centroid_rts\tmasstrace_centroid_mzs\tclustername\n";
//    out << "mono_mass\tclustername\tshared\tcentroid_mzs\trts\tmzs\tintys\n";
//
//    String s_out = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "_sharedMTs.tsv";
//    printSharedMassTraces(shared_m_traces, s_out, fgroups, input_mtraces);

    /// test writing ends here

    // *********************************************************** //
    // Step 1 constructing hypergraph from featuregroups
    //        node = mass traces
    //        hyperedge = feature groups
    // *********************************************************** //
    Size num_nodes = shared_m_traces.size();
    std::vector<bool> bfs_visited;
    bfs_visited.resize(num_nodes, false);
    std::queue<Size> bfs_queue;
    Size search_pos = 0; // keeping track of mass trace index to look for seed

    std::vector<FeatureGroup> out_features;
    out_features.reserve(fgroups.size());
//    std::vector<Size> out_feature_idx;
//    out_feature_idx.reserve(fgroups.size());

//    this->startProgress(0, shared_m_traces.size(), "clustering features based on the shared mass traces");
//    Size progress = 0;
    Size cluster_counter = 0;
    // BFS
    while (true)
    {
//      this->setProgress(progress);
//      ++progress;
      // finding a seed 'shared_mass_trace' to start with (for constructing a cluster)
      bool finished = true;
      for (Size i = search_pos; i < num_nodes; ++i)
      {
        if (!bfs_visited[i])
        {
          // check if this mass_trace is used to any FeatureGroup
          if (shared_m_traces[i].size()==0)
          {
            bfs_visited[i] = true;
            continue;
          }

          bfs_queue.push(i);
          bfs_visited[i] = true;
          finished = false;
          search_pos = i + 1;
          break;
        }
      }
      if (finished) // if no possible seed is left
        break;

      std::set<Size> fg_indices_in_current_cluster;

      while (!bfs_queue.empty())
      {
        Size i = bfs_queue.front(); // index of seed
        bfs_queue.pop();

        // get FeatureGroup indices sharing this seed
        for (vector<Size>::const_iterator it = shared_m_traces[i].begin();
             it != shared_m_traces[i].end();
             ++it)
        {
          // if this FeatureGroup was visited before
          if (fg_indices_in_current_cluster.find(*it) != fg_indices_in_current_cluster.end())
          {
            continue;
          }

          fg_indices_in_current_cluster.insert(*it);
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

      // this feature is not sharing any mass traces with others
      if (fg_indices_in_current_cluster.size() == 1){
        out_features.push_back(fgroups[*(fg_indices_in_current_cluster.begin())]);
//        out_feature_idx.push_back(*(fg_indices_in_current_cluster.begin()));
        continue;
      }

      // resolve the conflict among feature groups
      cluster_counter++;
      String cluster_name = "cluster" + to_string(cluster_counter);
      resolveConflictInCluster_(fgroups, input_mtraces, shared_m_traces, fg_indices_in_current_cluster, out_features, out, cluster_name);
    }
//    this->endProgress();

    out_features.shrink_to_fit();
    std::swap(out_features, fgroups);
//    out_feature_idx.shrink_to_fit();
//    out.close();
    OPENMS_LOG_INFO << "#final feature groups: " << out_features.size() << endl;
    OPENMS_LOG_INFO << "#cluster (shared) :" << cluster_counter << endl;
  }

  void FLASHDeconvQuant::writeMassTracesOfFeatureGroup(const std::vector<FeatureGroup> &feature_groups,
                                                       const std::vector<std::vector<Size>> &shared_m_traces_indices,
                                                       ofstream &out) const
  {
    OPENMS_LOG_INFO << "----- # FeatureGroup : " << feature_groups.size() << " -----" << endl;
    for (const auto& feat : feature_groups)
    {
      double quant(0.0);
      auto mt_idxs = feat.getTraceIndices();

      for (const auto &mt : feat)
      {

        // check if current mt is shared with other features
        bool isShared = false;
        if (shared_m_traces_indices[mt.getTraceIndex()].size() > 1)
        {
          isShared = true;
        }

        stringstream rts;
        stringstream mzs;
        stringstream intys;

        for (const auto &peak : *(mt.getMassTrace()))
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

        out << std::to_string(feat.getMonoisotopicMass()) << "\t"
            << mt.getCharge() << "\t"
            << mt.getIsotopeIndex() << "\t"
            << std::to_string(mt.getIntensity()) << "\t"
            << isShared << "\t" << std::to_string(mt.getCentroidMz()) << "\t"
            << peaks + "\n";
      }
    }
  }

  void FLASHDeconvQuant::resolveConflictInCluster_(std::vector<FeatureGroup>& feature_groups,
                                                   const std::vector<MassTrace> & input_masstraces,
                                                   const std::vector<std::vector<Size> >& shared_m_traces_indices,
                                                   const std::set<Size>& fg_indices_in_this_cluster,
                                                   std::vector<FeatureGroup>& out_featuregroups,
                                                   ofstream& out,
                                                   String& cluster_name) const
  {
    OPENMS_LOG_DEBUG << "-------- " << cluster_name <<  " -------- \n";

    /// conflict resolution is done in feature level (not feature group level)

    // prepare nodes for hypergraph
    std::vector<std::vector<Size>> mt_and_feature_idx(shared_m_traces_indices.size(), std::vector<Size>());
//                                                                      std::vector<std::pair<Size, Size>>()); // pair <fg_idx_in_cluster, charge_idx>
    std::vector<FeatureElement> feature_candidates; // index : (fg_indices_in_this_cluster+1)*charge_idx
    for (auto& fg_idx : fg_indices_in_this_cluster)
    {
      auto &fg = feature_groups[fg_idx];
      OPENMS_LOG_DEBUG << fg_idx << "\t" << fg.getMonoisotopicMass() << endl;

      auto cs_range = fg.getChargeRange();
      int min_cs = std::get<0>(cs_range);
      Size candidate_size = std::get<1>(cs_range) - min_cs + 1;
      Size feature_idx_starts = feature_candidates.size();
      for (Size tmp_idx = 0; tmp_idx < candidate_size; ++tmp_idx)
      {
        FeatureElement tmp;
        tmp.mass_traces.reserve(fg.size());
        tmp.mass_trace_indices.reserve(fg.size());
        tmp.charge = min_cs + tmp_idx;
        tmp.feature_group_index = fg_idx;
        feature_candidates.push_back(tmp);
      }

      for (auto lmt_iter = fg.begin(); lmt_iter != fg.end(); ++lmt_iter)
      {
        Size f_idx = feature_idx_starts + (lmt_iter->getCharge()-min_cs);
        mt_and_feature_idx[lmt_iter->getTraceIndex()].push_back(f_idx);
        feature_candidates[f_idx].mass_traces.push_back(&(*lmt_iter));
        feature_candidates[f_idx].mass_trace_indices.push_back(lmt_iter->getTraceIndex());
      }

    }

    Size num_nodes = shared_m_traces_indices.size();
    std::vector<bool> bfs_visited;
    bfs_visited.resize(num_nodes, false);
    std::queue<Size> bfs_queue;
    Size search_pos = 0; // keeping track of mass trace index to look for seed

    // BFS
    while (true)
    {
      bool finished = true;

      for (Size i = search_pos; i < num_nodes; ++i)
      {
        if (!bfs_visited[i])
        {
          // check if this mass trace is used in any feature
          if (mt_and_feature_idx[i].size() == 0)
          {
            bfs_visited[i] = true;
            continue;
          }

          bfs_queue.push(i);
          bfs_visited[i] = true;
          finished = false;
          search_pos = i + 1;
          break;
        }
      }
      if (finished) // if no possible seed is left
        break;

      std::set<Size> feature_idx_in_current_conflict_region;
      std::set<Size> conflicting_mt_indices;

      while(!bfs_queue.empty())
      {
        Size i = bfs_queue.front(); // index of seed
        bfs_queue.pop();

        // get feature indices sharing this seed
        for (std::vector<Size>::const_iterator it = mt_and_feature_idx[i].begin();
             it != mt_and_feature_idx[i].end();
             ++it)
        {
          // if this feature was visited before
          if (feature_idx_in_current_conflict_region.find(*it) != feature_idx_in_current_conflict_region.end())
          {
            continue;
          }

          feature_idx_in_current_conflict_region.insert(*it);
          FeatureElement &selected_feature = feature_candidates[*it];

          for (const auto &mt_index : selected_feature.mass_trace_indices)
          {
            if (!bfs_visited[mt_index])
            {
              bfs_queue.push(mt_index);
              bfs_visited[mt_index] = true;
            }
            else if (mt_and_feature_idx[mt_index].size() > 1){
              // check if this visited mass trace is the shared one (because the seed mass trace can end up here as well)
              conflicting_mt_indices.insert(mt_index);
            }
          }
        }
      }

      // this feature is not sharing any mass traces with others
      if (feature_idx_in_current_conflict_region.size() == 1){
        continue;
      }

      // collect conflicting mass traces (not LogMassTrace, originals)
      std::vector<const MassTrace*> conflicting_mts;
      for (auto &mt_idx : conflicting_mt_indices)
      {
        auto i = input_masstraces.begin() + mt_idx;
        conflicting_mts.push_back(&(*i));
      }

      // set isotope probabilities for feature_candidates && check if mass traces are only shared ones
      vector<Size> features_not_for_resolution;
      for (auto &feat_idx : feature_idx_in_current_conflict_region)
      {
        auto &feat = feature_candidates[feat_idx];
        auto &feat_group = feature_groups[feat.feature_group_index];
        auto tmp_iso = iso_model_.get(feat_group.getMonoisotopicMass());
        feat.isotope_probabilities.reserve(feat.mass_traces.size());
        for (auto &lmt : feat.mass_traces)
        {
          // if isotope index of lmt exceed tmp_iso length, give 0
          auto tmp_iso_idx = lmt->getIsotopeIndex();
          if (tmp_iso_idx < tmp_iso.size())
          {
            feat.isotope_probabilities.push_back(tmp_iso[tmp_iso_idx].getIntensity());
          }
          else
          {
            feat.isotope_probabilities.push_back(0.0);
          }
        }

        // check if this feature has enough mass traces to get a component
        if (feat.mass_traces.size() > conflicting_mts.size())
        {
          continue;
        }

        // check if this feature is composed of only shared mass traces
        bool composed_of_only_conflicts = true;
        for (auto& mt_idx : feat.mass_trace_indices)
        {
          if (shared_m_traces_indices[mt_idx].size() == 1)
          {
            composed_of_only_conflicts = false;
            break;
          }
        }
        if (!composed_of_only_conflicts)
        {
          continue;
        }

        // get the most abundant mass trace from FeatureGroup (except recruited ones)
        LogMassTrace* most_abundant_mt = nullptr;
        getMostAbundantMassTraceFromFeatureGroup(feat_group, feat.charge, most_abundant_mt, shared_m_traces_indices);
        if (most_abundant_mt == nullptr)
        {
          features_not_for_resolution.push_back(feat_idx);
          continue;
        }

        // get the most abundant mass trace from FeatureGroup (except recruited ones)
        feat.mass_traces.push_back(most_abundant_mt);
        feat.mass_trace_indices.push_back(most_abundant_mt->getTraceIndex());
        feat.isotope_probabilities.push_back(0.2);
      }

      // if any feature is not eligible for resolution
      if(features_not_for_resolution.size() > 0)
      {
        for (Size& idx : features_not_for_resolution)
        {
          feature_idx_in_current_conflict_region.erase(idx);
          for (auto& lmt : feature_candidates[idx].mass_traces)
          {
            lmt->setIntensity(0.0);
          }
        }

        if (feature_idx_in_current_conflict_region.size() < 2)
        {
          continue;
        }

        // check if any shared mass trace is not shared anymore
        vector<Size> is_mt_included(conflicting_mts.size(), 0); // element is equal to the size of feature_vec when it's shared
        for (auto &feat_idx : feature_idx_in_current_conflict_region)
        {
          auto &lmts = feature_candidates[feat_idx].mass_traces;
          for (auto &lmt : lmts)
          {
            for (Size mt_i = 0; mt_i < conflicting_mts.size(); ++mt_i)
            {
              if(conflicting_mts[mt_i] == lmt->getMassTrace())
              {
                is_mt_included[mt_i] += 1;
                break;
              }
            }
          }
        }
        for (int mt_i = conflicting_mts.size()-1; mt_i >= 0; --mt_i)
        {
          if (is_mt_included[mt_i] != feature_idx_in_current_conflict_region.size())
          {
            conflicting_mts.erase(conflicting_mts.begin() + mt_i);
          }
        }
      }

      // convert feature_idx_in_current_conflict_region from set to vector (for accessing with indices later)
      std::vector<Size> feature_indices_vec;
      feature_indices_vec.reserve(feature_idx_in_current_conflict_region.size());
      for (auto& i : feature_idx_in_current_conflict_region) feature_indices_vec.push_back(i);

      // resolve this conflict region
      resolveConflictRegion_(feature_candidates, feature_indices_vec, conflicting_mts);

      OPENMS_LOG_INFO << "---------conflicting masses----------" << endl;
      for (auto& i : feature_indices_vec)
      {
        auto fg_index = feature_candidates[i].feature_group_index;
        OPENMS_LOG_INFO << to_string(feature_groups[fg_index].getMonoisotopicMass()) << endl;
      }
    }

    // update feature group quantities
    for (auto& idx : fg_indices_in_this_cluster)
    {
      auto &fgroup = feature_groups[idx];

      // remove 0 intensity mass traces
      auto fg_iter = fgroup.begin();
      while (fg_iter != fgroup.end())
      {
        if (fg_iter->getIntensity() == 0)
        {
          fg_iter = fgroup.erase(fg_iter);
        }
        else {
          ++fg_iter;
        }
      }

      // check if enough number of mass traces are left
      if (fgroup.size() < min_nr_mtraces_)
      {
        continue;
      }

      // check if the FeatureGroup still contains more than three charges
      std::set<int> charges;
      for (auto &lmt : fgroup)
      {
        charges.insert(lmt.getCharge());
      }
      if (charges.size() < 3)
      {
        continue;
      }

      // update intensity
      feature_groups[idx].updateIntensity();
      out_featuregroups.push_back(feature_groups[idx]);
    }

  }

   void FLASHDeconvQuant::resolveConflictRegion_(std::vector<FeatureElement> &features,
                              const std::vector<Size> &feature_idx_in_current_conflict_region,
                              const std::vector<const MassTrace*> &conflicting_mts) const
  {
    /// Prepare Components per features (excluding conflicting mts)
    std::vector<std::vector<double>> components;
    Matrix<int> pointer_to_components; // row : index of conflicting_mts , column : index of feature_idx_in_current_conflict_region
    pointer_to_components.resize(conflicting_mts.size(), feature_idx_in_current_conflict_region.size(), -1);
    std::vector<std::vector<int>> conflicting_mt_idx_lookup_vec; // linker between mt_idx in feature & conflicting_mts
    for (Size i_of_f = 0; i_of_f < feature_idx_in_current_conflict_region.size(); ++i_of_f)
    {
      auto &tmp_feat = features[feature_idx_in_current_conflict_region[i_of_f]];

      std::vector<int> lookup_list_for_mt_idx; // to get theoretical intensities (isotope probabilities)
      lookup_list_for_mt_idx.resize(conflicting_mts.size(), -1);
      Size mass_traces_size = tmp_feat.mass_traces.size();

      // preparation for ElutionModelFit
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces traces_for_fitting;
      traces_for_fitting.reserve(mass_traces_size);
      traces_for_fitting.max_trace = 0;
      vector<Peak1D> peaks;
      // reserve space once, to avoid copying and invalidating pointers:
      peaks.reserve(tmp_feat.getPeakSizes());

      // get theoretical information from non-shared mass traces
      for(Size idx = 0; idx < mass_traces_size; ++idx)
      {
        // ignore this mt if it's the conflicting one
        auto tmp_iter = std::find(conflicting_mts.begin(), conflicting_mts.end(), tmp_feat.mass_traces[idx]->getMassTrace());
        if (tmp_iter != conflicting_mts.end())
        {
          lookup_list_for_mt_idx[std::distance(conflicting_mts.begin(), tmp_iter)] = idx;
          continue;
        }

        auto mass_trace_ptr = tmp_feat.mass_traces[idx]->getMassTrace();

        // prepare for fitting
        FeatureFinderAlgorithmPickedHelperStructs::MassTrace tmp_mtrace;
        Size m_trace_size = (*mass_trace_ptr).getSize();
        tmp_mtrace.peaks.reserve(m_trace_size);
        for (auto &p_2d : (*mass_trace_ptr))
        {
          auto tmp_inty = p_2d.getIntensity();
          if (tmp_inty > 0.0) // only use non-zero intensities for fitting
          {
            Peak1D peak;
            peak.setMZ(p_2d.getMZ());
            peak.setIntensity(tmp_inty);
            peaks.push_back(peak);
            tmp_mtrace.peaks.emplace_back(p_2d.getRT(), &peaks.back());
          }
        }
        tmp_mtrace.updateMaximum();
        tmp_mtrace.theoretical_int = tmp_feat.isotope_probabilities[idx];
        traces_for_fitting.push_back(tmp_mtrace);
      }

      if (traces_for_fitting.size() == 1 && traces_for_fitting[0].peaks.size() < 4)
      {
        while (traces_for_fitting[0].peaks.size() < 4)
        {
          auto& tmp_trace = traces_for_fitting[0];
          Peak1D peak;
          peak.setMZ(tmp_trace.peaks[0].second->getMZ());
          peak.setIntensity(0.0);
          peaks.push_back(peak);
          double region_start = tmp_trace.peaks.front().second->getPos();
          double region_end = tmp_trace.peaks.back().second->getPos();
          double offset = 0.2 * (region_end - region_start);
          tmp_trace.peaks.emplace_back(region_end + offset, &peaks.back());
        }
      }

      // ElutionModelFit
      EGHTraceFitter* fitter = new EGHTraceFitter();
      // TODO : is this necessary?
      Param params = fitter->getDefaults();
      params.setValue("weighted", "true");

      fitter->setParameters(params);
      runElutionModelFit_(traces_for_fitting, fitter);

      // store Component information
      for (Size row = 0; row < conflicting_mts.size(); ++row)
      {
        std::vector<double> component;

        int index_in_feature = lookup_list_for_mt_idx[row];
        // this mt is not included in current feature
        if (index_in_feature < 0)
        {
          continue;
        }

        double intensity_ratio = tmp_feat.isotope_probabilities[index_in_feature];
        for (auto &peak : *(conflicting_mts[row]))
        {
          double rt = peak.getRT();
          double theo_intensity = intensity_ratio * fitter->getValue(rt);

          component.push_back(theo_intensity);
        }

        pointer_to_components.setValue(row, i_of_f, components.size());
        components.push_back(component);
      }
      conflicting_mt_idx_lookup_vec.push_back(lookup_list_for_mt_idx);
    }

    // reconstruction of observed XICs (per conflicting mt)
    for (Size row = 0; row < conflicting_mts.size(); ++row)
    {
      auto components_indices = pointer_to_components.row(row);
      Size column_size = components_indices.size() - std::count(components_indices.begin(), components_indices.end(), -1);

      Size mt_size = conflicting_mts[row]->getSize();
      auto mt_iter = conflicting_mts[row]->begin();
      Matrix<double> obs;
      obs.resize(mt_size, 1);
      for (Size i = 0; i < mt_size; ++i)
      {
        obs.setValue(i, 0, mt_iter->getIntensity());
        mt_iter++;
      }

      Matrix<double> theo_matrix;
      theo_matrix.resize(mt_size, column_size);
      Size col = 0;
      for (Size comp_idx = 0; comp_idx < components_indices.size(); ++comp_idx)
      {
        // skipping no-feature column
        if (components_indices[comp_idx] == -1) continue;

        auto &temp_comp = components[pointer_to_components.getValue(row, comp_idx)];
        for (Size tmp_r = 0; tmp_r < temp_comp.size(); ++tmp_r)
        {
          theo_matrix.setValue(tmp_r, col, temp_comp[tmp_r]);
        }
        col++;
      }

      Matrix<double> out_quant;
      out_quant.resize(column_size, 1);
      NonNegativeLeastSquaresSolver::solve(theo_matrix, obs, out_quant);

      OPENMS_LOG_INFO << "-------------------------" << endl;
      OPENMS_LOG_INFO <<  "observed m/z : " <<  std::to_string(conflicting_mts[row]->getCentroidMZ()) << endl;
      for (auto& peak : (*conflicting_mts[row]))
      {
        OPENMS_LOG_INFO << std::to_string(peak.getRT()) << "\t" << std::to_string(peak.getIntensity()) << "\n";
      }

      // update lmt & feature intensity
      col = 0;
      for (Size i_of_f = 0; i_of_f < components_indices.size(); ++i_of_f)
      {
        // skipping no-feature column
        if (components_indices[i_of_f] == -1)
          continue;

        auto &feat = features[feature_idx_in_current_conflict_region[i_of_f]];
        int &lmt_index = conflicting_mt_idx_lookup_vec[i_of_f][row];
        auto lmt_ptr = feat.mass_traces[lmt_index];

        auto temp_component = components[pointer_to_components.getValue(row, i_of_f)];
        double theo_intensity = std::accumulate(temp_component.begin(), temp_component.end(), 0.0);

        lmt_ptr->setIntensity(out_quant.getValue(col, 0) * theo_intensity);

        OPENMS_LOG_INFO << "--- theo[" << col << "]--- "<< std::to_string(out_quant.getValue(col, 0)) << "\t" << std::to_string(lmt_ptr->getIntensity()) << "\n";
        for (auto& m : theo_matrix.col(col))
        {
          OPENMS_LOG_INFO << to_string(m) << "\t";
        }
        OPENMS_LOG_INFO << "\n";
        col++;
      }
    }
    components.clear();
  }

  // from ElutionModelFitter. TODO : this is test purpose. use the existing method insteaad.
  void FLASHDeconvQuant::runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces,
                                             EGHTraceFitter* fitter) const
  {
    // find the trace with maximal intensity:
    Size max_trace = 0;
    double max_intensity = 0;
    for (Size i = 0; i < m_traces.size(); ++i)
    {
      if (m_traces[i].max_peak->getIntensity() > max_intensity)
      {
        max_trace = i;
        max_intensity = m_traces[i].max_peak->getIntensity();
      }
    }
    m_traces.max_trace = max_trace;
    m_traces.baseline = 0.0;

    // TODO : add zeros?


    // Fitting
    bool fit_success = true;
    try
    {
      fitter->fit(m_traces);
    }
    catch (Exception::UnableToFit& except)
    {
      OPENMS_LOG_ERROR << "Error fitting model to feature '"
                       << except.getName()
                       << " - " << except.getMessage() << endl;
      fit_success = false;
    }

    // record model parameters:
//    double center = fitter->getCenter(), height = fitter->getHeight();
//    double sigma = fitter->getSigma();
//    double tau = fitter->getTau();
//    double width = sigma * 0.6266571 + abs(tau);
//    double asymmetry = abs(tau) / sigma;
//
//    double lower_rt_bound = fitter->getLowerRTBound();
//    double upper_rt_bound = fitter->getUpperRTBound();
  }

  void FLASHDeconvQuant::getMostAbundantMassTraceFromFeatureGroup(const FeatureGroup &fgroup,
                                                                  const int &skip_this_charge,
                                                                  LogMassTrace* &most_abundant_mt_ptr,
                                                                  const std::vector<std::vector<Size>>& shared_m_traces) const
  {
    // get intensities
    double max_intensity = .0; // maximum mt intensity in this feature group

    for (auto lmt = fgroup.begin(); lmt != fgroup.end(); ++lmt)
    {
      if (skip_this_charge > 0 && lmt->getCharge() == skip_this_charge)
      {
        continue;
      }

      // check if this mt is shared with other FeatureGroup
      if (shared_m_traces[lmt->getTraceIndex()].size() > 1)
      {
        continue;
      }

      if (lmt->getIntensity() > max_intensity)
      {
        max_intensity = lmt->getIntensity();
        most_abundant_mt_ptr = (LogMassTrace *) &(*lmt);
      }
    }
  }
}