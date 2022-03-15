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
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  FLASHDeconvQuant::FLASHDeconvQuant() :
      ProgressLogger(),
      DefaultParamHandler("FLASHDeconvQuant")
  {
    defaults_.setValue("local_rt_range", 15.0, "RT range where to look for coeluting mass traces", {"advanced"});
    defaults_.setValue("local_mz_range",
                       6.5,
                       "MZ range where to look for isotopic mass traces",
                       {"advanced"}); // (-> decides size of isotopes =(local_mz_range_ * lowest_charge))
    defaults_.setValue("charge_lower_bound", 5, "Lowest charge state to consider");
    defaults_.setValue("charge_upper_bound", 50, "Highest charge state to consider");
    defaults_.setValue("min_mass", 10000, "minimim mass");
    defaults_.setValue("max_mass", 70000, "maximum mass");
    defaults_.setValue("mz_tol", 20, "ppm tolerance for m/z");

    defaults_.setValue("use_smoothed_intensities",
                       "true",
                       "Use LOWESS intensities instead of raw intensities.",
                       {"advanced"});
    defaults_.setValidStrings("use_smoothed_intensities", {"false", "true"});

    defaultsToParam_();

    this->setLogType(CMD);
  }

  FLASHDeconvQuant::~FLASHDeconvQuant()
  {
  }

  void FLASHDeconvQuant::updateMembers_()
  {
//    local_rt_range_ = (double) param_.getValue("local_rt_range");
    local_mz_range_ = (double) param_.getValue("local_mz_range");

    charge_lower_bound_ = (Size) param_.getValue("charge_lower_bound");
    charge_upper_bound_ = (Size) param_.getValue("charge_upper_bound");
    charge_range_ = charge_upper_bound_ - charge_lower_bound_ + 1;

    min_mass_ = (double) param_.getValue("min_mass");
    max_mass_ = (double) param_.getValue("max_mass");

    mz_tolerance_ = (double) param_.getValue("mz_tol");
    mz_tolerance_ *= 1e-6;

    use_smoothed_intensities_ = param_.getValue("use_smoothed_intensities").toBool();
  }

  Param FLASHDeconvQuant::getFLASHDeconvParams_()
  {
    Param fd_defaults = FLASHDeconvAlgorithm().getDefaults();
    // overwrite algorithm default so we export everything (important for copying back MSstats results)
    fd_defaults.setValue("min_charge", (int) charge_lower_bound_);
    fd_defaults.setValue("max_charge", (int) charge_upper_bound_);
    fd_defaults.setValue("min_mass", min_mass_);
    fd_defaults.setValue("max_mass", max_mass_);
    fd_defaults.setValue("min_isotope_cosine", DoubleList{min_isotope_cosine_, .85});

    fd_defaults.setValue("min_qscore", .0);
    fd_defaults.setValue("tol", DoubleList{10.0, 10.0});
    fd_defaults.setValue("rt_window", 180.0);
    fd_defaults.setValue("min_peaks", IntList{(int) min_nr_mtraces_, 3});//

    // theoretical isotopical model (calculateAveragine)
//    iso_model_ = FLASHDeconvHelperStructs().calculateAveragines(max_mass_, false);

    return fd_defaults;
  }

  void FLASHDeconvQuant::writeFeatureGroupsInFile(std::vector<FeatureGroup> &fgroups,
                                                  std::vector<MassTrace> &input_mtraces)
  {
    // to get "shared" information
    std::vector<std::vector<Size>> shared_m_traces(input_mtraces.size(), std::vector<Size>());
    for (Size fg_index = 0; fg_index < fgroups.size(); ++fg_index)
    {
      for (auto &mt_i: fgroups[fg_index].getTraceIndices())
      {
        shared_m_traces[mt_i].push_back(fg_index);
      }
    }

    //    std::vector<bool> fg_with_shared(fgroups.size(), false);
    //    for (auto &sharing_fgs : shared_m_traces)
    //    {
    //      if(sharing_fgs.size() < 2)
    //        continue;
    //
    //      for (auto &fg_idx : sharing_fgs)
    //      {
    //        fg_with_shared[fg_idx] = true;
    //      }
    //    }

    std::fstream out_stream;
    OPENMS_LOG_INFO << "writing output..." << outfile_path << endl;
    out_stream.open(outfile_path, std::fstream::out);

    // header
    out_stream << "mono_mass\tmin_charge\tmax_charge\tquant_value\tsummed_inty\t" // most_abundant_charge\t"
                  // "quant_of_most_abundant_charge\tsummed_inty_of_most_abundant_cs\tarea_of_most_abundant_cs\t"
                  // "quant_of_apices_of_cs\tsummed_inty_of_apices_of_cs\t"
                  //           "quant_of_top_three_cs\tsummed_inty_of_top_three_cs\tarea_of_top_three_cs\t"
                  "all_area_under_the_curve\t"
                  // "charge_score\t
                  "iso_cosine\tweighted_iso_score\t"
                  "fwhm_start\tfwhm_end\n"; // fwhm_avg_length\trt_start\trt_end\tcentroid_rt_of_apices\trt_of_apex\t"
    // "iso_length\tperIsotopeIntensity\n"; //	is_shared

    int counter = -1;
    for (auto &fg: fgroups)
    {
      counter++;

      double quant_value = fg.getIntensity();

      // intys
      //      std::vector<double> per_iso_area = std::vector<double>(iso_model_.getMaxIsotopeIndex()+1, .0);
      std::vector<double> per_cs_intys = std::vector<double>(fg.getMaxCharge() + 1, .0);
      //      std::vector<double> per_cs_area = std::vector<double>(std::get<1>(fg.getChargeRange())+1 , .0);
      std::vector<double> per_cs_area_all = std::vector<double>(fg.getMaxCharge() + 1, .0);

      // rt range
      //      double min_rt = std::numeric_limits<double>::max();
      //      double max_rt = 0;

      // getting information while looping through mass traces in featuregroup
      for (auto &lmt: fg)
      {
        auto lmt_ptr = lmt.getMassTrace();
        //        double tmp_min = lmt_ptr->getConvexhull().getBoundingBox().minX();
        //        double tmp_max = lmt_ptr->getConvexhull().getBoundingBox().maxX();

        double int_before = (*lmt_ptr)[0].getIntensity();
        double rt_before = (*lmt_ptr)[0].getRT();

        double tmp_area = 0;
        for (auto &peaks: *lmt_ptr)
        {
          per_cs_intys[lmt.getCharge()] += peaks.getIntensity();
          tmp_area += (int_before + peaks.getIntensity()) / 2 * (peaks.getRT() - rt_before);
          int_before = peaks.getIntensity();
          rt_before = peaks.getRT();
        }
        per_cs_area_all[lmt.getCharge()] += tmp_area;
        //        per_iso_area[lmt.getIsotopeIndex()] += lmt.getIntensity();
        //        per_cs_area[lmt.getCharge()] += lmt.getIntensity();

        //        if (tmp_min < min_rt)
        //          min_rt = tmp_min;
        //        if (tmp_max > max_rt)
        //          max_rt = tmp_max;
      }

      // get centroid rt of Apices
      //      fg.setCentroidRtOfApices();

      // get adundances of CS apices
      //      auto abundances_per_cs = fg.getSummedIntensityOfMostAbundantMTperCS();

      // get most abundant charge
      //      int most_abundant_cs = std::distance(per_cs_intys.begin(), std::max_element(per_cs_intys.begin(), per_cs_intys.end()));

      // intensities
      auto summedIntensities = std::accumulate(per_cs_intys.begin(), per_cs_intys.end(), .0);
      //      double summedInty_most_abundant_cs = *std::max_element(per_cs_intys.begin(), per_cs_intys.end());
      //      double area_most_abundant_cs = *std::max_element(per_cs_area.begin(), per_cs_area.end());
      double area_under_the_curve = std::accumulate(per_cs_area_all.begin(), per_cs_area_all.end(), .0);
      //      double all_area_most_abundant_cs = *std::max_element(per_cs_area_all.begin(), per_cs_area_all.end());

      // top3
      //      std::sort(per_cs_area.begin(), per_cs_area.end(), std::greater<double>());
      //      std::sort(per_cs_intys.begin(), per_cs_intys.end(), std::greater<double>());
      //      std::sort(per_cs_area_all.begin(), per_cs_area_all.end(), std::greater<double>());
      //      double top3_area = per_cs_area[0] + per_cs_area[1] + per_cs_area[2];
      //      double top3_intys = per_cs_intys[0] + per_cs_intys[1] + per_cs_intys[2];
      //      double top3_area_all = per_cs_area_all[0] + per_cs_area_all[1] + per_cs_area_all[2];

      //      auto iso_start = std::distance(per_iso_area.begin(), find_if( per_iso_area.begin(), per_iso_area.end(), [](double x) { return x != 0; }));
      //      auto iso_end = std::distance(per_iso_area.rbegin(), find_if( per_iso_area.rbegin(), per_iso_area.rend(), [](double x) { return x != 0; }));
      //      iso_end = per_iso_area.size() - iso_end;

      //      stringstream iso_ss;
      //      for(int i = iso_start; i < iso_end; ++i)
      //      {
      //        iso_ss << std::to_string(per_iso_area[i]) << ";";
      //      }
      //      std::string iso_ss_str = iso_ss.str();
      //      iso_ss_str.pop_back();

      out_stream << std::to_string(fg.getMonoisotopicMass()) << "\t"
                 << fg.getMinCharge() << "\t"
                 << fg.getMaxCharge() << "\t"
                 << std::to_string(quant_value) << "\t"
                 << std::to_string(summedIntensities) << "\t"
                 //                  << most_abundant_cs << "\t"
                 //                  << std::to_string(area_most_abundant_cs) << "\t"
                 //                  << std::to_string(summedInty_most_abundant_cs) << "\t"
                 //                  << std::to_string(all_area_most_abundant_cs) << "\t"
                 //                  << std::to_string(abundances_per_cs.first) << "\t"
                 //                  << std::to_string(abundances_per_cs.second) << "\t"
                 //                  << std::to_string(top3_area) << "\t"
                 //                  << std::to_string(top3_intys) << "\t"
                 //                  << std::to_string(top3_area_all) << "\t"
                 << std::to_string(area_under_the_curve) << "\t"
                 //                  << std::to_string(fg.getChargeScore()) << "\t"
                 << std::to_string(fg.getIsotopeCosine()) << "\t"
                 << std::to_string(fg.getFeatureGroupScore()) << "\t"
                 << std::to_string(fg.getFwhmRange().first) << "\t"
                 << std::to_string(fg.getFwhmRange().second) << std::endl;
      //                  << std::to_string(fg.getAvgFwhmLength()) << "\t"
      //                  << std::to_string(min_rt) << "\t"
      //                  << std::to_string(max_rt) << "\t"
      //                  << std::to_string(fg.getCentroidRtOfApices()) << "\t"
      //                  << std::to_string(fg.getRtOfApex()) << "\t"
      //                  << (iso_end - iso_start) << "\t"
      //                  << iso_ss_str << std::endl;
      //                  << fg_with_shared[counter] << std::endl;

      out_stream.flush();
    }
    out_stream.close();
    OPENMS_LOG_INFO << "----- output writing done -----" << endl;

    /// only for writing purpose belowe here :
    // writing mass traces (when no resolution is done)
    //    String out_path = outfile_path.substr(0, outfile_path.find_last_of(".")-1) + "groups.tsv";
    //    ofstream groups_out;
    //    groups_out.open(out_path, ios::out);
    //    groups_out << "mono_mass\tcharge\tisotope_index\tquant_value\tshared\tcentroid_mzs\trts\tmzs\tintys\n"; // header
    //    writeMassTracesOfFeatureGroup(fgroups, shared_m_traces, groups_out);
    //    groups_out.close();
    return;
  }

  void printSharedMassTraces(std::vector<std::vector<Size>> &shared_m_traces,
                             String outputPath,
                             std::vector<FeatureGroup> &fgroups,
                             std::vector<MassTrace> &input_mtraces)
  {
    std::fstream out;
    out.open(outputPath, std::fstream::out);
    out
        << "mass_trace_idx\tmass_trace_centroid_mz\tmass_trace_centroid_rt\tfeature_indices\tfeature_masses\tfeature_cs_iso\trts_or_apices\n";

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
      for (auto &f_idx: shared_m_traces[mt_idx])
      {
        auto &feature = fgroups[f_idx];

        feature.setCentroidRtOfApices();
        f_indices << f_idx << ", ";
        f_masses << std::to_string(feature.getMonoisotopicMass()) << ", ";
        f_rts << std::to_string(feature.getCentroidRtOfApices()) << ", ";

        for (auto &mt: feature)
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
//    getFLASHDeconvConsensusResult();

    if (Constants::C13C12_MASSDIFF_U * max_nr_traces_ / charge_lower_bound_ < local_mz_range_)
    {
      local_mz_range_ = Constants::C13C12_MASSDIFF_U * max_nr_traces_ / charge_lower_bound_;
    }

    // initialize input & output
    std::vector<FeatureSeed> input_seeds;
    Size index = 0;
    double sum_intensity = 0;
    for (auto iter=input_mtraces.begin(); iter!=input_mtraces.end(); ++iter, ++index)
    {
      FeatureSeed tmp_seed(*iter);
      tmp_seed.setTraceIndex(index);
      input_seeds.push_back(tmp_seed);
      sum_intensity += tmp_seed.getIntensity();
    }
    // sort input mass traces in RT
    std::sort(input_seeds.begin(), input_seeds.end(), CmpFeatureSeedByRT());
    total_intensity_ = sum_intensity;
    std::vector<FeatureGroup> features;
    features.reserve(input_mtraces.size());

    // run deconvolution
    buildMassTraceGroups_(input_seeds, features);
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
    clusterFeatureGroups_(features, input_mtraces);
    writeFeatureGroupsInFile(features, input_mtraces);
  }

  void FLASHDeconvQuant::calculatePerChargeIsotopeIntensity_(std::vector<double> &per_isotope_intensity,
                                                             std::vector<double> &per_charge_intensity,
                                                             const int max_isotope_count,
                                                             FeatureGroup &fg) const
  {
    int min_pg_charge = INT_MAX;
    int max_pg_charge = INT_MIN;

    for (auto &p: fg)
    {
      if (p.getIsotopeIndex() < 0 || p.getIsotopeIndex() >= max_isotope_count)
      {
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

   void FLASHDeconvQuant::makeMSSpectrum_(std::vector<FeatureSeed *> &local_traces, MSSpectrum &spec, const double &rt) const
  {
    for (auto &tmp_trace : local_traces)
    {
      spec.emplace_back(tmp_trace->getCentroidMz(), tmp_trace->getIntensity());
    }
    spec.setMSLevel(1);
    spec.setName("");
    spec.setRT(rt);
    spec.sortByPosition();
  }

  void FLASHDeconvQuant::getFeatureFromSpectrum_(std::vector<FeatureSeed *> &local_traces,
                                                 std::vector<FeatureGroup> &local_fgroup,
                                                 const double &rt)
  {
     // convert local_traces_ to MSSpectrum
     MSSpectrum spec;
     makeMSSpectrum_(local_traces, spec, rt);

     // run deconvolution
     std::vector<DeconvolutedSpectrum> tmp;
     std::map<int, std::vector<std::vector<double>>> empty;
     DeconvolutedSpectrum deconv_spec = fd_.getDeconvolutedSpectrum(spec, tmp, 0, empty);

     if (deconv_spec.empty()) // if no result was found
     {
       return;
     }

     // convert deconvoluted result into FeatureGroup
     for (auto &deconv : deconv_spec)
     {
       FeatureGroup fg(deconv);
       for (auto &peak: deconv)
       {
         // find seed index
         auto it = std::find_if(local_traces.begin(),
                                local_traces.end(),
                                [=] (FeatureSeed* const& f){
           return (f->getCentroidMz() == peak.mz);
         });
         FeatureSeed tmp_seed(**it);
         tmp_seed.setCharge(peak.abs_charge);
         tmp_seed.setIsotopeIndex(peak.isotopeIndex);
         tmp_seed.setMass(peak.mass);

         fg.push_back(tmp_seed);
       }
       local_fgroup.push_back(fg);
     }
  }

  double FLASHDeconvQuant::getCosine_(const std::vector<double> &a,
                                      const int &a_start,
                                      const int &a_end,
                                      const IsotopeDistribution &b,
                                      const int &b_size,
                                      const int offset)
  {
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
      n += a[j] * b[i].getIntensity();//
    }
    if (a_norm <= 0)
    {
      return 0;
    }
    return n / sqrt(a_norm);
  }

  double FLASHDeconvQuant::scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const
  {
    double diff_mz(std::fabs(tr2.getCentroidMZ() - tr1.getCentroidMZ()));

    double mt_sigma1(tr1.getCentroidSD());
    double mt_sigma2(tr2.getCentroidSD());
    double mt_variances(std::exp(2 * std::log(mt_sigma1)) + std::exp(2 * std::log(mt_sigma2)));

    double mz_score(0.0);
    /// mz scoring by expected mean w/ C13
    double mu = (Constants::C13C12_MASSDIFF_U * iso_pos) / charge; // using '1.0033548378'
    double sd = .0;
    double sigma_mult(3.0);

    //standard deviation including the estimated isotope deviation
    double score_sigma(std::sqrt(std::exp(2 * std::log(sd)) + mt_variances));

    if ((diff_mz < mu + sigma_mult * score_sigma) && (diff_mz > mu - sigma_mult * score_sigma))
    {
      double tmp_exponent((diff_mz - mu) / score_sigma);
      mz_score = std::exp(-0.5 * tmp_exponent * tmp_exponent);
    }

    return mz_score;
  }

  double FLASHDeconvQuant::scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const
  {
    std::map<double, std::vector<double> > coinciding_rts;

    std::pair<Size, Size> tr1_fwhm_idx(tr1.getFWHMborders());
    std::pair<Size, Size> tr2_fwhm_idx(tr2.getFWHMborders());

    double tr1_length(tr1.getFWHM());
    double tr2_length(tr2.getFWHM());
    double max_length = (tr1_length > tr2_length) ? tr1_length : tr2_length;

    // Extract peak shape between FWHM borders for both peaks
    for (Size i = tr1_fwhm_idx.first; i <= tr1_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr1[i].getRT()].push_back(tr1[i].getIntensity());
    }
    for (Size i = tr2_fwhm_idx.first; i <= tr2_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr2[i].getRT()].push_back(tr2[i].getIntensity());
    }

    // Look at peaks at the same RT
    std::vector<double> x, y, overlap_rts;
    for (std::map<double, std::vector<double> >::const_iterator m_it = coinciding_rts.begin(); m_it != coinciding_rts.end(); ++m_it)
    {
      if (m_it->second.size() == 2)
      {
        x.push_back(m_it->second[0]);
        y.push_back(m_it->second[1]);
        overlap_rts.push_back(m_it->first);
      }
    }

    double overlap(0.0);
    if (overlap_rts.size() > 0)
    {
      double start_rt(*(overlap_rts.begin())), end_rt(*(overlap_rts.rbegin()));
      overlap = std::fabs(end_rt - start_rt);
    }

    double proportion(overlap / max_length);
    if (proportion < 0.7)
    {
      return 0.0;
    }
    return computeCosineSim_(x, y);
  }

  double FLASHDeconvQuant::computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const
  {
    if (x.size() != y.size())
    {
      return 0.0;
    }

    double mixed_sum(0.0);
    double x_squared_sum(0.0);
    double y_squared_sum(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
      mixed_sum += x[i] * y[i];
      x_squared_sum += x[i] * x[i];
      y_squared_sum += y[i] * y[i];
    }

    double denom(std::sqrt(x_squared_sum) * std::sqrt(y_squared_sum));
    return (denom > 0.0) ? mixed_sum / denom : 0.0;
  }


  bool FLASHDeconvQuant::doFWHMbordersOverlap(const std::pair<double, double> &border1,
                                              const std::pair<double, double> &border2) const
  {
    if ((border1.first > border2.second) || (border2.first > border1.second))
      return false;

    const double overlap_length = std::min(border1.second, border2.second) - std::max(border1.first, border2.first);
    if ((overlap_length / (border1.second - border1.first) < 0.5) &&
        (overlap_length / (border2.second - border2.first) < 0.5))
      return false;

    return true;
  }

  bool FLASHDeconvQuant::doMassTraceIndicesOverlap(const FeatureGroup &fg1, const FeatureGroup &fg2) const
  {
    // get overlapping charge states
    int min_overlapping_charge = std::max(fg1.getMinCharge(), fg2.getMinCharge());
    int max_overlapping_charge = std::min(fg1.getMaxCharge(), fg2.getMaxCharge());

    if (min_overlapping_charge > max_overlapping_charge) // no overlapping charge
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
    std::set_intersection(mt_indices_1.begin(),
                          mt_indices_1.end(),
                          mt_indices_2.begin(),
                          mt_indices_2.end(),
                          std::back_inserter(inters_vec));

    //    double overlap_p = inters_vec.size() / min_vec_size;
    //    double overlap_percentage = static_cast<double>(inters_vec.size()) / static_cast<double>(min_vec_size);
    // TODO : change this to overlapping only major cs?
    if (inters_vec.size() < 1)
    {
      return false;
    }
    return true;
  }

  bool FLASHDeconvQuant::rescoreFeatureGroup_(FeatureGroup &fg) const
  {
    // update private members in FeatureGroup based on the changed LogMassTraces

    if (!scoreFeatureGroup_(fg))
    {
      return false;
    }

    fg.setFwhmRange();
    fg.setTraceIndices();
    return true;
  }

  bool FLASHDeconvQuant::scoreFeatureGroup_(FeatureGroup &fg) const
  {
    // return false when scoring is not done (filtered out)

    // if this FeatureGroup is within the target, pass any filter
    bool isNotTarget = true;
    fg.setCentroidRtOfApices();
    if (with_target_masses_ && isThisMassOneOfTargets(fg.getMonoisotopicMass(), fg.getRtOfApex()))
    {
      isNotTarget = false;
    }

    auto per_isotope_intensities = std::vector<double>(iso_model_.getMaxIsotopeIndex(), 0);
    auto per_charge_intensities = std::vector<double>(charge_range_, 0);

    // TODO: why max from iso_model_.getMaxIsotopeIndex()???
    calculatePerChargeIsotopeIntensity_(per_isotope_intensities, per_charge_intensities,
                                        iso_model_.getMaxIsotopeIndex(), fg);

    // if the number of charges are not enough
    if (isNotTarget && (fg.getMaxCharge() - fg.getMinCharge() +1 < (int) min_nr_mtraces_))
    {
      return false;
    }

    /// isotope cosine calculation
    int offset = 0;
    double isotope_score =
        FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(fg.getMonoisotopicMass(),
                                                                       per_isotope_intensities,
                                                                       offset,
                                                                       iso_model_,
                                                                       true);
    fg.setIsotopeCosine(isotope_score);
    if ( isNotTarget && isotope_score < min_isotope_cosine_ )
    {
      return false;
    }

    // TODO : remove? does charge distribution make a big difference?
//    double cs = FLASHDeconvAlgorithm::getChargeFitScore_(per_charge_intensities);
//    fg.setChargeScore(cs);
    //    // NOTE : if not using cs dist, too many harmonics
    //    bool is_charge_well_distributed = checkChargeDistribution_(per_charge_intensities);
    //    //double tmp = getChargeFitScore_(per_abs_charge_intensities, charge_range);
    //    if (!is_charge_well_distributed) {
    //        return false;
    //    }
    //    // TODO_ends_here

    /// update monoisotopic mass of this FeatureGroup
    fg.updateMassesAndIntensity(offset, iso_model_.getMaxIsotopeIndex());
    if (isNotTarget && (fg.getMonoisotopicMass() < min_mass_ || fg.getMonoisotopicMass() > max_mass_))
    {
      return false;
    }

    /// setting per_isotope_score
    auto iso_dist = iso_model_.get(fg.getMonoisotopicMass());
    int iso_size = (int) iso_dist.size();

    // initialize private vectors
    fg.initializePerChargeVectors();

    // setting for feature group score
    double feature_score = .0;
//    double weighted_iso_score = fg.getIntensity() / total_intensity_;

    for (int abs_charge = fg.getMinCharge();
         abs_charge <= fg.getMaxCharge();
         ++abs_charge)
    {
      int j = abs_charge - charge_lower_bound_;
      if (per_charge_intensities[j] <= 0)
      {
        continue;
      }

      double max_intensity = .0;
      std::vector<FeatureSeed*> traces_in_this_charge;
      FeatureSeed* apex_trace;
      // find the apex trace in this charge
      for (auto &peak: fg)
      {
        if (peak.getCharge() != abs_charge)
        {
          continue;
        }

        if (peak.getIsotopeIndex() > iso_size)
        {
          continue;
        }

        if (max_intensity < peak.getIntensity())
        {
          max_intensity = peak.getIntensity();
          apex_trace = &peak;
        }
        traces_in_this_charge.push_back(&peak);
      }

      // if no trace is collected
      if (max_intensity == .0)
      {
        continue;
      }

      auto current_per_isotope_intensities = std::vector<double>(iso_model_.getMaxIsotopeIndex(), 0);

      int min_isotope_index = iso_model_.getMaxIsotopeIndex();
      int max_isotope_index = 0;
      double per_charge_score = .0;
      // loop over traces with this charge to collect scores
      for(auto &trace_ptr : traces_in_this_charge)
      {
        current_per_isotope_intensities[trace_ptr->getIsotopeIndex()] += trace_ptr->getIntensity();
        min_isotope_index = min_isotope_index < trace_ptr->getIsotopeIndex() ? min_isotope_index : trace_ptr->getIsotopeIndex();
        max_isotope_index = max_isotope_index < trace_ptr->getIsotopeIndex() ? trace_ptr->getIsotopeIndex() : max_isotope_index;

        // if apex trace, no further scoring
        if (trace_ptr == apex_trace)
        {
          per_charge_score += apex_trace->getIntensity() / fg.getIntensity();
          continue;
        }

        // score per pair between this trace and the apex
        double mz_score(scoreMZ_(*(trace_ptr->getMassTrace()), *(apex_trace->getMassTrace()),
                                 std::abs(trace_ptr->getIsotopeIndex()-apex_trace->getIsotopeIndex()), abs_charge));
        double rt_score(scoreRT_(*(trace_ptr->getMassTrace()), *(apex_trace->getMassTrace())));
        double inty_score(trace_ptr->getIntensity() / fg.getIntensity());
        double total_pair_score = std::exp(std::log(rt_score) + log(mz_score) + log(inty_score));

        per_charge_score += total_pair_score;
      }

      // isotope cosine score for only this charge
      double cos_score = getCosine_(current_per_isotope_intensities,
                                    min_isotope_index,
                                    max_isotope_index,
                                    iso_dist,
                                    iso_size,
                                    0);

//      weighted_iso_score += (per_charge_intensities[j] / total_intensity_) * cos_score;
      feature_score += per_charge_score;
      fg.setChargeIsotopeCosine(abs_charge, cos_score);
      fg.setChargeIntensity(abs_charge, per_charge_intensities[j]);
    }

    fg.setFeatureGroupScore(feature_score);

    return true;
  }

  void FLASHDeconvQuant::refineFeatureGroups_(std::vector<FeatureGroup> &in_features)
  {
    // change min, max charges based on built FeatureGroups (for later use in scoring)
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;

    // output features
    std::vector<FeatureGroup> out_feature;
    out_feature.reserve(in_features.size());

    // sort features by masses
    std::sort(in_features.begin(), in_features.end());

    for (auto &f: in_features)
    {
      min_abs_charge =
          min_abs_charge < f.getMinCharge() ? min_abs_charge : f.getMinCharge();
      max_abs_charge =
          max_abs_charge > f.getMaxCharge() ? max_abs_charge : f.getMaxCharge();
    }
    charge_lower_bound_ = min_abs_charge;
    charge_upper_bound_ = max_abs_charge;
    charge_range_ = charge_upper_bound_ - charge_lower_bound_ + 1;

    Size initial_size = in_features.size();
    this->startProgress(0, initial_size, "refining feature groups");
    // insert FeatureGroup with highest score to out_features, merge if other FeatureGroup exist within mass_tol
    while (!in_features.empty())
    {
      this->setProgress(initial_size - in_features.size());

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

      FeatureGroup lower_fg(candidate_fg->getMonoisotopicMass() - mass_tolerance_da_);
      FeatureGroup upper_fg(candidate_fg->getMonoisotopicMass() + mass_tolerance_da_);

      low_it = std::lower_bound(in_features.begin(), in_features.end(), lower_fg);
      up_it = std::upper_bound(in_features.begin(), in_features.end(), upper_fg);

      // no matching in features (found only candidate itself)
      if (up_it - low_it == 1)
      {
        // save it to out_features
        out_feature.push_back(*candidate_fg);

        // remove candidate from features
        in_features.erase(candidate_fg);
        continue;
      }

      // check if found features are overlapping with the candidate feature
      std::vector<int> v_indices_to_remove;
      v_indices_to_remove.reserve(up_it - low_it);
      std::set<Size> mt_indices_to_add;
      std::vector<FeatureSeed *> mts_to_add;
      mts_to_add.reserve((up_it - low_it) * candidate_fg->size());

      for (; low_it != up_it; ++low_it)
      {
        // if low_it is candidate feature, ignore
        if (candidate_fg == low_it)
        {
          v_indices_to_remove.push_back(low_it - in_features.begin());
          continue;
        }

        // check if fwhm overlaps
        if (!doFWHMbordersOverlap(low_it->getFwhmRange(), candidate_fg->getFwhmRange()))
        {
          continue;
        }

//        // check if masstrace overlaps
//        if (!doMassTraceIndicesOverlap(*low_it, *candidate_fg))
//        {
//          continue;
//        }

        // merge found feature to candidate feature
        auto trace_indices = candidate_fg->getTraceIndices();
        for (auto &new_mt: *low_it)
        {
          // if this mass trace is not used in candidate_fg
          if (std::find(trace_indices.begin(), trace_indices.end(), new_mt.getTraceIndex()) == trace_indices.end())
          {
            mt_indices_to_add.insert(new_mt.getTraceIndex());
            mts_to_add.push_back(&new_mt);
          }
        }
        // add index of found feature to "to_be_removed_vector"
        v_indices_to_remove.push_back(low_it - in_features.begin());
      }

      // sort mts_to_add by abundance
      std::sort(mts_to_add.begin(), mts_to_add.end(), CmpFeatureSeedByIntensity());

      // add extra masstraces to candidate_feature
      FeatureGroup final_candidate_fg = *candidate_fg; // copy of candidate_feature
      for (auto &new_mt: mts_to_add)
      {
        // check if this mt of this lmt was already added to candidate_fg
        if (mt_indices_to_add.find(new_mt->getTraceIndex()) == mt_indices_to_add.end())
        {
          continue;
        }
        mt_indices_to_add.erase(new_mt->getTraceIndex());

        FeatureSeed *apex_lmt_in_this_cs = final_candidate_fg.getApexLMTofCharge(new_mt->getCharge());
        // if this mt is introducing new charge
        if (apex_lmt_in_this_cs == nullptr)
        {
          final_candidate_fg.push_back(*new_mt);
          continue;
        }

        /// re-calculate isotope index (from FLASHDeconvAlgorithm::getCandidatePeakGroup
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
          // if new iso_idx is not within iso_range (too small)
          continue;
        }
        if (apex_lmt_in_this_cs->getIsotopeIndex() + tmp_iso_idx > iso_model_.getMaxIsotopeIndex())
        {
          // if new iso_idx is not within iso_range (too large)
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
      tmp_out_fgs.reserve(in_features.size() - v_indices_to_remove.size());
      for (Size i = 0; i < in_features.size(); ++i)
      {
        if (std::find(v_indices_to_remove.begin(), v_indices_to_remove.end(), i) == v_indices_to_remove.end())
        {
          tmp_out_fgs.push_back(in_features[i]);
        }
      }
      in_features.swap(tmp_out_fgs);
    }
    this->endProgress();

    in_features.swap(out_feature);
  }

  void FLASHDeconvQuant::buildMassTraceGroups_(std::vector<FeatureSeed> &mtraces,
                                               std::vector<FeatureGroup> &features)
  {
    /// FLASHDeconvAlgorithm setting
    Param fd_defaults = getFLASHDeconvParams_();
    fd_.setParameters(fd_defaults);
    fd_.calculateAveragine(false);
    iso_model_ = fd_.getAveragine();
//    std::vector<double> target_masses_; // monoisotope
//    fd_.setTargetMasses(target_masses_, ms_level);

    /// group mass traces to spectrum
    std::vector<std::pair<double, FeatureSeed*>> mt_rt_starts;
    std::vector<std::pair<double, FeatureSeed*>> mt_rt_ends;
    mt_rt_starts.reserve(mtraces.size());
    mt_rt_ends.reserve(mtraces.size());
    int counter = 0;

    // collect rt information from mtraces to generate spectrum
    double min_fwhm_length = numeric_limits<double>::max();
    for (auto &trace: mtraces)
    {
      mt_rt_starts.push_back(std::make_pair(trace.getFwhmStart(), &trace));
      mt_rt_ends.push_back(std::make_pair(trace.getFwhmEnd(), &trace));
      if (trace.getMassTrace()->getFWHM() < min_fwhm_length)
      {
        min_fwhm_length = trace.getMassTrace()->getFWHM();
      }
    }

    if (min_fwhm_length > rt_window_)
    {
      rt_window_ = min_fwhm_length;
    }
//    OPENMS_LOG_INFO << "spectrum rt window : " << std::to_string(rt_window_) << endl;

    // sorting mass traces in rt
    std::sort(mt_rt_starts.begin(), mt_rt_starts.end());
    std::sort(mt_rt_ends.begin(), mt_rt_ends.end());

    std::vector<std::pair<double, FeatureSeed *>>::const_iterator rt_s_iter = mt_rt_starts.begin();
    std::vector<std::pair<double, FeatureSeed *>>::const_iterator rt_e_iter = mt_rt_ends.begin();
    auto end_of_iter = mt_rt_starts.end();
    double end_of_current_rt_window = mt_rt_starts[0].first;
    double last_rt = mt_rt_ends[mt_rt_ends.size() - 1].first;

    // mass traces to be added in a spectrum
    std::vector<FeatureSeed *> local_traces;
    local_traces.reserve(mtraces.size());

    int possible_spec_size = int((mt_rt_starts[mt_rt_starts.size() - 1].first - end_of_current_rt_window) / rt_window_);
    this->startProgress(0, possible_spec_size, "assembling mass traces to features");

    while (rt_s_iter != end_of_iter && end_of_current_rt_window < last_rt)
    {
      this->setProgress(counter);

      // initial rt binning is 1 sec (for generating spectrum)
      end_of_current_rt_window += rt_window_;

      // add mass traces within rt range
      bool is_new_mt_added = false;
      for (; rt_s_iter != end_of_iter && rt_s_iter->first <= end_of_current_rt_window; ++rt_s_iter)
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
      for (; rt_e_iter != mt_rt_ends.end() && rt_e_iter->first < end_of_current_rt_window; ++rt_e_iter)
      {
        local_traces.erase(std::remove_if(local_traces.begin(), local_traces.end(),
                                          [&rt_e_iter](auto const &p) { return rt_e_iter->second == p; }));
      }
      if (local_traces.empty())
      {
        continue;
      }

      // sort local traces in mz
      sort(local_traces.begin(), local_traces.end(), CmpFeatureSeedByMZ());

      std::vector<FeatureGroup> local_fgroup;
      getFeatureFromSpectrum_(local_traces, local_fgroup, end_of_current_rt_window);
      ++counter; // to track the number of generated spectra
      if (local_fgroup.size() == 0)
      {
        continue;
      }

      for (auto &pg: local_fgroup)
      {
        sort(pg.begin(), pg.end());
        pg.setFwhmRange();
        pg.setTraceIndices();
        pg.setChargeVector();
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

    for (auto &mt: fg)
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
    String s_out = outfile_path.substr(0, outfile_path.find_last_of(".") - 1) + "_sharedMTs.tsv";
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
      for (auto &f_idx: shared_m_traces[mt_idx])
      {
        const auto &feature = fgroups[f_idx];

        // find charge of shared mass traces
        int feat_cs;
        Size lmt_idx_in_feat;
        for (Size lmt_idx = 0; lmt_idx < feature.size(); ++lmt_idx)
        {
          auto &lmt = feature[lmt_idx];
          if (mt_idx != lmt.getTraceIndex())
            continue;
          feat_cs = lmt.getCharge();
          lmt_idx_in_feat = lmt_idx;
          break;
        }

        // find max quant of feature where shared mass trace belongs (shared cannot be the max)
        double max_quant = .0;
        for (Size lmt_idx = 0; lmt_idx < feature.size(); ++lmt_idx)
        {
          auto &lmt = feature[lmt_idx];
          if (mt_idx == lmt.getTraceIndex())
            continue; // if this lmt is the shared one
          if (feat_cs != lmt.getCharge())
            continue;

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
                                            [](auto &a, auto &b) { return a + b.first; });

      for (auto &q_pair: feat_quant_pair)
      {
        double weight = q_pair.first / sum_of_quant;
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
    for (auto &f: fgroups)
    {
      f.updateIntensity();
    }

    /// test writing for all mass traces in FeatureGroup
    String out_path = outfile_path.substr(0, outfile_path.find_last_of(".") - 1) + "groups.tsv";
    std::fstream out;
    out.open(out_path, fstream::out);
    out << "mono_mass\tcharge\tisotope_index\tquant_value\tshared\tcentroid_mzs\trts\tmzs\tintys\n"; // header
    writeMassTracesOfFeatureGroup(fgroups, shared_m_traces, out);
    out.close();
  }


  void FLASHDeconvQuant::resolveSharedMassTraces(std::vector<FeatureGroup> &fgroups,
                                                 std::vector<std::vector<Size>> &shared_m_traces,
                                                 std::vector<MassTrace> &input_mtraces) const
  {
    /// generate array with indices of shared mass traces
    for (Size fg_index = 0; fg_index < fgroups.size(); ++fg_index)
    {
      for (auto &mt_i: fgroups[fg_index].getTraceIndices())
      {
        shared_m_traces[mt_i].push_back(fg_index);
      }
    }

    /// test writing for shared mass traces
    String s_out = outfile_path.substr(0, outfile_path.find_last_of(".") - 1) + "_sharedMTs.tsv";
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
      for (auto &f_idx: shared_m_traces[mt_idx])
      {
        const auto &feature = fgroups[f_idx];
        if (feature.getIsotopeCosine() == 0) // if this featureGroup is rejected in other competition
        {
          continue;
        }
        for (Size lmt_idx = 0; lmt_idx < feature.size(); ++lmt_idx)
        {
          auto &lmt = feature[lmt_idx];
          if (mt_idx != lmt.getTraceIndex())
            continue;

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
      for (auto &f_idx: shared_m_traces[mt_idx])
      {
        if (f_idx != feat_idx_with_max_score)
        {
          auto &feature = fgroups[f_idx];
          for (auto lmt_iter = feature.begin(); lmt_iter != feature.end(); ++lmt_iter)
          {
            if (mt_idx != lmt_iter->getTraceIndex())
              continue;
            feature.erase(lmt_iter);
            break;
          }
        }
        if (!rescoreFeatureGroup_(fgroups[f_idx]))
        {
          fgroups[f_idx].setIsotopeCosine(0.0);
        }
      }
    }

    // rescore
    std::vector<FeatureGroup> tmpFGroups;
    tmpFGroups.swap(fgroups);
    fgroups.reserve(tmpFGroups.size());
    for (auto &f: tmpFGroups)
    {
      if (!rescoreFeatureGroup_(f))
      {
        continue;
      }
      fgroups.push_back(f);
    }

    /// test writing for all mass traces in FeatureGroup
    String out_path = outfile_path.substr(0, outfile_path.find_last_of(".") - 1) + "groups.tsv";
    std::fstream out;
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
  void FLASHDeconvQuant::clusterFeatureGroups_(std::vector<FeatureGroup> &fgroups,
                                               std::vector<MassTrace> &input_mtraces) const
  {
    // *********************************************************** //
    // Step 1 preparation for hypergraph : collect feature idx with shared mass traces
    // *********************************************************** //
    std::vector<std::vector<Size>> shared_m_traces(input_mtraces.size(), std::vector<Size>());
    for (Size fg_index = 0; fg_index < fgroups.size(); ++fg_index)
    {
      for (auto &mt_i: fgroups[fg_index].getTraceIndices())
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
    //    std::fstream out;
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
          if (shared_m_traces[i].size() == 0)
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

          for (const auto &mt_index: current_fg.getTraceIndices())
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
      if (fg_indices_in_current_cluster.size() == 1)
      {
        out_features.push_back(fgroups[*(fg_indices_in_current_cluster.begin())]);
        //        out_feature_idx.push_back(*(fg_indices_in_current_cluster.begin()));
        continue;
      }

      // resolve the conflict among feature groups
      cluster_counter++;
      String cluster_name = "cluster" + to_string(cluster_counter);
      resolveConflictInCluster_(fgroups, input_mtraces, shared_m_traces, fg_indices_in_current_cluster, out_features);
    }
    //    this->endProgress();

    out_features.shrink_to_fit();
    std::swap(out_features, fgroups);
    //    out_feature_idx.shrink_to_fit();
    //    out.close();
    OPENMS_LOG_INFO << "#final feature groups: " << fgroups.size() << endl;
    OPENMS_LOG_INFO << "#cluster (shared) :" << cluster_counter << endl;
  }

  void FLASHDeconvQuant::writeMassTracesOfFeatureGroup(const std::vector<FeatureGroup> &feature_groups,
                                                       const std::vector<std::vector<Size>> &shared_m_traces_indices,
                                                       std::fstream &out) const
  {
    OPENMS_LOG_INFO << "----- # FeatureGroup : " << feature_groups.size() << " -----" << endl;
    for (const auto &feat: feature_groups)
    {
      double quant(0.0);
      auto mt_idxs = feat.getTraceIndices();

      for (const auto &mt: feat)
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

        for (const auto &peak: *(mt.getMassTrace()))
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

  void FLASHDeconvQuant::resolveConflictInCluster_(std::vector<FeatureGroup> &feature_groups,
                                                   const std::vector<MassTrace> &input_masstraces,
                                                   const std::vector<std::vector<Size> > &shared_m_traces_indices,
                                                   const std::set<Size> &fg_indices_in_this_cluster,
                                                   std::vector<FeatureGroup> &out_featuregroups) const
  {
    /// conflict resolution is done in feature level (not feature group level)

    // prepare nodes for hypergraph
    std::vector<std::vector<Size>> mt_and_feature_idx(shared_m_traces_indices.size(), std::vector<Size>());
    //                                                                      std::vector<std::pair<Size, Size>>()); // pair <fg_idx_in_cluster, charge_idx>
    std::vector<FeatureElement> feature_candidates; // index : (fg_indices_in_this_cluster+1)*charge_idx
    for (auto &fg_idx: fg_indices_in_this_cluster)
    {
      auto &fg = feature_groups[fg_idx];
      OPENMS_LOG_DEBUG << fg_idx << "\t" << fg.getMonoisotopicMass() << endl;

      // get charge vector
      Size feature_idx_starts = feature_candidates.size();
      fg.setChargeVector();
      std::vector<int> current_fg_cs = fg.getChargeVector();
      for (auto &tmp_cs : current_fg_cs)
      {
        FeatureElement tmp;
        tmp.mass_traces.reserve(fg.size());
        tmp.mass_trace_indices.reserve(fg.size());
        tmp.charge = tmp_cs;
        tmp.feature_group_index = fg_idx;
        feature_candidates.push_back(tmp);
      }

      for (auto lmt_iter = fg.begin(); lmt_iter != fg.end(); ++lmt_iter)
      {
        Size f_idx = feature_idx_starts + std::distance(current_fg_cs.begin(),
                                                        std::find(current_fg_cs.begin(), current_fg_cs.end(), lmt_iter->getCharge()));
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

      while (!bfs_queue.empty())
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

          for (const auto &mt_index: selected_feature.mass_trace_indices)
          {
            if (!bfs_visited[mt_index])
            {
              bfs_queue.push(mt_index);
              bfs_visited[mt_index] = true;
            }
            else if (mt_and_feature_idx[mt_index].size() > 1)
            {
              // check if this visited mass trace is the shared one (because the seed mass trace can end up here as well)
              conflicting_mt_indices.insert(mt_index);
            }
          }
        }
      }

      // this feature is not sharing any mass traces with others
      if (feature_idx_in_current_conflict_region.size() == 1)
      {
        continue;
      }

      // collect conflicting mass traces (not LogMassTrace, originals)
      std::vector<const MassTrace *> conflicting_mts;
      for (auto &mt_idx: conflicting_mt_indices)
      {
        auto i = input_masstraces.begin() + mt_idx;
        conflicting_mts.push_back(&(*i));
      }

      // set isotope probabilities for feature_candidates && check if mass traces are only shared ones
      vector<Size> features_not_for_resolution;
      for (auto &feat_idx: feature_idx_in_current_conflict_region)
      {
        auto &feat = feature_candidates[feat_idx];
        auto &feat_group = feature_groups[feat.feature_group_index];
        auto tmp_iso = iso_model_.get(feat_group.getMonoisotopicMass());
        feat.isotope_probabilities.reserve(feat.mass_traces.size());
        for (auto &lmt: feat.mass_traces)
        {
          // if isotope index of lmt exceed tmp_iso length, give 0
          int tmp_iso_idx = lmt->getIsotopeIndex();
          if (tmp_iso_idx < (int) tmp_iso.size())
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
        for (auto &mt_idx: feat.mass_trace_indices)
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
        FeatureSeed *most_abundant_mt = nullptr;
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
      if (features_not_for_resolution.size() > 0)
      {
        for (Size &idx: features_not_for_resolution)
        {
          feature_idx_in_current_conflict_region.erase(idx);
          for (auto &lmt: feature_candidates[idx].mass_traces)
          {
            lmt->setIntensity(0.0);
          }
        }

        if (feature_idx_in_current_conflict_region.size() < 2)
        {
          continue;
        }

        // check if any shared mass trace is not shared anymore
        vector<Size> is_mt_included(conflicting_mts.size(),
                                    0); // element is equal to the size of feature_vec when it's shared
        for (auto &feat_idx: feature_idx_in_current_conflict_region)
        {
          auto &lmts = feature_candidates[feat_idx].mass_traces;
          for (auto &lmt: lmts)
          {
            for (Size mt_i = 0; mt_i < conflicting_mts.size(); ++mt_i)
            {
              if (conflicting_mts[mt_i] == lmt->getMassTrace())
              {
                is_mt_included[mt_i] += 1;
                break;
              }
            }
          }
        }
        for (int mt_i = conflicting_mts.size() - 1; mt_i >= 0; --mt_i)
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
      for (auto &i: feature_idx_in_current_conflict_region)
        feature_indices_vec.push_back(i);

      // resolve this conflict region
      resolveConflictRegion_(feature_candidates, feature_indices_vec, conflicting_mts);

      OPENMS_LOG_DEBUG << "---------conflicting masses----------" << endl;
      for (auto &i: feature_indices_vec)
      {
        auto fg_index = feature_candidates[i].feature_group_index;
        OPENMS_LOG_DEBUG << to_string(feature_groups[fg_index].getMonoisotopicMass()) << endl;
      }
    }

    // update feature group quantities
    for (auto &idx: fg_indices_in_this_cluster)
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
        else
        {
          ++fg_iter;
        }
      }

      // check if enough number of mass traces are left
      if (fgroup.size() < min_nr_mtraces_)
      {
        continue;
      }

      // check if the FeatureGroup still contains more than three charges
      fgroup.setChargeVector();
      if (fgroup.getChargeVector().size() < 3)
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
                                                const std::vector<const MassTrace *> &conflicting_mts) const
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
      for (Size idx = 0; idx < mass_traces_size; ++idx)
      {
        // ignore this mt if it's the conflicting one
        auto tmp_iter = std::find(conflicting_mts.begin(),
                                  conflicting_mts.end(),
                                  tmp_feat.mass_traces[idx]->getMassTrace());
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
        for (auto &p_2d: (*mass_trace_ptr))
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
          auto &tmp_trace = traces_for_fitting[0];
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
      EGHTraceFitter *fitter = new EGHTraceFitter();
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
        for (auto &peak: *(conflicting_mts[row]))
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
      Size column_size =
          components_indices.size() - std::count(components_indices.begin(), components_indices.end(), -1);

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
        if (components_indices[comp_idx] == -1)
          continue;

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

      OPENMS_LOG_DEBUG << "-------------------------" << endl;
      OPENMS_LOG_DEBUG << "observed m/z : " << std::to_string(conflicting_mts[row]->getCentroidMZ()) << endl;
      for (auto &peak: (*conflicting_mts[row]))
      {
        OPENMS_LOG_DEBUG << std::to_string(peak.getRT()) << "\t" << std::to_string(peak.getIntensity()) << "\n";
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

        OPENMS_LOG_DEBUG << "--- theo[" << col << "]--- " << std::to_string(out_quant.getValue(col, 0)) << "\t"
                         << std::to_string(lmt_ptr->getIntensity()) << "\n";
        for (auto &m: theo_matrix.col(col))
        {
          OPENMS_LOG_DEBUG << to_string(m) << "\t";
        }
        OPENMS_LOG_DEBUG << "\n";
        col++;
      }
    }
    components.clear();
  }

  // from ElutionModelFitter. TODO : this is test purpose. use the existing method insteaad.
  void FLASHDeconvQuant::runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces,
                                             EGHTraceFitter *fitter) const
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
    catch (Exception::UnableToFit &except)
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
                                                                  FeatureSeed *&most_abundant_mt_ptr,
                                                                  const std::vector<std::vector<Size>> &shared_m_traces) const
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
        most_abundant_mt_ptr = (FeatureSeed *) &(*lmt);
      }
    }
  }

  void FLASHDeconvQuant::getFLASHDeconvConsensusResult()
  {
    std::fstream fin;
//    fin.open("/Users/jeek/Documents/A4B/FDQ/Kiel-Human/FD_consensus.csv", std::ios::in);
    fin.open("/Users/jeek/Documents/A4B/FDQ/Kiel-Human/ProSightPD_concensus.csv", std::ios::in);

    std::string line;
    std::getline(fin, line); // skip the first line
    while (std::getline(fin, line))
    {
      Size comma_loc = line.find(',');
      target_masses_.push_back(std::make_pair(stod(line.substr(0, comma_loc)), stod(line.substr(comma_loc + 1))));
//      target_masses_.push_back(stod(line.substr(0, comma_loc)));
//      rt_of_target_masses_.push_back(stod(line.substr(comma_loc + 1)));
    }

    std::sort(target_masses_.begin(), target_masses_.end());
    fin.close();
  }

  bool FLASHDeconvQuant::isThisMassOneOfTargets(const double &candi_mass, const double &candi_rt) const
  {
    auto low_it = std::lower_bound(target_masses_.begin(), target_masses_.end(), std::make_pair(candi_mass - 1.5, .0));
    auto up_it = std::upper_bound(target_masses_.begin(), target_masses_.end(), std::make_pair(candi_mass + 1.5, .0));

    bool is_it_in_the_list = false;
    // if only one mass is found
    if (low_it == up_it)
    {
      // check if any mass is within range
      if (abs(low_it->first-candi_mass) <= 1.5 && abs(low_it->second-candi_rt) < 180)
      {
        return true;
      }
    }
    for (auto tmp_it = low_it; tmp_it != up_it; ++tmp_it)
    {
      // check if any mass is within range
      if (abs(tmp_it->first-candi_mass) <= 1.5 && abs(tmp_it->second-candi_rt) < 180)
      {
        is_it_in_the_list = true;
      }
    }
    return is_it_in_the_list;
  }
}