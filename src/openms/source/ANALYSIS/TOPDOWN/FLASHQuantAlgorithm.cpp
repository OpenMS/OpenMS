// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <include/OpenMS/ANALYSIS/TOPDOWN/FLASHQuantAlgorithm.h>
#include <include/OpenMS/CONCEPT/LogStream.h>
#include <include/OpenMS/DATASTRUCTURES/Matrix.h>
#include <include/OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>
#include <include/OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <queue>

namespace OpenMS
{
  FLASHQuantAlgorithm::FLASHQuantAlgorithm() :
      ProgressLogger(),
      DefaultParamHandler("FLASHQuantAlgorithm")
  {
    defaults_.setValue("min_charge", 5, "Minimum charge state to consider");
    defaults_.setValue("max_charge", 50, "Maximum charge state to consider");
    defaults_.setValue("min_mass", 10000, "Minimum mass (Da)");
    defaults_.setValue("max_mass", 70000, "Maximum mass (Da)");
    defaults_.setValue("mz_tol", 10, "Ppm tolerance for m/z values in deconvolution");
    defaults_.setValue("mass_tol", 3, "Mass tolerance in Dalton for integrating similar feature groups into a single one");
    defaults_.setValue("min_isotope_cosine", 0.85, "Cosine threshold between averagine and observed isotope pattern. Note that 0.8 is used for deconvolution", {"advanced"});

    defaults_.setValue("use_smoothed_intensities",
                       "true",
                       "Use LOWESS intensities instead of raw intensities.",
                       {"advanced"});
    defaults_.setValidStrings("use_smoothed_intensities", {"false", "true"});
    defaults_.setValue("out_shared_details", "false", "Outputs a tsv file including detailed information about the resolved signals (filename = <out_file_name>_shared.tsv", {"advanced"});
    defaults_.setValidStrings("out_shared_details", {"false", "true"});
    defaults_.setValue("resolving_shared_signal", "true", "Resolve shared signals between feature groups (i.e., co-elution)", {"advanced"});
    defaults_.setValidStrings("resolving_shared_signal", {"true", "false"});
    defaultsToParam_();

    this->setLogType(CMD);
  }

  void FLASHQuantAlgorithm::updateMembers_()
  {
    charge_lower_bound_ = (Size) param_.getValue("min_charge");
    charge_upper_bound_ = (Size) param_.getValue("max_charge");

    min_mass_ = (double) param_.getValue("min_mass");
    max_mass_ = (double) param_.getValue("max_mass");

    mz_tolerance_ = (double) param_.getValue("mz_tol");
    mass_tolerance_da_ = (double) param_.getValue("mass_tol");
    min_isotope_cosine_ = (double) param_.getValue("min_isotope_cosine");

    use_smoothed_intensities_ = param_.getValue("use_smoothed_intensities").toBool();
    shared_output_requested_ = param_.getValue("out_shared_details").toBool();
    resolving_shared_signal_ = param_.getValue("resolving_shared_signal").toBool();
  }

  Param FLASHQuantAlgorithm::getFLASHDeconvParams_() const
  {
    Param fd_defaults = SpectralDeconvolution().getDefaults();
    // overwrite algorithm default so we export everything (important for copying back MSstats results)
    fd_defaults.setValue("min_charge", (int) charge_lower_bound_);
    fd_defaults.setValue("max_charge", (int) charge_upper_bound_);
    fd_defaults.setValue("min_mass", min_mass_);
    fd_defaults.setValue("max_mass", max_mass_);
    fd_defaults.setValue("min_cos", DoubleList{.8, .8});
    fd_defaults.setValue("tol", DoubleList{mz_tolerance_, 10.0});
    return fd_defaults;
  }

  void FLASHQuantAlgorithm::run(std::vector<MassTrace> &input_mtraces, std::vector<FeatureGroup>& out_fgs)
  {
    // *********************************************************** //
    // Step 1 deconvolute mass traces
    // *********************************************************** //
    // initialize input & output
    std::vector<FeatureSeed> input_seeds;
    input_seeds.reserve(input_mtraces.size());
    std::vector<MassTrace> updated_masstraces;
    updated_masstraces.reserve(input_mtraces.size());
    Size index = 0;
    for (auto iter = input_mtraces.begin(); iter != input_mtraces.end(); ++iter) {
      if (iter->getSize() < min_nr_peaks_in_mtraces_)
      {
        continue;
      }
      FeatureSeed tmp_seed(*iter);
      tmp_seed.setTraceIndex(index);
      input_seeds.push_back(tmp_seed);
      updated_masstraces.push_back(*iter);
      ++index;
    }
    input_mtraces.swap(updated_masstraces);
    updated_masstraces.clear();

    // sort input mass traces in RT
    std::sort(input_seeds.begin(), input_seeds.end(), FLASHQuantHelper::CmpFeatureSeedByRT());
    std::vector<FeatureGroup> features;
    features.reserve(input_mtraces.size());

    // run deconvolution
    buildMassTraceGroups_(input_seeds, features);
    features.shrink_to_fit();

    // *********************************************************** //
    // Step 2 mass artifact removal & post processing...
    // *********************************************************** //
    refineFeatureGroups_(features);
    OPENMS_LOG_INFO << "#Detected feature groups : " << features.size() << endl;

    // *********************************************************** //
    // Step 3 clustering features
    // *********************************************************** //
    if (resolving_shared_signal_)
    {
      if (shared_output_requested_)
      {
        setOptionalDetailedOutput_();
      }

      clusterFeatureGroups_(features, input_mtraces);
      if (shared_output_requested_)
      {
        shared_out_stream_.close();
      }
    }
    OPENMS_LOG_INFO << "#Final feature groups: " << features.size() << endl;

    // output
    setFeatureGroupMembersForResultWriting_(features);
    out_fgs.swap(features);
    features.clear();
  }

   void FLASHQuantAlgorithm::makeMSSpectrum_(std::vector<FeatureSeed *> &local_traces, MSSpectrum &spec, const double &rt) const
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

  void FLASHQuantAlgorithm::setFeatureGroupMembersForResultWriting_(std::vector<FeatureGroup> &f_groups) const
  {
    // this cannot be done in FeatureGroup (in FLASHQuantHelper) due to some methods to be used only in here

    for (auto &fgroup : f_groups)
    {
      // initialize
      auto per_charge_int = std::vector<float>(1 + fgroup.getMaxCharge(), .0);
      auto per_charge_cos = std::vector<float>(1 + fgroup.getMaxCharge(), .0);

      std::vector<std::vector<float>> per_cs_isos(1 + fgroup.getMaxCharge(), std::vector<float>(fgroup.getIsotopeIntensities().size(), .0));

      // iterate all FeatureSeeds to collect per feature values
      for(auto& seed : fgroup)
      {
        if (seed.getIsotopeIndex() < 0)
        {
          continue;
        }
        per_charge_int[seed.getCharge()] += seed.getIntensity();
        per_cs_isos[seed.getCharge()][seed.getIsotopeIndex()] += seed.getIntensity();
      }

      for(auto& cs : fgroup.getChargeSet())
      {
        auto &this_cs_isos = per_cs_isos[cs];
        int min_isotope_index = std::distance( this_cs_isos.begin(), std::find_if( this_cs_isos.begin(), this_cs_isos.end(), [](auto &x) { return x != 0; }));
        int max_isotope_index = std::distance( this_cs_isos.begin(), (std::find_if( this_cs_isos.rbegin(), this_cs_isos.rend(), [](auto &x) { return x != 0; }) + 1).base());

        auto iso_dist = iso_model_.get(fgroup.getMonoisotopicMass());
        float cos_score = SpectralDeconvolution::getCosine(per_cs_isos[cs], min_isotope_index, max_isotope_index + 1,
                                                           iso_dist, iso_dist.size(), 0, 2); // min_iso_size_: based on FLASHDeconv
        per_charge_cos[cs] = cos_score;
      }
      // calculate average mass
      auto avg_mass = iso_model_.getAverageMassDelta(fgroup.getMonoisotopicMass()) + fgroup.getMonoisotopicMass();

      // setting values
      fgroup.setPerChargeIntensities(per_charge_int);
      fgroup.setPerChargeCosineScore(per_charge_cos);
      fgroup.setAverageMass(avg_mass);
    }
  }

  void FLASHQuantAlgorithm::getFeatureFromSpectrum_(std::vector<FeatureSeed *> &local_traces,
                                                          std::vector<FeatureGroup> &local_fgroup,
                                                          const double &rt)
  {
     // convert local_traces_ to MSSpectrum
     MSSpectrum spec;
     makeMSSpectrum_(local_traces, spec, rt);

     // run deconvolution
     PeakGroup empty_pg; // empty one, since only MS1s are considered.
     deconv_.performSpectrumDeconvolution(spec, 0, empty_pg);
     DeconvolvedSpectrum& deconv_spec = deconv_.getDeconvolvedSpectrum();

     if (deconv_spec.empty()) // if no result was found
     {
       return;
     }

     // convert deconvolved result into FeatureGroup
     for (auto &deconv : deconv_spec)
     {
       // filter out if deconv results are not sufficient
       if (deconv.size() < min_nr_mtraces_)
       {
         continue;
       }

       FeatureGroup fg(deconv);
       fg.setMaxIsotopeIndex(iso_model_.get(deconv.getMonoMass()).size());

       // Add individual FeatureSeeds to FeatureGroup
       for (auto &peak: deconv)
       {
         // if isotope index of this peak is out of threshold, don't include this
         if (peak.isotopeIndex >= (int) fg.getMaxIsotopeIndex())
         {
           continue;
         }

         // find seed index
         auto it = std::find_if(local_traces.begin(), local_traces.end(),
                                [peak] (FeatureSeed* const& f) -> bool {
                                  return (f->getCentroidMz() == peak.mz); // intensity changes in performSpectrumDeconvolution, thus cannot be used for filtering
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

  bool FLASHQuantAlgorithm::doFWHMbordersOverlap_(const std::pair<double, double> &border1,
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

  bool FLASHQuantAlgorithm::doMassTraceIndicesOverlap(const FeatureGroup &fg1, const FeatureGroup &fg2,
                                                            const double overlap_percentage_threshold, const bool charge_specific) const
  {
    // get overlapping charge states
    int min_overlapping_charge = std::max(fg1.getMinCharge(), fg2.getMinCharge());
    int max_overlapping_charge = std::min(fg1.getMaxCharge(), fg2.getMaxCharge());

    if (charge_specific && (min_overlapping_charge > max_overlapping_charge)) // no overlapping charge
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
      if (!charge_specific || (fg1[fg1_idx].getCharge() >= min_overlapping_charge && fg1[fg1_idx].getCharge() <= max_overlapping_charge))
      {
        mt_indices_1.push_back(fg1[fg1_idx].getTraceIndex());
      }
    }
    if (mt_indices_1.empty())
    {
      return false;
    }
    for (Size fg2_idx = 0; fg2_idx < fg2.size(); ++fg2_idx)
    {
      if (!charge_specific || (fg2[fg2_idx].getCharge() >= min_overlapping_charge && fg2[fg2_idx].getCharge() <= max_overlapping_charge))
      {
        mt_indices_2.push_back(fg2[fg2_idx].getTraceIndex());
      }
    }
    if (mt_indices_2.empty())
    {
      return false;
    }

    std::sort(mt_indices_1.begin(), mt_indices_1.end());
    std::sort(mt_indices_2.begin(), mt_indices_2.end());

    Size min_vec_size = std::min(mt_indices_1.size(), mt_indices_2.size());
    std::vector<Size> inters_vec;
    inters_vec.reserve(min_vec_size);
    std::set_intersection(mt_indices_1.begin(),
                          mt_indices_1.end(),
                          mt_indices_2.begin(),
                          mt_indices_2.end(),
                          std::back_inserter(inters_vec));

    double overlap_percentage = static_cast<double>(inters_vec.size()) / static_cast<double>(min_vec_size);
    // TODO : change this to overlapping only major cs?
    if (overlap_percentage < overlap_percentage_threshold)
    {
      return false;
    }
    return true;
  }

  bool FLASHQuantAlgorithm::rescoreFeatureGroup_(FeatureGroup &fg) const
  {
    if ((!scoreAndFilterFeatureGroup_(fg)))
    {
      return false;
    }
    // update private members in FeatureGroup based on the changed FeatureSeeds
    fg.updateMembers();
    return true;
  }

  bool FLASHQuantAlgorithm::scoreAndFilterFeatureGroup_(FeatureGroup &fg, double min_iso_score) const
  {
    /// return false when scoring is not done (filtered out)
    // update monoisotopic mass, isotope_intensities_ and charge vector to use for filtering
    fg.updateMembersForScoring();

    /// if this FeatureGroup is within the target, pass any filter
    bool isNotTarget = true;
//    if (with_target_masses_ && isThisMassOneOfTargets_(fg.getMonoisotopicMass(), fg.getRtOfMostAbundantMT()))
//    {
//      isNotTarget = false;
//    }

    /// filter if the number of charges are not enough
    if (isNotTarget && (fg.getChargeSet().size() < min_nr_mtraces_))
    {
      return false;
    }

    /// filter if the mass is out of range
    if (isNotTarget && (fg.getMonoisotopicMass() < min_mass_ || fg.getMonoisotopicMass() > max_mass_))
    {
      return false;
    }

    /// isotope cosine calculation
    /// based on: SpectralDeconvolution::scoreAndFilterPeakGroups_()
    if (min_iso_score == -1)
    {
      min_iso_score = min_isotope_cosine_;
    }
    int offset = 0;
    float isotope_score = SpectralDeconvolution::getIsotopeCosineAndIsoOffset(fg.getMonoisotopicMass(),
                                                                              fg.getIsotopeIntensities(),
                                                                              offset, iso_model_, 1);
    fg.setIsotopeCosine(isotope_score);
    if (isNotTarget && isotope_score < min_iso_score)
    {
      return false;
    }
    // update values based on the calculated offset
    if (offset > 0)
    {
      fg.setMonoisotopicMass(fg.getMonoisotopicMass() + iso_da_distance_ * offset);
      fg.updateIsotopeIndices(offset); // TODO: change?
    }

    return true;
  }

  void FLASHQuantAlgorithm::refineFeatureGroups_(std::vector<FeatureGroup> &in_features)
  {
    // change min, max charges based on built FeatureGroups (for later use in scoring)
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;

    // minimum isotope cosine score for the refining FeatureGroups
    double min_iso_score = 0.5;
    if (!resolving_shared_signal_) // if not resolving shared signal, use the final min_isotope_cosine_ (because this is the last step)
    {
      min_iso_score = min_isotope_cosine_;
    }

    // output features
    std::vector<FeatureGroup> out_feature;
    out_feature.reserve(in_features.size());

    // sort features by masses
    std::sort(in_features.begin(), in_features.end());

    // setting pointer vector (to reduce time complexity for push_back & erase)
    std::vector<FeatureGroup*> in_feature_ptrs;
    in_feature_ptrs.reserve(in_features.size());

    // set variables according to the detected features
    for (auto &f: in_features)
    {
      in_feature_ptrs.push_back(&f);

      min_abs_charge =
          min_abs_charge < f.getMinCharge() ? min_abs_charge : f.getMinCharge();
      max_abs_charge =
          max_abs_charge > f.getMaxCharge() ? max_abs_charge : f.getMaxCharge();
    }
    charge_lower_bound_ = min_abs_charge;
    charge_upper_bound_ = max_abs_charge;

    Size initial_size = in_features.size();
    this->startProgress(0, initial_size, "refining feature groups");
    // insert FeatureGroup with the highest score to out_features, merge if other FeatureGroup exist within mass_tol
    while (!in_feature_ptrs.empty())
    {
      this->setProgress(initial_size - in_feature_ptrs.size());

      // get a feature with the highest Intensity
      auto candidate_fg_pointer = std::max_element(in_feature_ptrs.begin(), in_feature_ptrs.end(), FLASHQuantHelper::CmpFeatureGroupByScore());
      auto candidate_fg = *candidate_fg_pointer;

      // get all features within mass_tol from candidate FeatureGroup
      std::vector<FeatureGroup*>::iterator low_it, up_it;

      // open up the search space (10 Da) to check the mass trace overlap. mass_tolerance_da_ will be checked later
      FeatureGroup lower_fg(candidate_fg->getMonoisotopicMass() - 10);
      FeatureGroup upper_fg(candidate_fg->getMonoisotopicMass() + 10);

      low_it = std::lower_bound(in_feature_ptrs.begin(), in_feature_ptrs.end(), &lower_fg, FLASHQuantHelper::CmpFeatureGroupPointersByMass());
      up_it = std::upper_bound(in_feature_ptrs.begin(), in_feature_ptrs.end(), &upper_fg, FLASHQuantHelper::CmpFeatureGroupPointersByMass());

      // no matching in features (found only candidate itself)
      if (up_it - low_it == 1)
      {
        // save it to out_features
        if (scoreAndFilterFeatureGroup_(*candidate_fg, min_iso_score)) // rescoring is needed to set scores in FeatureGroup
        {
          candidate_fg->updateMembers();
          out_feature.push_back(*candidate_fg);
        }

        // remove candidate from features
        in_feature_ptrs.erase(candidate_fg_pointer);
        continue;
      }

      // check if found features are overlapping with the candidate feature
      std::vector<int> v_indices_to_remove;
      v_indices_to_remove.reserve(up_it - low_it);
      std::set<Size> mt_indices_to_add;
      std::vector<FeatureSeed *> mts_to_add; // One unique mt can be included in different FGs -> different FeatureSeed
      mts_to_add.reserve((up_it - low_it) * candidate_fg->size());

      for (; low_it != up_it; ++low_it)
      {
        // if low_it is candidate feature, ignore
        if (candidate_fg_pointer == low_it)
        {
          v_indices_to_remove.push_back(low_it - in_feature_ptrs.begin());
          continue;
        }

        // if the mass difference is larger than mass_tolerance_da_, check the mass trace overlap
        if ((std::abs(candidate_fg->getMonoisotopicMass() - (*low_it)->getMonoisotopicMass()) > mass_tolerance_da_)
            && !doMassTraceIndicesOverlap(**low_it, *candidate_fg))
        {
          continue;
        }

        // check if fwhm overlaps
        if (!doFWHMbordersOverlap_((*low_it)->getFwhmRange(), candidate_fg->getFwhmRange()))
        {
          continue;
        }

        // merge found feature to candidate feature
        auto &trace_indices = candidate_fg->getTraceIndices();
        for (auto &new_mt: **low_it)
        {
          // if this mass trace is not used in candidate_fg
          if (std::find(trace_indices.begin(), trace_indices.end(), new_mt.getTraceIndex()) == trace_indices.end())
          {
            mt_indices_to_add.insert(new_mt.getTraceIndex());
            mts_to_add.push_back(&new_mt);
          }
        }
        // add index of found feature to "to_be_removed_vector"
        v_indices_to_remove.push_back(low_it - in_feature_ptrs.begin());
      }

      // sort mts_to_add by abundance
      std::sort(mts_to_add.begin(), mts_to_add.end(), FLASHQuantHelper::CmpFeatureSeedByIntensity());

      // add extra masstraces to candidate_feature
      FeatureGroup final_candidate_fg(*candidate_fg); // copy of candidate_feature
      double mono_mass = final_candidate_fg.getMonoisotopicMass();
      int max_isotope = (int)iso_model_.getLastIndex(mono_mass);
      int min_isotope = (int)(iso_model_.getApexIndex(mono_mass) - iso_model_.getLeftCountFromApex(mono_mass));
      FeatureSeed* apex_of_this_fg = final_candidate_fg.getApexLMT();
      double rt_of_apex = (apex_of_this_fg->getFwhmEnd() + apex_of_this_fg->getFwhmStart())/2;
      for (auto &new_mt: mts_to_add)
      {
        // to skip duplicated masstraces that are included
        if (mt_indices_to_add.find(new_mt->getTraceIndex()) == mt_indices_to_add.end())
        {
          continue;
        }
        mt_indices_to_add.erase(new_mt->getTraceIndex());

        // if this mt is introducing new charge
        if (!final_candidate_fg.doesThisChargeExist(new_mt->getCharge()))
        {
          final_candidate_fg.push_back(*new_mt);
          continue;
        }

        /// re-calculate isotope index (from PeakGroup::recruitAllPeaksInSpectrum)
        double cmz = (mono_mass) / new_mt->getCharge() + Constants::PROTON_MASS_U; // mono mz
        const double iso_delta = iso_da_distance_ / (double)new_mt->getCharge();
        int iso_index = (int)round((new_mt->getCentroidMz() - cmz) / iso_delta);
        if (iso_index > max_isotope)
        {
          break;
        }
        if (iso_index < min_isotope)
        {
          continue;
        }

        new_mt->setIsotopeIndex(iso_index);

        // checking if this mt is colliding with the existing mt (iso index)
        bool introduce_new_mt = true;
        for (Size i = 0; i < final_candidate_fg.size(); ++i)
        {
          const FeatureSeed& seed = final_candidate_fg[i];
          if (seed.getCharge() != new_mt->getCharge())
          {
            continue;
          }
          if (seed.getIsotopeIndex() == iso_index) // if this isotope exist in feature group,
          {
            // the one with closer apex retention time to the max peak of feature group wins.
            double diff_to_new_mt = abs((new_mt->getFwhmEnd()+new_mt->getFwhmStart())/2 - rt_of_apex);
            double diff_to_org_mt = abs((seed.getFwhmEnd()+seed.getFwhmStart())/2 - rt_of_apex);
            if (diff_to_new_mt >= diff_to_org_mt)
            {
              introduce_new_mt = false;
            }
            else
            {
              final_candidate_fg.erase(final_candidate_fg.begin() + i);
            }
            break;
          }
        }
        if (introduce_new_mt)
        {
          final_candidate_fg.push_back(*new_mt);
        }
      }

      // don't merge when it failed to exceed filtering threshold // TODO: change threshold to the original iso cosine
      if (!scoreAndFilterFeatureGroup_(final_candidate_fg, min_iso_score))
      {
        if (scoreAndFilterFeatureGroup_(*candidate_fg, min_iso_score))
        {
          candidate_fg->updateMembers();
          out_feature.push_back(*candidate_fg);
        }
      }
      else // if to be merged, save the updated one to out_feature
      {
        final_candidate_fg.updateMembers();
        out_feature.push_back(final_candidate_fg);
      }

      // remove candidate from features
      std::sort(v_indices_to_remove.begin(), v_indices_to_remove.end(), std::greater<>());

      for (auto &idx : v_indices_to_remove)
      {
        in_feature_ptrs.erase(in_feature_ptrs.begin() + idx);
      }
    }
    this->endProgress();

    in_features.swap(out_feature);
  }

  void FLASHQuantAlgorithm::buildMassTraceGroups_(std::vector<FeatureSeed> &mtraces, std::vector<FeatureGroup> &features)
  {
    /// FLASHDeconvAlgorithm setting
    Param fd_defaults = getFLASHDeconvParams_();
    deconv_.setParameters(fd_defaults);
    deconv_.calculateAveragine(false);
    iso_model_ = deconv_.getAveragine();
//    std::vector<double> target_masses_; // monoisotope
//    fd_.setTargetMasses(target_masses_, ms_level);

    /// group mass traces to spectrum
    std::vector<std::pair<double, FeatureSeed*>> mt_rt_starts;
    std::vector<std::pair<double, FeatureSeed*>> mt_rt_ends;
    mt_rt_starts.reserve(mtraces.size());
    mt_rt_ends.reserve(mtraces.size());
    int counter = 0;

    // collect rt information from mtraces to generate spectrum
    double rt_binning_size = numeric_limits<double>::max();
    for (auto &trace: mtraces)
    {
      mt_rt_starts.push_back(std::make_pair(trace.getFwhmStart(), &trace));
      mt_rt_ends.push_back(std::make_pair(trace.getFwhmEnd(), &trace));
      if (trace.getMassTrace().getFWHM() < rt_binning_size)
      {
        rt_binning_size = trace.getMassTrace().getFWHM();
      }
    }
    // minimum bin size is 1
    if (rt_binning_size < 1)
    {
      rt_binning_size = 1;
    }

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

    int possible_spec_size = int((mt_rt_starts[mt_rt_starts.size() - 1].first - end_of_current_rt_window) / rt_binning_size);
    this->startProgress(0, possible_spec_size, "assembling mass traces to features");

    while (rt_s_iter != end_of_iter && end_of_current_rt_window < last_rt)
    {
      this->setProgress(counter);  // TODO; why finish at 2n%?

      // initial rt binning is 1 sec (for generating spectrum)
      end_of_current_rt_window += rt_binning_size;

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
      sort(local_traces.begin(), local_traces.end(), FLASHQuantHelper::CmpFeatureSeedByMZ());

      std::vector<FeatureGroup> local_fgroup;
      getFeatureFromSpectrum_(local_traces, local_fgroup, end_of_current_rt_window);
      ++counter; // to track the number of generated spectra
      // no feature has been detected
      if (local_fgroup.size() == 0)
      {
        continue;
      }

      for (auto &tmp_fg: local_fgroup)
      {
        sort(tmp_fg.begin(), tmp_fg.end());

        tmp_fg.updateMembersForScoring();
        tmp_fg.updateMembers();
      }

      features.insert(features.end(), local_fgroup.begin(), local_fgroup.end());
    }

    this->endProgress();
    OPENMS_LOG_INFO << "# generated spec from mass traces : " << counter << endl;
    OPENMS_LOG_INFO << "# generated feature groups from mass traces : " << features.size() << endl;
  }

  /// cluster FeatureGroups with shared mass traces, and resolve the shared ones. If not, report as output
  void FLASHQuantAlgorithm::clusterFeatureGroups_(std::vector<FeatureGroup> &fgroups,
                                                        std::vector<MassTrace> &input_mtraces)
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

    // *********************************************************** //
    // Step 2 constructing hypergraph from featuregroups
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

    // BFS
    while (true)
    {
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
        // re-scoring (score threshold changed from this method)
        auto &out_feat = fgroups[*(fg_indices_in_current_cluster.begin())];
        if (rescoreFeatureGroup_(out_feat))
        {
          out_features.push_back(out_feat);
        }
        continue;
      }

      // resolve the conflict among feature groups
      resolveConflictInCluster_(fgroups, input_mtraces, shared_m_traces, fg_indices_in_current_cluster, out_features);
    }

    out_features.shrink_to_fit();
    std::swap(out_features, fgroups);
  }

  void FLASHQuantAlgorithm::setOptionalDetailedOutput_()
  {
    String out_path = output_file_path_.substr(0, output_file_path_.find_last_of(".")) + "_shared.tsv";
    shared_out_stream_.open(out_path, std::fstream::out);
    shared_out_stream_ << "FeatureGroupID\tTraceType\tMass\tCharge\tIsotopeIndex\tQuantValue\tCentroidMz\tRTs\tMZs\tIntensities\n"; // header
    // shared = 0 (false), 1 (true, before resolution), 2(true, after resolution), 3(theoretical shape)
  }

  void FLASHQuantAlgorithm::writeMassTracesOfFeatureGroup_(const FeatureGroup &fgroup,
                                                                 const Size &fgroup_idx,
                                                                 const std::vector<std::vector<Size>> &shared_m_traces_indices,
                                                                 const bool &is_before_resolution)
  {
    /// shared_type = 0 (before resolution), 2 (after resolution)
    auto mt_idxs = fgroup.getTraceIndices();
    Size shared_tag_for_output = 0;
    if (!is_before_resolution)
    {
      shared_tag_for_output = 2;
    }

    for (const auto &mt: fgroup)
    {
      // check if current mt is shared with other features
      if (is_before_resolution) // before resolution
      {
        if (shared_m_traces_indices[mt.getTraceIndex()].size() == 1) // unique mass trace
        {
          shared_tag_for_output = 0;
        }
        else
        {
          shared_tag_for_output = 1;
        }
      }

      stringstream rts;
      stringstream mzs;
      stringstream intys;

      for (const auto &peak: mt.getMassTrace())
      {
        mzs << std::to_string(peak.getMZ()) << ",";
        rts << std::to_string(peak.getRT()) << ",";
        intys << std::to_string(peak.getIntensity()) << ",";
      }
      std::string peaks = rts.str();
      peaks.pop_back();
      peaks = peaks + "\t" + mzs.str();
      peaks.pop_back();
      peaks = peaks + "\t" + intys.str();
      peaks.pop_back();

      shared_out_stream_ << fgroup_idx << "\t"
                         << shared_tag_for_output << "\t"
                         << std::to_string(fgroup.getMonoisotopicMass()) << "\t"
                         << mt.getCharge() << "\t"
                         << mt.getIsotopeIndex() << "\t"
                         << std::to_string(mt.getIntensity()) << "\t"
                         << std::to_string(mt.getCentroidMz()) << "\t"
                         << peaks + "\n";
    }
  }

  void FLASHQuantAlgorithm::writeTheoreticalShapeForConflictResolution_(const Size &fgroup_idx,
                                                                              const FeatureSeed &shared_mt,
                                                                              const std::vector<double> &theo_intensities,
                                                                              const double &calculated_ratio)
  {
    stringstream rts;
    stringstream mzs;
    stringstream intys;

    for (const auto &peak: shared_mt.getMassTrace())
    {
      mzs << std::to_string(peak.getMZ()) << ",";
      rts << std::to_string(peak.getRT()) << ",";
    }
    for (const auto &inty : theo_intensities)
    {
      intys << std::to_string(inty) << ",";
    }

    std::string peaks = rts.str();
    peaks.pop_back();
    peaks = peaks + "\t" + mzs.str();
    peaks.pop_back();
    peaks = peaks + "\t" + intys.str();
    peaks.pop_back();

    shared_out_stream_ << fgroup_idx << "\t"
                       << 3 << "\t" // shared tag = 3
                       << std::to_string(shared_mt.getMass()) << "\t"
                       << shared_mt.getCharge() << "\t"
                       << shared_mt.getIsotopeIndex() << "\t"
                       << std::to_string(calculated_ratio) << "\t" // ratio
                       << std::to_string(shared_mt.getCentroidMz()) << "\t"
                       << peaks + "\n";
  }

  void FLASHQuantAlgorithm::filterOutIneligibleFeatureGroupsInCluster(std::vector<FeatureGroup>& feature_groups,
                                                                            std::vector<std::vector<Size>>& shared_m_traces_indices,
                                                                            std::set<Size>& candidate_fg_indices) const
  {
    /// Ineligible feature groups = harmonics OR the ones without unique feature groups
    std::set<Size> fg_idx_to_remove;
    for (auto fg_idx = candidate_fg_indices.begin(); fg_idx != candidate_fg_indices.end(); ++fg_idx)
    {
      if (std::find(fg_idx_to_remove.begin(), fg_idx_to_remove.end(), *fg_idx) != fg_idx_to_remove.end()) // if this is already removed
      {
        continue;
      }

      auto& fgroup = feature_groups[*fg_idx];

      // check if any masstrace is unique
      bool is_this_all_sharing = true;
      for (auto& seed_iter : fgroup)
      {
        auto trace_idx = seed_iter.getTraceIndex();
        if (shared_m_traces_indices[trace_idx].size() == 1) // this masstrace is unique
        {
          is_this_all_sharing = false;
          break;
        }
      }
      if (is_this_all_sharing)
      {
        fg_idx_to_remove.insert(*fg_idx);
        continue;
      }

      // check if this is harmonics to the other
      std::vector<double> harmonic_masses; // harmonic mass candidates
      const double fg_mass = fgroup.getMonoisotopicMass();
      for (int harmonic_ratio = 2; harmonic_ratio < 8; ++harmonic_ratio)
      {
        for (int iso = -2; iso < 3; ++iso)
        {
          // high harmonics
          double h_mass = (fg_mass + iso*iso_da_distance_) / harmonic_ratio;
          if (h_mass > min_mass_)
          {
            harmonic_masses.push_back(h_mass);
          }
          // low harmonic
          h_mass = (fg_mass + iso*iso_da_distance_) * harmonic_ratio;
          if (h_mass < max_mass_)
          {
            harmonic_masses.push_back(h_mass);
          }
        }
      }
      std::sort(harmonic_masses.begin(), harmonic_masses.end());

      for (auto fg_idx_to_compare = std::next(fg_idx); fg_idx_to_compare != candidate_fg_indices.end(); ++fg_idx_to_compare)
      {
        if (std::find(fg_idx_to_remove.begin(), fg_idx_to_remove.end(), *fg_idx_to_compare) != fg_idx_to_remove.end()) // if this is already removed
        {
          continue;
        }

        auto& fgroup_to_compare = feature_groups[*fg_idx_to_compare];
        const double mass_tolerance = fgroup_to_compare.getMonoisotopicMass() * 2 *  1e-5; // 20 ppm

        // is harmonic or not
        auto low_it = std::lower_bound(harmonic_masses.begin(), harmonic_masses.end(), fgroup_to_compare.getMonoisotopicMass() - mass_tolerance); // inclusive
        auto up_it = std::upper_bound(harmonic_masses.begin(), harmonic_masses.end(), fgroup_to_compare.getMonoisotopicMass() + mass_tolerance); // exclusive

        if (low_it == up_it) // if nothing is found, continue
        {
          continue;
        }

        Size harmonic_fg_idx = fgroup.getIntensity() < fgroup_to_compare.getIntensity() ? *fg_idx : *fg_idx_to_compare;
        if (doMassTraceIndicesOverlap(fgroup, fgroup_to_compare, 0.5, false))
        {
          fg_idx_to_remove.insert(harmonic_fg_idx);
        }
        else
        {
          // how many mass traces are shared?
          Size num_of_shared_mt = 0;
          auto &fg_to_remove = feature_groups[harmonic_fg_idx];
          for (auto &mt : fg_to_remove)
          {
            auto trace_idx = mt.getTraceIndex();
            if (shared_m_traces_indices[trace_idx].size() > 1) // this masstrace is shared
            {
              num_of_shared_mt++;
            }
          }

          // if the number of unique mt are less than the range of charges, remove this feature group
          if ((fg_to_remove.size()-num_of_shared_mt) < fg_to_remove.getChargeSet().size())
          {
            fg_idx_to_remove.insert(harmonic_fg_idx);
          }
        }
      }
    }

    // remove feature groups for further work & update the mass trace link
    for (auto &idx : fg_idx_to_remove)
    {
      // remove feature groups for further work
      candidate_fg_indices.erase(idx);
//      OPENMS_LOG_INFO << idx << " ";

      // update the mass link
      for (auto& seed_iter : feature_groups[idx])
      {
        auto trace_idx = seed_iter.getTraceIndex();
        if (shared_m_traces_indices[trace_idx].size() > 1) // update the link
        {
          auto &link = shared_m_traces_indices[trace_idx];
          link.erase(std::find(link.begin(), link.end(), idx));
        }
      }
    }
  }

  void FLASHQuantAlgorithm::resolveConflictInCluster_(std::vector<FeatureGroup> &feature_groups,
                                                            std::vector<MassTrace> &input_masstraces,
                                                            std::vector<std::vector<Size>> &shared_m_traces_indices,
                                                            std::set<Size> &fg_indices_in_this_cluster,
                                                            std::vector<FeatureGroup> &out_featuregroups)
  {
    /// conflict resolution is done in feature level (not feature group level), starting from the most abundant shared signal
    /// input masstraces needs to be changed, if any resolution is done per FeatureGroup -> add new input_masstrace with modified intensities, and change the pointer to it.

    /// remove FeatureGroups without unique masstraces (+ harmonics)
    filterOutIneligibleFeatureGroupsInCluster(feature_groups, shared_m_traces_indices, fg_indices_in_this_cluster);
    if (fg_indices_in_this_cluster.size() == 1)
    {
      auto &fg = feature_groups[*(fg_indices_in_this_cluster.begin())];
      if (rescoreFeatureGroup_(fg))
      {
        out_featuregroups.push_back(fg);
      }
      return;
    }

    /// 1. find shared MTs and corresponding "Features" (not FeatureGroup)
    std::map<Size, std::vector<std::pair<Size, int>>> shared_mt_N_features; // key: mt_idx, value: pair<feature_group_idx, charge>
    std::vector<std::pair<Size, double>> shared_mt_idx_inty_pairs; // first: mt_idx, second: mt intensity
    for (auto &fg_idx : fg_indices_in_this_cluster)
    {
      auto &fgroup = feature_groups[fg_idx];
      for (auto &seed_iter : fgroup)
      {
        auto trace_idx = seed_iter.getTraceIndex();
        // collect only shared one. Skip unique traces
        if (shared_m_traces_indices[trace_idx].size() < 2)
        {
          continue;
        }

        if (shared_mt_N_features.find(trace_idx) == shared_mt_N_features.end()) // if this seed is not added to the map
        {
          shared_mt_N_features.emplace(trace_idx, std::vector<std::pair<Size, int>>());
          shared_mt_idx_inty_pairs.push_back(std::make_pair(trace_idx, seed_iter.getIntensity()));
        }
        shared_mt_N_features[trace_idx].push_back(std::make_pair(fg_idx, seed_iter.getCharge()));
      }

      if(shared_output_requested_)
      {
        writeMassTracesOfFeatureGroup_(fgroup, fg_idx, shared_m_traces_indices, true);
      }
    }

    /// 2. sort shared MT by intensity (descending)
    std::sort(shared_mt_idx_inty_pairs.begin(), shared_mt_idx_inty_pairs.end(), [](const auto &left, const auto &right) {
      return left.second > right.second;
    });

    /// 3. start from the highest intensity, iterate shared MT and get conflicting region (collecting conflicting Features)
    while(!shared_mt_idx_inty_pairs.empty())
    {
      /// 3.1 collecting all shared mts and features in this conflicting region
      std::unordered_map<Size, bool> shared_mts_in_this_region; // to track all detected shared mass traces. key = shared_mt_idx, value = visited?
      std::set<String> visited_features_in_this_region; // to track all detected features. key = [FeatGroupID]&[CS]
      std::vector<Feature> conflicting_features;
      std::vector<Feature> features_not_for_resolution;
      std::vector<MassTrace> conflicting_mts; // collect shared traces for resolution (not FeatureSeed, but original MassTrace. FeatureSeed has conflicting information)
      std::vector<Size> conflicting_mt_indices; // trace_index of collect shared traces for resolution

      shared_mts_in_this_region.emplace(shared_mt_idx_inty_pairs[0].first, false); // add the first shared mt
      while(true)
      {
        // find unvisited shared mt
        auto found_mt = std::find_if(std::begin(shared_mts_in_this_region), std::end(shared_mts_in_this_region),
                               [](auto&& p) { return p.second == false; });
        // if no shared mt left
        if (found_mt == std::end(shared_mts_in_this_region))
        {
          break;
        }

        Size this_shared_mt_index = found_mt->first;
        shared_mts_in_this_region[this_shared_mt_index] = true; // mark as visited

        // per conflicting features, make Feature
        for (auto &feat : shared_mt_N_features[this_shared_mt_index])
        {
          FeatureGroup &org_featgroup = feature_groups[feat.first];
          int cs_of_feat = feat.second;

          // check if FeatureElement is not created yet
          String feat_label = std::to_string(feat.first) + '&' + std::to_string(cs_of_feat);
          if (visited_features_in_this_region.find(feat_label) != visited_features_in_this_region.end())
          {
            continue;
          }
          visited_features_in_this_region.insert(feat_label); // mark as this Feature is visited

          // initialize Feature
          Feature new_feature;
          new_feature.prepareVectors(org_featgroup.size());
          new_feature.charge = cs_of_feat;
          new_feature.feature_group_index = feat.first;

          // add seeds in element
          for (auto &feat_seed: org_featgroup)
          {
            if (feat_seed.getCharge() != cs_of_feat)
            {
              continue;
            }

            if (shared_m_traces_indices[feat_seed.getTraceIndex()].size() > 1) // this seed is shared
            {
              new_feature.shared_traces.push_back(feat_seed);
              new_feature.shared_trace_indices.push_back(feat_seed.getTraceIndex());
            }
            else // this seed is not shared
            {
              new_feature.unique_traces.push_back(feat_seed);
              new_feature.unique_trace_indices.push_back(feat_seed.getTraceIndex());
            }
          }

          // check if this Feature is eligible for resolution
          if (!isEligibleFeatureForConflictResolution_(new_feature, shared_m_traces_indices, org_featgroup))
          {
            // put this feature in the list
            features_not_for_resolution.push_back(new_feature);
            continue;
          }
          new_feature.shrinkVectors();

          // if any seed is shared, mark for visiting in the next round
          for (auto &mt_idx : new_feature.shared_trace_indices)
          {
            if (shared_mts_in_this_region.find(mt_idx) == shared_mts_in_this_region.end()) // if this shared trace is not visited yet
            {
              shared_mts_in_this_region.emplace(mt_idx, false);
            }
          }
          conflicting_features.push_back(new_feature);
        }
      }

      /// 3.2 remove all conflicting MTs from "look-up table"
      for (auto &shared : shared_mts_in_this_region)
      {
        shared_mt_idx_inty_pairs.erase(std::remove_if(shared_mt_idx_inty_pairs.begin(), shared_mt_idx_inty_pairs.end(),
                                                      [shared](auto const &p) {return p.first == shared.first;}));
      }

      /// 3.3 remove invalid feature from FeatureGroup and then collect conflicting mass traces
      for (auto &feat : features_not_for_resolution)
      {
        auto &feat_group = feature_groups[feat.feature_group_index];
        auto &feat_cs = feat.charge;
        auto seed_iter = feat_group.begin();
        while (seed_iter != feat_group.end()) {
          if (seed_iter->getCharge() == feat_cs) {
            // if this seed is from the "feat_not_for_resolution", remove it.
            seed_iter = feat_group.erase(seed_iter);
          }
          else {
            ++seed_iter;
          }
        }
        // only update members for scoring. (FeatureGroup can be empty at this point)
        feat_group.updateMembersForScoring();

        // update the lookup table to collect shared masstrace to be updated later
        auto pair_to_remove = std::make_pair(feat.feature_group_index, feat_cs);
        for (auto &shared_idx : feat.shared_trace_indices)
        {
          auto &feat_vec = shared_mt_N_features[shared_idx];
          feat_vec.erase(std::find(feat_vec.begin(), feat_vec.end(), pair_to_remove));
        }
      }

      // collect shared mass trace to be resolved
      for (auto &candidate_mt : shared_mts_in_this_region)
      {
        Size shared_mt_idx = candidate_mt.first;
        auto &conflict_feat_tags = shared_mt_N_features[shared_mt_idx];

        // if this is not used by any features
        if (conflict_feat_tags.empty())
        {
          // mark the shared mt as not shared
          shared_m_traces_indices[shared_mt_idx] = std::vector<Size>();
          continue;
        }
        // check if this mt is not shared one anymore
        else if (conflict_feat_tags.size() == 1)
        {
          Size fg_id = conflict_feat_tags[0].first; // feature group id
          int f_charge = conflict_feat_tags[0].second;

          // mark this MT is not shared one
          shared_m_traces_indices[shared_mt_idx] = std::vector<Size>{fg_id};

          // if the corresponding feature is already added to conflict_features, update it
          auto found_feat = find_if(conflicting_features.begin(), conflicting_features.end(),
                                    [fg_id, f_charge](auto &&x) { return (x.feature_group_index == fg_id) && (x.charge == f_charge); });
          if (found_feat == conflicting_features.end())
          {
            continue;
          }
          // find the shared trace
          auto it = std::find(found_feat->shared_trace_indices.begin(), found_feat->shared_trace_indices.end(), shared_mt_idx);
          Size trace_loc = it - found_feat->shared_trace_indices.begin();

          // move the shared trace to unique
          found_feat->unique_trace_indices.push_back(shared_mt_idx);
          FeatureSeed the_trace = found_feat->shared_traces[trace_loc];
          found_feat->unique_traces.push_back(the_trace);

          // delete the trace from the shared trace vectors
          found_feat->shared_trace_indices.erase(found_feat->shared_trace_indices.begin()+trace_loc);
          found_feat->shared_traces.erase(found_feat->shared_traces.begin()+trace_loc);

          // check if any unique_trace is left
          if (found_feat->shared_traces.empty())
          {
            // remove from conflicting features
            conflicting_features.erase(found_feat);
          }
          continue;
        }
        // this trace is shared one -> add this shared mt as a conflicting one
        auto i = input_masstraces.begin() + shared_mt_idx;
        conflicting_mts.push_back(*i);
        conflicting_mt_indices.push_back(shared_mt_idx);
      }
      // if no resolving conflict is needed
      if (conflicting_features.size() < 2)
      {
        continue;
      }

      /// 3.4 set isotope probabilities for features (to be used as theoretical intensity when modeling a theoretical shape)
      double minimum_probability = 1.0;
      for (auto &feat : conflicting_features)
      {
        auto &feat_group = feature_groups[feat.feature_group_index];
        auto theo_isodist = iso_model_.get(feat_group.getMonoisotopicMass());
        feat.isotope_probabilities.reserve(feat.unique_traces.size()); // isotope prob for each mass_trace, based on their iso position
        for (auto &lmt: feat.unique_traces)
        {
          // if this trace is taken from different Feature
          if (lmt.getCharge() != feat.charge)
          {
            feat.isotope_probabilities.push_back(-1); // to be set later (as the minimum value)
            continue;
          }

          // if isotope index of lmt exceed tmp_iso length, give 0 (to not use for the modeling)
          Size tmp_iso_idx = (Size) lmt.getIsotopeIndex();
          if (tmp_iso_idx < theo_isodist.size() && tmp_iso_idx >= 0) // To avoid bad access error. If index is larger than the isodist length, no intensity to use.
          {
            // give weight per masstrace based on their isotope position, not the real intensity.
            double tmp_prob = theo_isodist[tmp_iso_idx].getIntensity();
            feat.isotope_probabilities.push_back(tmp_prob);
            if (tmp_prob < minimum_probability)
            {
              minimum_probability = tmp_prob;
            }
          }
          else
          {
            feat.isotope_probabilities.push_back(0.0);
          }
        }
      }
      // switch -1 probability to minimum prob
      minimum_probability = minimum_probability < 0 ? 0 : minimum_probability;
      for (auto &feat : conflicting_features)
      {
        for (auto &p: feat.isotope_probabilities)
        {
          if (p == -1)
          {
            p = minimum_probability;
          }
        }
      }

      /// 3.5 resolve the conflict in this region
      resolveConflictRegion_(conflicting_features, conflicting_mts, conflicting_mt_indices);

      /// 3.6 Update FeatureGroup (only FeatureSeed vectors, other members will be updated later)
      for (auto &feat : conflicting_features)
      {
        // get original FeatureGroup, and update it
        FeatureGroup &feature_group = feature_groups[feat.feature_group_index];

        for (Size t_index = 0; t_index < feat.shared_trace_indices.size(); ++t_index)
        {
          auto &updated_seed = feat.shared_traces[t_index];
          Size &index_of_updated_seed = feat.shared_trace_indices[t_index];

          // find corresponding FeatureSeed in FeatureGroup
          auto found_seed = find_if(feature_group.begin(), feature_group.end(),
                                    [index_of_updated_seed](auto &&x) { return x.getTraceIndex() == index_of_updated_seed; });

          if (updated_seed.getIntensity() == found_seed->getIntensity()) // not changed
          {
            continue;
          }

          // if updated seed is empty, add remove found seed
          if (updated_seed.getIntensity() == 0)
          {
            feature_group.erase(found_seed);
            continue;
          }

          *found_seed = updated_seed;
          found_seed->setTraceIndex(input_masstraces.size()); // add new trace index
          // add new FeatureSeed
          input_masstraces.push_back(updated_seed.getMassTrace());
        }

        // add theoretical shape information
        feature_group.updateTheoreticalShapes(feat.theoretical_shapes);
      }
    }

    /// 4. update FeatureGroup members and quantities
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

      if (rescoreFeatureGroup_(fgroup))
      {
        out_featuregroups.push_back(fgroup);
        if (shared_output_requested_)
        {
          writeMassTracesOfFeatureGroup_(fgroup, idx, shared_m_traces_indices, false);
        }
      }
    }
  }

  bool FLASHQuantAlgorithm::isEligibleFeatureForConflictResolution_(Feature &feat, std::vector<std::vector<Size>> &shared_m_traces_indices, FeatureGroup &feat_group) const
  {
    // check if this feature is composed of only shared mass traces
    if (feat.unique_traces.size() > 0)
    {
      return true;
    }

    // get the most abundant mass trace from FeatureGroup (except recruited ones)
    FeatureSeed *most_abundant_mt = nullptr;
    getMostAbundantMassTraceFromFeatureGroup_(feat_group, feat.charge, most_abundant_mt, shared_m_traces_indices);
    if (most_abundant_mt == nullptr)
    {
      return false;
    }

    // get the most abundant mass trace from FeatureGroup (except recruited ones)
    feat.unique_traces.push_back(*most_abundant_mt);
    feat.unique_trace_indices.push_back(most_abundant_mt->getTraceIndex());
    return true;
  }

  void FLASHQuantAlgorithm::resolveConflictRegion_(std::vector<Feature> &conflicting_features,
                                                         const std::vector<MassTrace> &conflicting_mts,
                                                         const std::vector<Size> &conflicting_mt_indices)
  {
    // if only one feature has been passed, skip.
    if (conflicting_features.size() < 2)
    {
      return;
    }

    /// Prepare Components per features (excluding conflicting mts)
    std::vector<std::vector<double>> components;
    Matrix<int> pointer_to_components; // row : index of conflicting_mts , column : index of conflicting_features
    pointer_to_components.resize(conflicting_mts.size(), conflicting_features.size(), -1);
    for (Size i_of_f = 0; i_of_f < conflicting_features.size(); ++i_of_f)
    {
      auto &tmp_feat = conflicting_features[i_of_f];
      EGHTraceFitter* fitted_model = new EGHTraceFitter();
      fitTraceModelFromUniqueTraces_(tmp_feat, fitted_model);

      // the model is not valid
      double area = fitted_model->getArea();
      if (area != area) // x != x: test for NaN
      {
        continue;
      }

      // store Component information (Component = calculated model of the trace)
      for (Size row = 0; row < conflicting_mts.size(); ++row)
      {
        // if the feature doesn't own this conflicting_mt, skip
        auto &trace_indices = tmp_feat.shared_trace_indices;
        if (trace_indices.end() == std::find(trace_indices.begin(), trace_indices.end(), conflicting_mt_indices[row]))
        {
          continue;
        }

        // normalize fitted value
        std::vector<double> fit_intensities;
        double summed_intensities = .0;
        for (auto &peak: conflicting_mts[row])
        {
          double rt = peak.getRT();
          double fitted_value = fitted_model->getValue(rt);
          summed_intensities += fitted_value;
          fit_intensities.push_back(fitted_value);
        }

        // if fit_model's RT range doesn't overlap with the conflict_mt -> fit_intensities vec is all composed of 0.
        if (summed_intensities == 0)
        {
          continue;
        }

        // save normalized intensities into component
        std::vector<double> component;
        for (auto &inty: fit_intensities)
        {
          component.push_back(inty/summed_intensities);
        }

        pointer_to_components.setValue(row, i_of_f, components.size()); // set the index to this component
        components.push_back(component);
      }
    }

    /// reconstruction of observed XICs (per conflicting mt)
    for (Size row = 0; row < conflicting_mts.size(); ++row)
    {
      Size org_index_of_this_trace = conflicting_mt_indices[row];
      auto &obs_masstrace = conflicting_mts[row];

      auto components_indices = pointer_to_components.row(row);
      Size column_size =
          components_indices.size() - std::count(components_indices.begin(), components_indices.end(), -1); // number of features involved with this mt

      // if no resolution is needed, skip the next part
      if (column_size < 2) // only one or no Feature has valid model
      {
        for (Size i_of_f = 0; i_of_f < components_indices.size(); ++i_of_f) // iterate feature
        {
          if (components_indices[i_of_f] == -1) // if this feature is not the one with the valid model
          {
            Feature &feat = conflicting_features[i_of_f];
            auto lmt_ptr = std::find_if(feat.shared_traces.begin(), feat.shared_traces.end(),
                                        [org_index_of_this_trace](auto &&x) { return x.getTraceIndex() == org_index_of_this_trace; });
            if (lmt_ptr != feat.shared_traces.end())
            {
              lmt_ptr->setIntensity(0); // to be removed later
              if (shared_output_requested_)
              {
                std::vector<double> zero_vec = std::vector<double>(lmt_ptr->getMassTrace().getSize(), 0);
                writeTheoreticalShapeForConflictResolution_(feat.feature_group_index, *lmt_ptr, zero_vec, 0);
                feat.setTheoreticalShapes(*lmt_ptr, zero_vec, 0);
              }
            }
          }
          else // valid Feature
          {
            if (shared_output_requested_)
            {
              Feature &feat = conflicting_features[i_of_f];
              auto lmt_ptr = std::find_if(feat.shared_traces.begin(), feat.shared_traces.end(),
                                          [org_index_of_this_trace](auto &&x) { return x.getTraceIndex() == org_index_of_this_trace; });
              auto &temp_comp = components[pointer_to_components.getValue(row, i_of_f)];
              writeTheoreticalShapeForConflictResolution_(feat.feature_group_index, *lmt_ptr, temp_comp, 1);
              feat.setTheoreticalShapes(*lmt_ptr, temp_comp, 1);
            }
          }
        }

        continue;
      }
      updateFeatureWithFitModel(conflicting_features, row, obs_masstrace, org_index_of_this_trace, pointer_to_components, components);
    }
    components.clear();
  }

  void FLASHQuantAlgorithm::updateFeatureWithFitModel(std::vector<Feature>& conflicting_features, Size mt_index,
                                                            const MassTrace& obs_masstrace, const Size& org_index_of_obs_mt,
                                                            Matrix<int>& pointer_to_components, vector<std::vector<double>>& components)
  {
    auto components_indices = pointer_to_components.row(mt_index);
    Size column_size =
      components_indices.size() - std::count(components_indices.begin(), components_indices.end(), -1); // number of features involved with this mt

    // prepare observed XIC vector
    Size mt_size = obs_masstrace.getSize();
    double obs_total_intensity = .0;
    for (auto &peak : obs_masstrace)
    {
      obs_total_intensity += peak.getIntensity();
    }
    Matrix<double> obs; // TODO: smooth intensity?
    obs.resize(mt_size, 1);
    for (Size i = 0; i < mt_size; ++i)
    {
      obs.setValue(i, 0, obs_masstrace[i].getIntensity()/obs_total_intensity); // save normalized value
    }

    // prepare theoretical matrix (include only related features)
    Matrix<double> theo_matrix;
    theo_matrix.resize(mt_size, column_size);
    Size col = 0;
    for (Size comp_idx = 0; comp_idx < components_indices.size(); ++comp_idx)
    {
      // skipping no-feature column
      if (components_indices[comp_idx] == -1)
        continue;

      auto &temp_comp = components[pointer_to_components.getValue(mt_index, comp_idx)];
      for (Size tmp_r = 0; tmp_r < temp_comp.size(); ++tmp_r)
      {
        theo_matrix.setValue(tmp_r, col, temp_comp[tmp_r]);
      }
      col++;
    }

    Matrix<double> out_quant;
    out_quant.resize(column_size, 1);
    NonNegativeLeastSquaresSolver::solve(theo_matrix, obs, out_quant);

    // if any out_quant is zero, give the other group all.
    vector<double> calculated_ratio;
    Size zero_ratio_counter = 0;
    for (Size i = 0; i < out_quant.rows(); ++i)
    {
      calculated_ratio.push_back(out_quant.getValue(i, 0));
      if (out_quant.getValue(i, 0) == 0)
      {
        ++zero_ratio_counter;
      }
    }
    if (calculated_ratio.size() - zero_ratio_counter == 1) // only one out_quant was non-zero
    {
      for (auto &ratio:calculated_ratio)
      {
        if (ratio != 0)
        {
          ratio = 1;
        }
      }
    }

    /// update Features based on the calculated ratio
    col = 0;
    for (Size i_of_f = 0; i_of_f < components_indices.size(); ++i_of_f) // iterate feature
    {
      // skipping not-this-feature column
      if (components_indices[i_of_f] == -1)
        continue;

      auto &feat = conflicting_features[i_of_f];
      auto lmt_ptr = std::find_if(feat.shared_traces.begin(), feat.shared_traces.end(),
                                  [org_index_of_obs_mt](auto &&x) { return x.getTraceIndex() == org_index_of_obs_mt; });

      // Kyowon's advice! ratio should be applied to real intensity, not theoretical one
      lmt_ptr->setIntensity(calculated_ratio[col] * obs_total_intensity);

      if (this->shared_output_requested_)
      {
        this->writeTheoreticalShapeForConflictResolution_(feat.feature_group_index, *lmt_ptr, theo_matrix.col(col), calculated_ratio[col]);
        feat.setTheoreticalShapes(*lmt_ptr, theo_matrix.col(col), calculated_ratio[col]);
      }

      /// Update MassTrace itself
      MassTrace updated = this->updateMassTrace_(lmt_ptr->getMassTrace(), calculated_ratio[col]);
      lmt_ptr->setMassTrace(updated);

      col++;
    }
  }

  // modified ElutionModelFitter::fitElutionModels
  void FLASHQuantAlgorithm::fitTraceModelFromUniqueTraces_(FLASHQuantHelper::Feature const& tmp_feat, EGHTraceFitter* fitter) const
  {
    Size mass_traces_size = tmp_feat.unique_traces.size();

    // preparation for ElutionModelFit
    FeatureFinderAlgorithmPickedHelperStructs::MassTraces traces_for_fitting;
    traces_for_fitting.reserve(mass_traces_size);
    traces_for_fitting.max_trace = 0;
    vector<Peak1D> peaks;
    // reserve space once, to avoid copying and invalidating pointers:
    peaks.reserve(tmp_feat.getPeakSizes());

    // get theoretical information from unique mass traces
    for (Size idx = 0; idx < mass_traces_size; ++idx)
    {
      auto &mass_trace_ptr = tmp_feat.unique_traces[idx].getMassTrace();

      // prepare for fitting
      FeatureFinderAlgorithmPickedHelperStructs::MassTrace tmp_mtrace;
      Size m_trace_size = mass_trace_ptr.getSize();
      tmp_mtrace.peaks.reserve(m_trace_size);
      for (auto &p_2d: mass_trace_ptr)
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

    // ElutionModelFit
    // TODO : is this necessary? giving isotope probability as a weight
    Param params = fitter->getDefaults();
    params.setValue("weighted", "true");
    fitter->setParameters(params);
    runElutionModelFit_(traces_for_fitting, fitter);
  }

  // update MassTrace
  MassTrace FLASHQuantAlgorithm::updateMassTrace_(const MassTrace& ref_trace, const double &ratio) const
  {
    /// reference = ElutionPeakDetection::detectElutionPeaks_ (step3)
    // if ratio = 0, return empty MassTrace
    if (ratio == 0)
    {
      return MassTrace();
    }

    std::vector<PeakType> tmp_mt;
    std::vector<double> smoothed_tmp;
    for (Size p_index = 0; p_index < ref_trace.getSize(); ++p_index)
    {
      auto &peak = ref_trace[p_index];
      PeakType new_peak(peak);
      new_peak.setIntensity(peak.getIntensity() * ratio);
      tmp_mt.push_back(new_peak);
      smoothed_tmp.push_back(ref_trace.getSmoothedIntensities()[p_index] * ratio);
    }

    // create new mass trace
    MassTrace new_trace(tmp_mt);
    new_trace.setSmoothedIntensities(smoothed_tmp);
    // set label as if it's a sub-trace
    new_trace.setLabel(ref_trace.getLabel() + "_" + String(ratio));
    new_trace.updateSmoothedMaxRT();
    new_trace.updateWeightedMeanMZ();
    new_trace.updateWeightedMZsd();
    new_trace.setQuantMethod(ref_trace.getQuantMethod());
    new_trace.estimateFWHM(use_smoothed_intensities_);

    return new_trace;
  }

  // from ElutionModelFitter.
  void FLASHQuantAlgorithm::runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces,
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

    // Fitting
    try
    {
      fitter->fit(m_traces);
    }
    catch (Exception::UnableToFit &except)
    {
      OPENMS_LOG_ERROR << "Error fitting model to feature '"
                       << except.getName()
                       << " - " << except.getMessage() << endl;
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

  void FLASHQuantAlgorithm::getMostAbundantMassTraceFromFeatureGroup_(const FeatureGroup &fgroup,
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
      if (lmt->getTraceIndex() < shared_m_traces.size() && shared_m_traces[lmt->getTraceIndex()].size() > 1)
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

  void FLASHQuantAlgorithm::getFLASHDeconvConsensusResult()
  {
    with_target_masses_ = true;
    std::fstream fin;
    fin.open("/Users/jeek/Documents/A4B/FDQ/Kiel-Human/FDQ_target_from_rep9.csv", std::ios::in);
//    fin.open("/Users/jeek/Documents/A4B/FDQ/Kiel-Human/ProSightPD_concensus.csv", std::ios::in);

    std::string line;
    std::getline(fin, line); // skip the first line
    while (std::getline(fin, line))
    {
      Size comma_loc = line.find(',');
      target_masses_.push_back(std::make_pair(stod(line.substr(0, comma_loc)), stod(line.substr(comma_loc + 1))));
    }

    std::sort(target_masses_.begin(), target_masses_.end());
    fin.close();
  }

  bool FLASHQuantAlgorithm::isThisMassOneOfTargets_(const double &candi_mass, const double &candi_rt) const
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