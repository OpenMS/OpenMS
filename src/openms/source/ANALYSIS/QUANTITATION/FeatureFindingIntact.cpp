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

#include <OpenMS/ANALYSIS/QUANTITATION/FeatureFindingIntact.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <queue>

// test purpose
#include <iostream>
#include <fstream>
#include <OpenMS/FORMAT/HANDLERS/GraphMLHandler.h>

namespace OpenMS
{
  FeatureFindingIntact::FeatureFindingIntact():
      ProgressLogger()
  {
    this->setLogType(CMD);
  }

  FeatureFindingIntact::~FeatureFindingIntact(){}

  void getCoordinatesOfFeaturesForPython(std::vector<FeatureFindingIntact::FeatureHypothesis> hypo)
  {
    String out_path = "/Users/jeek/Documents/A4B_UKE/filg/rp/centroid/190904_Proswift50mm_7min_Filg_500ng.pphr.featHypo.tsv";
    ofstream out;
    out.open(out_path, ios::out);

    // header
    out << "feature_label\tcs\tscore\tbounding_box_pos\tbounding_box_width\tbounding_box_height\t"
           "iso_position\tmasstrace_centroid_rts\tmasstrace_centroid_mzs\n";

    for (const auto& feat : hypo)
    {
      double mz_upper_limit = 0.0;
      double rt_upper_limit  = 0.0;
      double mz_lower_limit = std::numeric_limits<double>::max();
      double rt_lower_limit = std::numeric_limits<double>::max();
      stringstream rts;
      stringstream mzs;
      for (const auto& mt : feat.getMassTraces())
      {
        const auto& boundingbox = mt->getConvexhull().getBoundingBox();

        if (boundingbox.minX() < rt_lower_limit)
          rt_lower_limit = boundingbox.minX();
        if (boundingbox.maxX() > rt_upper_limit)
          rt_upper_limit = boundingbox.maxX();
        if (boundingbox.minY() < mz_lower_limit)
          mz_lower_limit = boundingbox.minY();
        if (boundingbox.maxY() > mz_upper_limit)
          mz_upper_limit = boundingbox.maxY();

        mzs << mt->getCentroidMZ() << ",";
        rts << mt->getCentroidRT() << ",";
      }
      std::string centroids = rts.str();
      centroids.pop_back();
      centroids = centroids + "\t" + mzs.str();
      centroids.pop_back();

      stringstream isos;
      for (const auto& pair : feat.getIndicesOfMassTraces())
      {
        isos << pair.first << ",";
      }
      std::string iso_str = isos.str();
      iso_str.pop_back();

      String label = std::to_string(feat.getFeatureMass()) + "\t" + std::to_string(feat.getCharge());
      out << label << "\t" << to_string(feat.getScore()) << "\t"
          << to_string(rt_lower_limit) << "," << to_string(mz_lower_limit) << "\t"
          << to_string((rt_upper_limit-rt_lower_limit)) << "\t" << to_string((mz_upper_limit-mz_lower_limit)) << "\t"
          << iso_str << "\t" << centroids + "\n";
    }
    out.close();
  }

  void drawSharedMasstracesBetweenFeatures(vector<FeatureFindingIntact::FeatureHypothesis>& hypotheses,
                                           const std::vector<std::vector<Size>>& shared_m_traces_indices)
  {
    String out_path = "/Users/jeek/Documents/A4B_UKE/filg/rp/centroid/190904_Proswift50mm_7min_Filg_500ng.pphr.featHypoCluster.graphML";

    std::vector<Size> hypo_nodes;
//    std::vector<std::pair<Size, Size>> edges;
    std::map<std::pair<Size, Size>, Size> edge_map;

    for (Size i = 0; i < hypotheses.size(); ++i)
    {
      hypo_nodes.push_back(i);
    }

    for (auto& feat_vec : shared_m_traces_indices)
    {
      if (feat_vec.size()==0) continue;

      for (Size i = 0; i < feat_vec.size(); ++i)
      {
        for (Size j = i+1; j < feat_vec.size(); ++j)
        {
          auto key = std::make_pair(feat_vec[i], feat_vec[j]); // feat_vec[i] is always smaller than feat_vec[j]
          auto iter = edge_map.find(key);
          if (iter == edge_map.end())
          {
            edge_map.insert(std::make_pair(key, 0));
            iter = edge_map.find(key);
          }
          ++(iter->second);
        }
      }

    }

    OpenMS::Internal::GraphMLHandler gm_handler(hypo_nodes, edge_map, out_path);
    ofstream out;
    out.open(out_path, ios::out);
    gm_handler.writeTo(out);
    out.close();

  }

  void getFeatureMassDistributions(std::map<double, FeatureFindingIntact::DeconvMassStruct>& deconv_masses,
                                   vector<FeatureFindingIntact::FeatureHypothesis>& hypotheses)
  {
    // header
    OPENMS_LOG_INFO << "feature_mass\t#features\t#charges\tcharges\tmasses\trts" << endl;

    for (auto& node : deconv_masses)
    {
      std::string charges_str = "";
      std::string masses_str = "";
      std::string rts_str = "";

      for (auto& f : node.second.feature_idx)
      {
        auto& tmp_f = hypotheses[f];
        charges_str += to_string(tmp_f.getCharge()) + ", ";
        masses_str += to_string(tmp_f.getFeatureMass()) + ", ";
        rts_str += "[" + to_string(tmp_f.getRTRange().first) + ", " + to_string(tmp_f.getRTRange().second) + "], ";
      }
      charges_str.erase(charges_str.length()-2);
      masses_str.erase(masses_str.length()-2);
      rts_str.erase(rts_str.length()-2);

      OPENMS_LOG_INFO << to_string(node.first) << "\t" << node.second.feature_masses.size() << "\t" << node.second.charges.size() << "\t"
                      << charges_str << "\t" << masses_str << "\t" << rts_str << "\n";
    }
  }

  void FeatureFindingIntact::run(std::vector<MassTrace> &input_mtraces, FeatureMap &output_featmap)
  {
    // *********************************************************** //
    // Step 1 building FeatureHypotheses
    // *********************************************************** //
    // output feature hypotheses
    std::vector<FeatureHypothesis> feat_hypos;

    // for tracking shared masstraces between features with mass traces indices
    // link masstrace_index -> hypothesis index
    std::vector<std::vector<Size>> shared_m_traces_indices(input_mtraces.size(), std::vector<Size>());

    // prepare isotope model for scoring method
    setAveragineModel_();
    if (Constants::C13C12_MASSDIFF_U * max_nr_traces_ / charge_lower_bound_ < local_mz_range_ )
    {
      local_mz_range_ = Constants::C13C12_MASSDIFF_U * max_nr_traces_ / charge_lower_bound_;
    }

    buildFeatureHypotheses_(input_mtraces, feat_hypos, shared_m_traces_indices);
//    drawSharedMasstracesBetweenFeatures(feat_hypos, shared_m_traces_indices);

    // *********************************************************** //
    // Step 2 clustering FeatureHypotheses
    // *********************************************************** //
//    getCoordinatesOfFeaturesForPython(feat_hypos);
    clusterFeatureHypotheses_(feat_hypos, shared_m_traces_indices);

    // *********************************************************** //
    // Step 3 resolving conflicts
    // *********************************************************** //
  }

  void FeatureFindingIntact::buildFeatureHypotheses_(std::vector<MassTrace>& input_mtraces,
                                                     std::vector<FeatureHypothesis>& output_hypotheses,
                                                     std::vector<std::vector<Size>>& shared_m_traces_indices) const
  {
    if (input_mtraces.empty())
    {
      return;
    }
    this->startProgress(0, input_mtraces.size(), "assembling mass traces to candidate features");

    // *********************************************************** //
    // Step 1 Preparation
    // *********************************************************** //
    // mass traces must be sorted by their centroid MZ
    std::sort(input_mtraces.begin(), input_mtraces.end(), CmpMassTraceByMZ());

    // total_intensity is needed for calculating feature hypothesis scores
    double total_intensity(0.0);
    for (Size i = 0; i < input_mtraces.size(); ++i)
    {
      total_intensity += input_mtraces[i].getIntensity(use_smoothed_intensities_);
    }

    // temporary vector
    std::vector<FeatureHypothesis> candidate_hypotheses;
    candidate_hypotheses.reserve( input_mtraces.size() * (charge_upper_bound_-charge_lower_bound_+1) );

    // vector to connect related features together
    std::map<double, DeconvMassStruct> deconv_masses; // key is the calculated median mass from the features

    // *********************************************************** //
    // Step 2 Iterate through all mass traces to find likely matches
    // and generate isotopic / charge hypotheses
    // *********************************************************** //
    Size progress(0);
    for (SignedSize i = 0; i < (SignedSize)input_mtraces.size(); ++i)
    {
        this->setProgress(progress);

      ++progress;

      std::vector<std::pair<const MassTrace*, Size>> local_traces;
      double ref_trace_mz(input_mtraces[i].getCentroidMZ());
      double ref_trace_rt(input_mtraces[i].getCentroidRT());

      local_traces.push_back(std::make_pair(&input_mtraces[i], i));

      for (Size ext_idx = i + 1; ext_idx < input_mtraces.size(); ++ext_idx)
      {
        // traces are sorted by m/z, so we can break when we leave the allowed window
        double diff_mz = std::fabs(input_mtraces[ext_idx].getCentroidMZ() - ref_trace_mz);
        if (diff_mz > local_mz_range_) break;

        double diff_rt = std::fabs(input_mtraces[ext_idx].getCentroidRT() - ref_trace_rt);
        if (diff_rt <= local_rt_range_)
        {
          // std::cout << " accepted!" << std::endl;
          local_traces.push_back(std::make_pair(&input_mtraces[ext_idx], ext_idx));
        }
      }
      findLocalFeatures_(local_traces, total_intensity, candidate_hypotheses, deconv_masses);
    }
    this->endProgress();
    candidate_hypotheses.shrink_to_fit();

    // *********************************************************** //
    // Step 3 filter out mass artifacts
    // *********************************************************** //
    // set sorted vector of deconv_masses for mass artifact removal
//    std::vector<std::pair<double, DeconvMassStruct>> deconv_mass_vec(deconv_masses.begin(), deconv_masses.end());
//    std::sort(deconv_mass_vec.begin(), deconv_mass_vec.end(), greater_equal<DeconvMassStruct>());

    OPENMS_LOG_INFO << "feature hypotheses size:" << candidate_hypotheses.size() << std::endl;
    OPENMS_LOG_INFO << "feature masses size : " << deconv_masses.size() << std::endl;

    // reduce deconv_mass
    for( auto curr_it = deconv_masses.begin(); curr_it != deconv_masses.end();  )
    {
      filterDeconvMassStruct(deconv_masses, candidate_hypotheses, curr_it);
    }
    OPENMS_LOG_INFO << "feature masses size : " << deconv_masses.size() << std::endl;
    removeMassArtifacts_(candidate_hypotheses, deconv_masses);
    OPENMS_LOG_INFO << "mass artifacts removed!" << std::endl;

    OPENMS_LOG_INFO << "feature masses size : " << deconv_masses.size() << std::endl;

    for (auto& m : deconv_masses)
    {
      if (m.second.charges.size()<2) continue;

      OPENMS_LOG_INFO << to_string(m.first) << "\t" << m.second.feature_idx.size() << "\t" << to_string(m.second.combined_score);

      std::string charges= "";
      for (auto& c : m.second.charges)
      {
        charges += to_string(c) + ",";
      }
      charges.pop_back();

      OPENMS_LOG_INFO << "\t" << charges << "\t" << to_string(m.second.fwhm_border.first) << "\t" << to_string(m.second.fwhm_border.second) << "\n";
    }

    return;

//    getFeatureMassDistributions(deconv_masses, candidate_hypotheses);

    // *********************************************************** //
    // Step 4 connect hypothesis and corresponding mass traces - for clustering afterwards
    // *********************************************************** //
    // TODO : move this section to clustering?
    // 1. set shared_m_traces_indices based on filtered hypothesis (for later, in clustering)
    const Size candi_size = candidate_hypotheses.size();
    for (Size h_index=0; h_index< candi_size; ++h_index)
    {
      const auto& current_hypo = candidate_hypotheses[h_index];
      for (const auto& mt_index : current_hypo.getIndicesOfMassTraces() )
      {
        shared_m_traces_indices[mt_index.second].push_back(h_index);
      }
    }

    output_hypotheses = candidate_hypotheses;
  }

  void FeatureFindingIntact::findLocalFeatures_(const std::vector<std::pair<const MassTrace*, Size>>& candidates,
                                                const double total_intensity,
                                                std::vector<FeatureHypothesis>& output_hypotheses,
                                                std::map<double, DeconvMassStruct>& deconv_masses) const
  {
    // not storing hypothesis with only one mass trace (with only mono), while FeatureFindingMetabo does

    // pre-save all features using multiple charge states
//    std::vector<FeatureHypothesis> hypo_for_cur_candi;

    for (int charge = charge_lower_bound_; charge <= charge_upper_bound_; ++charge)
    {
      FeatureHypothesis fh_tmp;
      fh_tmp.addMassTrace(*candidates[0].first); // ref_mtrace (which is mono here)
      double mono_mt_score = (candidates[0].first->getIntensity(use_smoothed_intensities_)) / total_intensity;
      fh_tmp.setScore(mono_mt_score);
      fh_tmp.setRTScore(0.0);
      fh_tmp.setMZScore(0.0);
      fh_tmp.setIntyScore(0.0);
      fh_tmp.addMassTraceScore(mono_mt_score);
      fh_tmp.setCharge(charge);
      double mol_weight = (candidates[0].first->getCentroidMZ()-Constants::PROTON_MASS_U) * charge;
      if (mol_weight > mass_upper_bound_)
        break;
      fh_tmp.setFeatureMass(mol_weight);

      // calculate averagine isotope distribution here (based on mono trace)
      auto iso_dist = iso_model_.get(mol_weight);
      Size iso_size = iso_dist.size();
      double iso_norm = iso_model_.getNorm(mol_weight); // l2 norm
      int iso_index_in_avg_model = 0;

      // expected m/z window for iso_pos -> 13C isotope peak position
      double mz_window = Constants::C13C12_MASSDIFF_U * max_nr_traces_ / charge;

      // pre-save all features out of current charge state
      std::vector<FeatureHypothesis> hypo_for_curr_charge;
      hypo_for_curr_charge.reserve(iso_size);

      // for shared mass traces tracker
      std::vector<std::pair<Size, Size>> used_mass_trace_indices; // first: iso index of current feature, second: mt index
      used_mass_trace_indices.reserve(iso_size);
      used_mass_trace_indices.push_back(std::make_pair(0, candidates[0].second));

      std::vector<double> feature_iso_intensities; // collection of selected isos (if nothing's in such iso index, intensity of it is 0)
      feature_iso_intensities.reserve(iso_size);
      feature_iso_intensities.push_back(candidates[0].first->getIntensity(use_smoothed_intensities_));

      Size last_iso_idx(0); // largest index of found iso index
      for (Size iso_pos = 1; iso_pos < iso_size; ++iso_pos)
      {
        // Find mass trace that best agrees with current hypothesis of charge & isotopic position
        double best_so_far(0.0);
        Size best_idx(0);
        double best_mz_score(0.0);
        double best_rt_score(0.0);
        double best_inty_score(0.0);
        for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
        {
          // if out of mz_window, pass this mass trace
          if(std::fabs(candidates[0].first->getCentroidMZ() - candidates[mt_idx].first->getCentroidMZ()) >= mz_window)
            break; // candidates is sorted. So if reached here, next candidate definitely exceed mz_window

          // Score current mass trace candidates against hypothesis
          double rt_score(scoreRT_(*candidates[0].first, *candidates[mt_idx].first));
          double mz_score(scoreMZ_(*candidates[0].first, *candidates[mt_idx].first, iso_pos, charge));

          if (rt_score <= 0.0 | mz_score <= 0)
          {
            continue;
          }

          double int_score(1.0);
          int offset; // getting isotope index against iso_dist
          std::vector<double> tmp_ints(feature_iso_intensities); // intensities up to the last isotope
          tmp_ints.push_back(candidates[mt_idx].first->getIntensity(use_smoothed_intensities_));
          int_score = computeAveragineCosineSimScore_(tmp_ints, iso_dist, iso_size, iso_norm, offset);

          double total_pair_score(0.0);
          if (int_score > 0.0)
          {
            total_pair_score = std::exp(std::log(rt_score) + log(mz_score) + log(int_score));
          }
          if (total_pair_score > best_so_far)
          {
//            iso_index_in_avg_model = -offset;
            best_so_far = total_pair_score;
            best_idx = mt_idx;
            best_rt_score = rt_score;
            best_mz_score = mz_score;
            best_inty_score = int_score;
          }
        } // end of mt_idx

        // Store mass trace that best agrees with current hypothesis of charge and isotopic position
        if (best_so_far > 0.0)
        {
          double weighted_score(((candidates[best_idx].first->getIntensity(use_smoothed_intensities_)) * best_so_far) / total_intensity);
          fh_tmp.setScore(fh_tmp.getScore() + weighted_score);
          fh_tmp.addMassTraceScore(weighted_score);
          fh_tmp.addMassTrace(*candidates[best_idx].first);
          fh_tmp.setMZScore(fh_tmp.getMZScore() + best_mz_score);
          fh_tmp.setRTScore(fh_tmp.getRTScore() + best_rt_score);
          fh_tmp.setIntyScore(fh_tmp.getIntyScore() + best_inty_score);

          // save up information for trackers
          used_mass_trace_indices.push_back(std::make_pair(iso_pos, candidates[best_idx].second));
          feature_iso_intensities.push_back(candidates[best_idx].first->getIntensity(use_smoothed_intensities_));

          // modify isotope index for this feature
          auto tmp_indices(used_mass_trace_indices);
          for (Size i = 0; i< used_mass_trace_indices.size(); i++)
          {
            tmp_indices[i] = std::make_pair(tmp_indices[i].first + iso_index_in_avg_model, tmp_indices[i].second);
          }
          fh_tmp.setIndicesOfMassTraces(tmp_indices);

          // set rt and mz range of current hypothesis
          fh_tmp.setRTRange();
          fh_tmp.setMZRange(mz_window);
          fh_tmp.setFwhmRange();

          // update feature mass with modified isotope index
          if (iso_index_in_avg_model != 0)
          {
            fh_tmp.updateFeatureMass();
          }

//          output_hypotheses.push_back(fh_tmp);
          hypo_for_curr_charge.push_back(fh_tmp);
          last_iso_idx = best_idx;
        }
        else
        {
          feature_iso_intensities.push_back(0.);
          // if too many missing iso positions, break; TODO: find an elegant way...
          if (feature_iso_intensities.size() > 3 && feature_iso_intensities.rbegin()[0] == 0
              && feature_iso_intensities.rbegin()[1] == 0 && feature_iso_intensities.rbegin()[2] == 0)
          {
            break;
          }
        }
      } // end for iso_pos
      // using only one hypothesis per charge
      if (hypo_for_curr_charge.size() == 0)
      {
        continue;
      }
      double max_score = 0;
      Size max_score_index = 0;
      for (Size i = 0; i<hypo_for_curr_charge.size(); ++i)
      {
        if (hypo_for_curr_charge[i].getScore() > max_score)
        {
          max_score = hypo_for_curr_charge[i].getScore();
          max_score_index = i;
        }
      }

//      hypo_for_cur_candi.push_back(hypo_for_curr_charge[max_score_index]);
      output_hypotheses.push_back(hypo_for_curr_charge[max_score_index]);
      addFeature2DeconvMassStruct(output_hypotheses.back(), output_hypotheses.size()-1, deconv_masses);

      hypo_for_curr_charge.empty();
    } // end of charge

//    if (hypo_for_cur_candi.empty()) return;
    // to output the hypothesis with the highest score
//    std::sort(hypo_for_cur_candi.begin(), hypo_for_cur_candi.end(), greater<FeatureHypothesis>());
//    for (auto& h : hypo_for_cur_candi){
//      output_hypotheses.push_back(h);
//      addFeature2DeconvMassStruct(output_hypotheses.back(), output_hypotheses.size()-1, deconv_masses);
//    }

  }  // end of findLocalFeatures_(...)

  void FeatureFindingIntact::addFeature2DeconvMassStruct(FeatureHypothesis &in_feature,
                                                         Size feature_idx,
                                                         std::map<double, DeconvMassStruct> &deconv_masses) const
  {// add to deconv_mass_map
    map<double, DeconvMassStruct>::const_iterator low_it;
    map<double, DeconvMassStruct>::const_iterator up_it;

    const double &current_mass = in_feature.getFeatureMass();
    //    const double &tol = current_mass * mass_tolerance_ * 1e-6;
    low_it = deconv_masses.lower_bound(current_mass - mass_tolerance_); // less than or equal to
    up_it = deconv_masses.upper_bound(current_mass + mass_tolerance_); // greater than

    // no matching in deconv_masses map
    if (low_it == up_it)
    {
      DeconvMassStruct dms;
      dms.initialize(current_mass, in_feature.getCharge(), feature_idx, in_feature.getFwhmRange(), in_feature.getScore());
      deconv_masses.insert(make_pair(current_mass, dms));
      return;
    }

    // find nearest deconv mass
    double smallest_diff(mass_tolerance_);
    auto selected_it = low_it;
    bool not_selected = true;
    for (; low_it != up_it; ++low_it)
    {
      double diff = abs(low_it->first-current_mass);
      if ( diff < smallest_diff && doFWHMbordersOverlap(in_feature.getFwhmRange(), low_it->second.fwhm_border) )
      {
        selected_it = low_it;
        smallest_diff = diff;
        not_selected=false;
      }
    }

    // if no deconv_mass is within tolerance and overlapping FWHM border
    if(not_selected)
    {
      DeconvMassStruct dms;
      dms.initialize(current_mass, in_feature.getCharge(), feature_idx, in_feature.getFwhmRange(), in_feature.getScore());
      deconv_masses.insert(make_pair(current_mass, dms));
      return;
    }

    // add current feature to deconv_massses map
    double key_mass = selected_it->first;
    deconv_masses[key_mass].addFeatureHypothesis(current_mass, in_feature.getCharge(), feature_idx,
                                                 in_feature.getFwhmRange(), in_feature.getScore());

    // update key of map if needed
    if (deconv_masses[key_mass].updateDeconvMass())
    {
      auto curr_node = deconv_masses.extract(key_mass);
      curr_node.key() = curr_node.mapped().deconv_mass;
      deconv_masses.insert(move(curr_node));
    }
  }

  double FeatureFindingIntact::scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const
  {
    std::map<double, std::pair<double, double> > coinciding_rts;

    std::pair<Size, Size> tr1_fwhm_idx(tr1.getFWHMborders());
    std::pair<Size, Size> tr2_fwhm_idx(tr2.getFWHMborders());

    int coinciding_counter = 0;

    // Extract peak shape between FWHM borders for both peaks
    for (Size i = tr1_fwhm_idx.first; i <= tr1_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr1[i].getRT()] = std::make_pair(tr1[i].getIntensity(), 0.0);
    }
    for (Size i = tr2_fwhm_idx.first; i <= tr2_fwhm_idx.second; ++i)
    {
      auto it = coinciding_rts.find(tr2[i].getRT());
      if (it != coinciding_rts.end())
      {
        it->second.second = tr2[i].getIntensity();
        coinciding_counter++;
      }
      else
      {
        coinciding_rts[tr2[i].getRT()] = std::make_pair(0.0, tr2[i].getIntensity());
      }
    }

    if (coinciding_counter==0)
    {
      return 0;
    }

    std::vector<double> x, y;
    for (auto& rts : coinciding_rts)
    {
      x.push_back(rts.second.first);
      y.push_back(rts.second.second);
    }
    return computeCosineSim_(x, y);
  }

  double FeatureFindingIntact::scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, int charge) const
  {
    double diff_mz(std::fabs(tr2.getCentroidMZ() - tr1.getCentroidMZ()));

    double mt_sigma1(tr1.getCentroidSD());
    double mt_sigma2(tr2.getCentroidSD());
    double mt_variances(std::exp(2 * std::log(mt_sigma1)) + std::exp(2 * std::log(mt_sigma2)));

    double mz_score(0.0);
    /// mz scoring by expected mean w/ C13
    double mu = (Constants::C13C12_MASSDIFF_U * iso_pos) / charge; // using '1.0033548378'
//    double sd = (0.0016633 * iso_pos - 0.0004751) / charge;
    double sigma_mult(3.0);

    //standard deviation including the estimated isotope deviation
//    double score_sigma(std::sqrt(std::exp(2 * std::log(sd)) + mt_variances));
    double score_sigma(std::sqrt(mt_variances));

    if ((diff_mz < mu + sigma_mult * score_sigma) && (diff_mz > mu - sigma_mult * score_sigma))
    {
      double tmp_exponent((diff_mz - mu) / score_sigma);
      mz_score = std::exp(-0.5 * tmp_exponent * tmp_exponent);
    }

    return mz_score;
  }

  bool FeatureFindingIntact::doFWHMbordersOverlap(const std::pair<double, double>& border1,
                                                  const std::pair<double, double>& border2) const
  {
    if ( (border1.first > border2.second) || (border2.first > border1.second))
      return false;

    const double overlap_length = std::min(border1.second, border2.second)-std::max(border1.first, border2.first);
    if ( (overlap_length/(border1.second-border1.first) < 0.7) &&
         (overlap_length/(border2.second-border2.first) < 0.7) ) return false;

    return true;
  }

  // from FLASHDeconvAlgorithm::getCosine
  double FeatureFindingIntact::computeCosineSimOfDiffSizedVector_(const std::vector<double>& a,
                                         const IsotopeDistribution& b,
                                         const int& b_size,
                                         const double& b_norm,
                                         const int offset) const
  {
    double n = .0, a_norm= .0;
    for (Size j = 0; j < a.size(); j++)
    {
      a_norm += a[j] * a[j];
      int i = j - offset;
      if (i < 0 || i >= b_size)
      {
        continue;
      }
      n += a[j] * b[i].getIntensity(); //
    }
    double d = (a_norm * b_norm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  // modified based on FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex
  // isotopic intensity distribution comparison using the cosine similarity
  double FeatureFindingIntact::computeAveragineCosineSimScore_(const std::vector<double>& hypo_ints,
                                                               const IsotopeDistribution& iso_dist,
                                                               const Size& iso_size,
                                                               const double& iso_norm,
                                                               int& offset) const
  {
    // Unlike FLASHDeconv, candidate should be subset of averagine.
    // iso_size is always larger than hypo_size
    offset = 0; // initialize
    double max_cosine = -1;
    int hypo_size = hypo_ints.size();

    // find largest intensity in hypo_ints for normalization
    double max_inty_of_hypo = 0;
    for (int i = 0; i < hypo_size; ++i)
    {
      if (hypo_ints[i] > max_inty_of_hypo)
      {
        max_inty_of_hypo = hypo_ints[i];
      }
    }

    // normalize hypo_ints
    std::vector<double> hypo_norm_ints;
    for (int i = 0; i < hypo_size; ++i)
    {
      hypo_norm_ints.push_back(hypo_ints[i] / max_inty_of_hypo);
    }

    // determining isotope index within averagine guaranteeing highest score
    int max_cntr = 0;
    vector<double> tmp_cos_vec;
    for (int tmp_offset = -(iso_size-hypo_size); tmp_offset <= 0; tmp_offset++)
    {
      double tmp_cos = computeCosineSimOfDiffSizedVector_(hypo_norm_ints,
                                                          iso_dist,
                                                          iso_size,
                                                          iso_norm,
                                                          tmp_offset);
      tmp_cos_vec.push_back(tmp_cos);

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

  double FeatureFindingIntact::computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const
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

  void FeatureFindingIntact::setAveragineModel_()
  {
    auto generator = new CoarseIsotopePatternGenerator();
    auto maxIso = generator->estimateFromPeptideWeight(mass_upper_bound_);
    maxIso.trimRight(0.01 * maxIso.getMostAbundant().getIntensity());

    generator->setMaxIsotope(maxIso.size());
    iso_model_ = PrecalculatedAveragine(50, mass_upper_bound_, 25, generator);
    iso_model_.setMaxIsotopeIndex(maxIso.size() - 1);

    max_nr_traces_ = iso_model_.getMaxIsotopeIndex();
  }

  void FeatureFindingIntact::removeMassArtifacts_(const std::vector<FeatureHypothesis>& feat_hypotheses,
                                                  std::map<double, DeconvMassStruct>& deconv_masses) const
  {
    this->startProgress(0, deconv_masses.size(), "detecting mass artifacts");
    Size progress(0);

    double m_tol = mass_tolerance_;
    for( auto curr_it = deconv_masses.begin(); curr_it != deconv_masses.end();  )
    {
      this->setProgress(progress);
      progress++;

      bool is_artifact = false;
      bool is_current_struct_changed = false;
      // TODO : try 20 ppm tolerance?
//      double m_tol = 1e-5 * cur_fh.getFeatureMass(); // 10ppm

      // 1. check if harmonic artifacts
      for (int iso_off = -2; iso_off <= 2 && !is_artifact; ++iso_off) // up to 2 isotope error
      {
        double curr_mass =  curr_it->first + iso_off*Constants::C13C12_MASSDIFF_U;
        for (int harmonic = 2; harmonic <= 6 && !is_artifact; ++harmonic) // harmonics : 2, 4, 5, 6
        {
          const std::vector<double> h_masses = {curr_mass / harmonic, curr_mass * harmonic}; // low harmonic, high harmonic

          for (auto& h_mass : h_masses)
          {
            std::map<double, DeconvMassStruct>::const_iterator low_it;
            std::map<double, DeconvMassStruct>::const_iterator up_it;

            low_it = deconv_masses.lower_bound(h_mass - m_tol); // less than or equal to
            up_it = deconv_masses.upper_bound(h_mass + m_tol); // greater than

            // harmonic mass is found. check if score is higher & fwhms overlap
            for (; low_it != up_it; ++low_it)
            {
              if ( (low_it->second.combined_score > curr_it->second.combined_score) &&
                  doFWHMbordersOverlap(low_it->second.fwhm_border, curr_it->second.fwhm_border) )
              {
                is_artifact = true;
                break;
              }
            }
          }
        }
      }

      // 2. check if charge-off-by-n artifacts
      double current_mass = curr_it->first;
      for (auto cs_itr = curr_it->second.charges.begin(); cs_itr != curr_it->second.charges.end() && !is_artifact; )
      {
        int curr_cs = *cs_itr;
        bool is_curr_cs_removed = false;

        for ( int coff = 1; coff <= 4 && !is_artifact && !is_curr_cs_removed; ++coff)
        {
          for (int f = -1; f <= 1 && !is_artifact && !is_curr_cs_removed; f += 2) // for converting minus / plus
          {
            int tmp_cs = curr_cs + f * coff;
            if (tmp_cs <= 1)
            {
              continue;
            }

            // if this artifact is possible (distance between peaks of two consecutive charges are smaller than 4 Th)
//            if (std::fabs((current_mass/curr_cs)-(current_mass/tmp_cs)) > 4 ) continue;

            // check if harmonic mass exist in the map
            double h_mass = current_mass / curr_cs * tmp_cs;

            std::map<double, DeconvMassStruct>::const_iterator low_it;
            std::map<double, DeconvMassStruct>::const_iterator up_it;

            low_it = deconv_masses.lower_bound(h_mass - m_tol); // less than or equal to
            up_it = deconv_masses.upper_bound(h_mass + m_tol); // greater than

            // harmonic mass is found. check if score is higher & fwhms overlap
            for (; low_it != up_it; ++low_it)
            {
              if ( (low_it->second.combined_score > curr_it->second.combined_score) &&
                   doFWHMbordersOverlap(low_it->second.fwhm_border, curr_it->second.fwhm_border) )
              {
                // if current deconv_mass contains only one feature, remove it from deconv_mass
                if (curr_it->second.feature_idx.size()==1)
                {
                  is_artifact=true;
                  break;
                }

                // remove current feature hypothesis from deconv_mass map
                for (auto& f_index : curr_it->second.feature_idx)
                {
                  // all features having current charge should be removed. (all of them are harmonics)
                  if (feat_hypotheses[f_index].getCharge() == curr_cs)
                  {
                    curr_it->second.removeFeatureHypothesis(feat_hypotheses[f_index].getFeatureMass(),
                                                            feat_hypotheses[f_index].getScore(), f_index);
                  }
                }
                // remove current charge from deconv_mass map
                is_curr_cs_removed = true;
                is_current_struct_changed = true;

                // if most of features in current mass struct are removed, remove mass struct from map.
                if (curr_it->second.feature_idx.size() < 2)
                {
                  is_artifact=true;
                  break;
                }
              }
            }
          }
        }
        if (is_curr_cs_removed)
        {
          cs_itr = curr_it->second.charges.erase(cs_itr);
        }
        else
        {
          ++cs_itr;
        }
      }
      // if most of features in current deconv_mass are removed, it's artifact.

      // if artifact, remove it from the map
      if (is_artifact)
      {
        curr_it = deconv_masses.erase(curr_it);
      }
      else
      {
        // TODO; update current mass (if needed - should this have flag too?)
        filterDeconvMassStruct(deconv_masses, feat_hypotheses, curr_it, is_current_struct_changed);
      }
    }
    this->endProgress();


    //    /// remove charge artifacts
    //    for (Size x=0; x<hypo_for_cur_candi.size(); x++)
    //    {
    //      auto& cur_candi = hypo_for_cur_candi[x];
    //      bool is_artifact = false;
    //      for (Size y=x+1; y<hypo_for_cur_candi.size(); y++)
    //      {
    //        double m_tol = 1e-5 * cur_candi.getFeatureMass(); // 10ppm
    //        auto& comp_candi = hypo_for_cur_candi[y];
    //
    //        // 1. check if harmonic artifacts
    //        for (const Size& cs : harmonic_charges_)
    //        {
    //          for (int i = -10; i <= 10; ++i) // up to 10 isotope error
    //          {
    //            double l_mass = (i*Constants::C13C12_MASSDIFF_U + cur_candi.getFeatureMass()) * cs; // low harmonic
    //            double h_mass = (i*Constants::C13C12_MASSDIFF_U + cur_candi.getFeatureMass()) / cs; // high harmonic
    //            if (std::fabs(l_mass-comp_candi.getFeatureMass()) < m_tol ||
    //                std::fabs(h_mass-comp_candi.getFeatureMass()) < m_tol )
    //            {
    //              is_artifact = true;
    //              break;
    //            }
    //          }
    //          if (is_artifact) break;
    //        }
    //        if (is_artifact) break;
    //
    //        // 2. check if charge-off-by-n artifacts
    //        Size curr_charge = cur_fh.getCharge();
    //        if (std::abs(int(curr_charge)-int(overlap_iter->getCharge())) != 1) continue;
    //        // if this artifact is possible (distance between peaks of two consecutive charges are smaller than 4 Th)
    //        if ((cur_fh.getFeatureMass()/curr_charge)-(cur_fh.getFeatureMass()/(curr_charge+1)) > 4 ) continue;
    //        for (int i = -10; i <= 10; ++i) // up to 10 isotope error
    //        {
    //          double mass_cs_up = (i*Constants::C13C12_MASSDIFF_U + cur_fh.getFeatureMass())/curr_charge * (curr_charge+1);
    //          double mass_cs_down = (i*Constants::C13C12_MASSDIFF_U + cur_fh.getFeatureMass())/curr_charge * (curr_charge-1);
    //          if (std::fabs(mass_cs_up-overlap_iter->getFeatureMass()) < m_tol ||
    //              std::fabs(mass_cs_down-overlap_iter->getFeatureMass()) < m_tol )
    //          {
    //            is_artifact = true;
    //            break;
    //          }
    //        }
    //        if (is_artifact) break;
    //      }
    //
    //      // if not artifacts, add it to output hypothesis
    //      if (!is_artifact)
    //      {
    //        output_hypotheses.push_back(cur_candi);
    //      }
    //    }


  }

  void FeatureFindingIntact::setChargeScoreForFeatureHypothesis(std::vector<FeatureHypothesis>& candidate_hypotheses,
                                                          std::vector<std::pair<double, int>>& feat_and_charges) const
  {
    Size candi_size = candidate_hypotheses.size();

    // calculate mean_masses within tolerance
    std::vector<std::pair<double, std::set<int>>> collected_masses; // first : average of feature masses, second : vector of charges
    collected_masses.reserve(candi_size);
    const double tolerance = 1.5; // Da
    Size max_nr_charges = 0; // save the maximum size of charge set
    auto&& end_of_feature = feat_and_charges.end();

    this->startProgress(0, feat_and_charges.size(), "charge scoring");
    Size progress = 0;
    std::vector<std::pair<double, int>>::iterator f_iter=feat_and_charges.begin();
    while( f_iter != end_of_feature )
    {
      setProgress(progress);
      progress++;

      // TODO : just to prevent the error - couldn't figure out why tho...
      if(std::distance(f_iter, end_of_feature) < 0) break;

      auto f_iter_end = f_iter;

      double average = 0.0;
      int counter = 1;
      for(; f_iter_end!=end_of_feature; ++f_iter_end)
      {
        average = 0.0;
        auto next_end = std::next(f_iter_end);
        for (auto i=f_iter; i!=next_end; ++i) average += i->first;
        average /= counter;
        if (next_end->first-average > tolerance)
        {
          break;
        }
        ++counter;
      }

      // save found charge states
      std::set<int> tmp_charges;
      f_iter_end++; // starting position of next f_iter
      for(auto it=f_iter; it!=f_iter_end; ++it)
      {
        tmp_charges.insert(it->second);
      }
      collected_masses.push_back(std::make_pair(average, tmp_charges));
      if (tmp_charges.size() > max_nr_charges)
      {
        max_nr_charges = tmp_charges.size();
      }

      f_iter = std::next(f_iter, std::distance(f_iter, f_iter_end));
      //      OPENMS_LOG_INFO << std::distance(feat_and_charges.end(), f_iter) << endl;
    }
    collected_masses.shrink_to_fit();
    this->endProgress();

    // give charge score for each feature
    for (auto& feature : candidate_hypotheses)
    {
      double cs_score=0;
      for (auto f_it = collected_masses.begin(); f_it != collected_masses.end(); ++f_it)
      {
        if (std::fabs(f_it->first-feature.getFeatureMass()) > tolerance)
          continue;
        cs_score = f_it->second.size() / max_nr_charges;
        break;
      }
      feature.setChargeScore(cs_score);
    }
  }

  void FeatureFindingIntact::clusterFeatureHypotheses_(vector<FeatureHypothesis>& hypotheses,
                                                       const std::vector<std::vector<Size>>& shared_m_traces_indices) const
  {

    // test writing
    String out_path = "/Users/jeek/Documents/A4B_UKE/filg/rp/centroid/190904_Proswift50mm_7min_Filg_500ng.pphr.featHypoClusterNameHighest.tsv";
    ofstream out;
    out.open(out_path, ios::out);

    // header
    out << "feature_label\tcs\tscore\tquant\tbounding_box_pos\tbounding_box_width\tbounding_box_height\t"
           "iso_position\tmasstrace_centroid_rts\tmasstrace_centroid_mzs\tclustername\n";

    // *********************************************************** //
    // Step 1 constructing hypergraph from featurehypotheses
    //        node = mass traces
    //        hyperedge = hypotheses
    // *********************************************************** //
    Size num_nodes = shared_m_traces_indices.size();
    std::vector<bool> bfs_visited;
    bfs_visited.resize(num_nodes, false);
    std::queue<Size> bfs_queue;
    Size search_pos = 0; // keeping track of mass trace index to look for seed

    std::vector<FeatureHypothesis> out_features;

    // BFS
    // TODO : progress logger
    this->startProgress(0, shared_m_traces_indices.size(), "clustering features based on the shared mass traces");
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
        for (vector<Size>::const_iterator it = shared_m_traces_indices[i].begin();
             it != shared_m_traces_indices[i].end();
             ++it)
        {
          hypo_indices_in_current_cluster.insert(*it);
          FeatureHypothesis &current_hypo = hypotheses[*it];

          for (const auto &iso_mt_pair : current_hypo.getIndicesOfMassTraces())
          {
            Size j = iso_mt_pair.second;
            if (!bfs_visited[j])
            {
              bfs_queue.push(j);
              bfs_visited[j] = true;
            }
          }
        }
      }
      // current cluster
      if (hypo_indices_in_current_cluster.empty()) continue; // no cluster out of current mt
      cluster_counter++;
      if (hypo_indices_in_current_cluster.size() == 1){
//         no conflict, but cannot happen.
        out_features.push_back(hypotheses[*(hypo_indices_in_current_cluster.begin())]);
        continue;
      }
      String cluster_name = "cluster " + to_string(progress);
      resolveConflictInCluster_(hypotheses, shared_m_traces_indices, hypo_indices_in_current_cluster, out_features, out, cluster_name);
    }
    this->endProgress();

    out.close();
    OPENMS_LOG_INFO << "#cluster :" << cluster_counter << endl;
  }

  void getCoordinatesOfClusterForPython(std::vector<FeatureFindingIntact::FeatureHypothesis> hypo, ofstream& out, String& cluster_name)
  {
    for (const auto& feat : hypo)
    {
      double mz_upper_limit = 0.0;
      double rt_upper_limit  = 0.0;
      double mz_lower_limit = std::numeric_limits<double>::max();
      double rt_lower_limit = std::numeric_limits<double>::max();
      stringstream rts;
      stringstream mzs;
      double quant(0.0);
      for (const auto& mt : feat.getMassTraces())
      {
        const auto& boundingbox = mt->getConvexhull().getBoundingBox();

        if (boundingbox.minX() < rt_lower_limit)
          rt_lower_limit = boundingbox.minX();
        if (boundingbox.maxX() > rt_upper_limit)
          rt_upper_limit = boundingbox.maxX();
        if (boundingbox.minY() < mz_lower_limit)
          mz_lower_limit = boundingbox.minY();
        if (boundingbox.maxY() > mz_upper_limit)
          mz_upper_limit = boundingbox.maxY();

        mzs << mt->getCentroidMZ() << ",";
        rts << mt->getCentroidRT() << ",";

        quant += mt->computePeakArea();
      }
      std::string centroids = rts.str();
      centroids.pop_back();
      centroids = centroids + "\t" + mzs.str();
      centroids.pop_back();

      stringstream isos;
      for (const auto& pair : feat.getIndicesOfMassTraces())
      {
        isos << pair.first << ",";
      }
      std::string iso_str = isos.str();
      iso_str.pop_back();

      String label = std::to_string(feat.getFeatureMass()) + "\t" + std::to_string(feat.getCharge());
      out << label << "\t" << to_string(feat.getScore()) << "\t" << std::to_string(quant) << "\t"
          << to_string(rt_lower_limit) << "," << to_string(mz_lower_limit) << "\t"
          << to_string((rt_upper_limit-rt_lower_limit)) << "\t" << to_string((mz_upper_limit-mz_lower_limit)) << "\t"
          << iso_str << "\t" << centroids
          << "\t" << cluster_name << "\n";
    }
  }

  void FeatureFindingIntact::resolveConflictInCluster_(const std::vector<FeatureHypothesis>& feat_hypo,
                                 const std::vector<std::vector<Size> >& shared_m_traces_indices,
                                 const std::set<Size>& hypo_indices,
                                 std::vector<FeatureHypothesis>& out_features,
                                 ofstream& out,
                                 String& cluster_name) const
  {
    // compute scores for each features
//    vector<double> feat_scores(hypo_indices.size());
    std::vector<FeatureHypothesis> tmp_cluster;
    for (const auto& h_index : hypo_indices)
    {
//      feat_scores.push_back( feat_hypo[h_index].getScore() );
      tmp_cluster.push_back(feat_hypo[h_index]);

    }
    std::sort(tmp_cluster.begin(), tmp_cluster.end(), greater<FeatureHypothesis>());

    // remove until no mass traces are left
    std::vector<FeatureHypothesis> output_feats;
//    for (Size i=0; i<tmp_cluster.size(); i++)
//    {
//      Size mt_size = 0;
//      std::vector<Size> mt_indices;
//      // check if this feature contains enough mass traces
//      for(auto& tmp_pair : tmp_cluster[i].getIndicesOfMassTraces())
//      {
//        bool is_removed = false;
//        for (auto& removed : removed_mt_indices)
//        {
//          if (tmp_pair.second == removed)
//          {
//            is_removed = true;
//            break;
//          }
//        }
//        if (!is_removed)
//        {
//          mt_size++;
//          mt_indices.push_back(tmp_pair.second);
//        }
//      }
//
//      // if less than 3 mass traces, ignore
//      if (mt_size < 3) continue;
//
//      // add to output, and remove used mass traces
//      output_feats.push_back(tmp_cluster[i]);
//      removed_mt_indices.insert(removed_mt_indices.end(), mt_indices.begin(), mt_indices.end());
//    }

    while (tmp_cluster.size() > 0)
    {
      output_feats.push_back(tmp_cluster[0]);

      std::vector<Size> removed_mt_indices;
      for(auto& tmp_pair : tmp_cluster[0].getIndicesOfMassTraces())
      {
        removed_mt_indices.push_back(tmp_pair.second);
      }
      tmp_cluster.erase(tmp_cluster.begin());

      auto feature = std::begin(tmp_cluster);
      while(feature != std::end(tmp_cluster))
      {
        Size feat_mt_size = feature->getIndicesOfMassTraces().size();
        std::vector<Size> indices_to_remove;
        for (Size index = 0; index < feat_mt_size; index++)
        {
          Size curr_mt_index = feature->getIndicesOfMassTraces()[index].second;
          for (auto& removed : removed_mt_indices)
          {
            if (removed==curr_mt_index)
            {
              indices_to_remove.push_back(index);
              break;
            }
          }
        }

        // removing mass traces, with indices. so backward
        for (auto it=indices_to_remove.rbegin(); it!=indices_to_remove.rend(); ++it)
        {
          feature->removeMassTrace(*it);
        }

        if (feature->getIndicesOfMassTraces().size() < 3)
        {
          feature = tmp_cluster.erase(feature);
        }
        else
        {
          feature->reset_feature_score();
          ++feature;
        }
      }
      std::sort(tmp_cluster.begin(), tmp_cluster.end(), greater<FeatureHypothesis>());
    }

    getCoordinatesOfClusterForPython(output_feats, out, cluster_name);
  }

  void FeatureFindingIntact::filterDeconvMassStruct(std::map<double, DeconvMassStruct> &deconv_masses,
                                                    const vector<FeatureHypothesis> &feat_hypotheses,
                                                    std::map<double, DeconvMassStruct>::iterator &curr_it,
                                                    const bool struct_needs_update) const
  {// if struct contains only one feature or only one charge state, remove it.
    if (curr_it->second.feature_idx.size() == 1 || curr_it->second.charges.size() == 1)
    {
      // check if similar mass exist in the map within tolerance
      auto new_it = curr_it; // reference for new location to save current mass
      auto prev_it = std::prev(curr_it);
      if( curr_it != deconv_masses.begin() && (curr_it->first - prev_it->first) <= mass_tolerance_
          && doFWHMbordersOverlap(curr_it->second.fwhm_border, prev_it->second.fwhm_border))
      {
        new_it = prev_it;
      }

      auto next_it = std::next(curr_it);
      if ( next_it != deconv_masses.end() && (next_it->first - curr_it->first) <= mass_tolerance_
           && doFWHMbordersOverlap(curr_it->second.fwhm_border, next_it->second.fwhm_border) )
      {
        if (new_it!=curr_it) // if both prev and next have similar mass to curr_it
        {
          // the new location is up to which one is closer to curr_it
          new_it = ( (curr_it->first-new_it->first) < (next_it->first-curr_it->first) ) ? new_it : next_it;
        }
        else
        {
          new_it = next_it;
        }
      }

      // add features in current mass to new location
      if (new_it != curr_it)
      {
        for (auto& fidx : curr_it->second.feature_idx)
        {
          auto& feature = feat_hypotheses[fidx];
          new_it->second.addFeatureHypothesis(feature.getFeatureMass(), feature.getCharge(), fidx,
                                              feature.getFwhmRange(), feature.getScore());
        }
        // update key of map if needed
        if (new_it->second.updateDeconvMass())
        {
          auto curr_node = deconv_masses.extract(new_it->first);
          curr_node.key() = curr_node.mapped().deconv_mass;
          deconv_masses.insert(move(curr_node));
        }
      }

      curr_it = deconv_masses.erase(curr_it);
      return;
    }
    // if DeconvMassStruct needs update
    if (struct_needs_update)
    {
      // update fwhm_border
      curr_it->second.fwhm_border = feat_hypotheses[curr_it->second.feature_idx[0]].getFwhmRange();
      for (Size i=1; i<curr_it->second.feature_idx.size(); i++)
      {
        curr_it->second.updateFwhmBorder(feat_hypotheses[curr_it->second.feature_idx[i]].getFwhmRange());
      }

      // update key of map if needed
      if (curr_it->second.updateDeconvMass())
      {
        auto curr_node = deconv_masses.extract(curr_it->first);
        curr_node.key() = curr_node.mapped().deconv_mass;
        deconv_masses.insert(move(curr_node));
      }
    }

    ++curr_it;
  }

}