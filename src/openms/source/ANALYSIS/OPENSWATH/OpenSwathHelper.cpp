// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

#include <random>
#include <algorithm>
#include <unordered_set>

namespace OpenMS
{
  void OpenSwathHelper::selectSwathTransitions(const OpenMS::TargetedExperiment& targeted_exp,
                                               OpenMS::TargetedExperiment& transition_exp_used, double min_upper_edge_dist,
                                               double lower, double upper)
  {
    transition_exp_used.setPeptides(targeted_exp.getPeptides());
    transition_exp_used.setProteins(targeted_exp.getProteins());
    for (Size i = 0; i < targeted_exp.getTransitions().size(); i++)
    {
      ReactionMonitoringTransition tr = targeted_exp.getTransitions()[i];
      if (lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < upper &&
          std::fabs(upper - tr.getPrecursorMZ()) >= min_upper_edge_dist)
      {

         OPENMS_LOG_DEBUG << "Adding Precursor with m/z " << tr.getPrecursorMZ() <<  " to swath with mz lower of " << lower << " m/z upper of " << upper;
        transition_exp_used.addTransition(tr);
      }
    }
  }

  // For PASEF experiments it is possible to have DIA windows with the same m/z however different IM.
  // Extract from the DIA window in which the precursor is more centered across its IM.
  // Unlike the function above, current implementation may not be parrelization safe
  void OpenSwathHelper::selectSwathTransitionsPasef(const OpenSwath::LightTargetedExperiment& transition_exp, std::vector<int>& tr_win_map,
                                               double min_upper_edge_dist, const std::vector< OpenSwath::SwathMap > & swath_maps)
  {
      OPENMS_PRECONDITION(std::any_of(transition_exp.transitions.begin(), transition_exp.transitions.end(), [](auto i){return i.getPrecursorIM()!=-1;}), "All transitions must have a valid IM value (not -1)");

      tr_win_map.resize(transition_exp.transitions.size(), -1);
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
      {
        for (Size k = 0; k < transition_exp.transitions.size(); k++)
        {
          const OpenSwath::LightTransition& tr = transition_exp.transitions[k];

          // If the transition falls inside the current DIA window (both in IM and m/z axis), check
          // if the window is potentially a better match for extraction than
          // the one previously stored in the map:
          if (
             swath_maps[i].imLower < tr.getPrecursorIM() && tr.getPrecursorIM() < swath_maps[i].imUpper &&
             swath_maps[i].lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < swath_maps[i].upper &&
             std::fabs(swath_maps[i].upper - tr.getPrecursorMZ()) >= min_upper_edge_dist )
          {
            if (tr_win_map[k] == -1)
            {
              tr_win_map[k] = i;
            }
            else
            {
              // Check if the current window is better than the previously assigned window (across IM)
              double imOld = std::fabs(((swath_maps[ tr_win_map[k] ].imLower + swath_maps [ tr_win_map[k] ].imUpper) / 2) - tr.getPrecursorIM() );
              double imNew = std::fabs(((swath_maps[ i ].imLower + swath_maps [ i ].imUpper) / 2) - tr.getPrecursorIM() );
              if (imOld > imNew)
              {
                // current DIA window "i" is a better match
                OPENMS_LOG_DEBUG << "For Precursor " << tr.getPrecursorIM() << " Replacing Swath Map with IM center of " <<
                  imOld << " with swath map of im center " << imNew << std::endl;
                tr_win_map[k] = i;
              }
            }
          }
        }
      }
    }

  void OpenSwathHelper::checkSwathMap(const OpenMS::PeakMap& swath_map,
                                      double& lower, double& upper, double& center)
  {
    if (swath_map.empty() || swath_map[0].getPrecursors().empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Swath map has no Spectra");
    }
    const std::vector<Precursor>& first_prec = swath_map[0].getPrecursors();
    lower = first_prec[0].getMZ() - first_prec[0].getIsolationWindowLowerOffset();
    upper = first_prec[0].getMZ() + first_prec[0].getIsolationWindowUpperOffset();
    center = first_prec[0].getMZ();
    UInt expected_mslevel = swath_map[0].getMSLevel();

    for (Size index = 0; index < swath_map.size(); index++)
    {
      const std::vector<Precursor>& prec = swath_map[index].getPrecursors();
      if (prec.size() != 1)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Scan " + String(index) + " does not have exactly one precursor.");
      }
      if (swath_map[index].getMSLevel() != expected_mslevel)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Scan " + String(index) + " if of a different MS level than the first scan.");
      }
      if (
        fabs(prec[0].getMZ() - first_prec[0].getMZ()) > 0.1 ||
        fabs(prec[0].getIsolationWindowLowerOffset() - first_prec[0].getIsolationWindowLowerOffset()) > 0.1 ||
        fabs(prec[0].getIsolationWindowUpperOffset() - first_prec[0].getIsolationWindowUpperOffset()) > 0.1
        )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Scan " + String(index) + " has a different precursor isolation window than the first scan.");
      }
    }
  }

  void OpenSwathHelper::selectSwathTransitions(const OpenSwath::LightTargetedExperiment& targeted_exp,
                                               OpenSwath::LightTargetedExperiment& transition_exp_used, double min_upper_edge_dist,
                                               double lower, double upper)
  {
    std::set<std::string> matching_compounds;
    for (Size i = 0; i < targeted_exp.transitions.size(); i++)
    {
      const OpenSwath::LightTransition& tr = targeted_exp.transitions[i];
      if (lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < upper &&
          std::fabs(upper - tr.getPrecursorMZ()) >= min_upper_edge_dist)
      {
        transition_exp_used.transitions.push_back(tr);
        matching_compounds.insert(tr.getPeptideRef());
      }
    }
    std::set<std::string> matching_proteins;
    for (Size i = 0; i < targeted_exp.compounds.size(); i++)
    {
      if (matching_compounds.find(targeted_exp.compounds[i].id) != matching_compounds.end())
      {
        transition_exp_used.compounds.push_back( targeted_exp.compounds[i] );
        for (Size j = 0; j < targeted_exp.compounds[i].protein_refs.size(); j++)
        {
          matching_proteins.insert(targeted_exp.compounds[i].protein_refs[j]);
        }
      }
    }
    for (Size i = 0; i < targeted_exp.proteins.size(); i++)
    {
      if (matching_proteins.find(targeted_exp.proteins[i].id) != matching_proteins.end())
      {
        transition_exp_used.proteins.push_back( targeted_exp.proteins[i] );
      }
    }
  }

  std::pair<double,double> OpenSwathHelper::estimateRTRange(const OpenSwath::LightTargetedExperiment & exp)
  {
    if (exp.getCompounds().empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Input list of targets is empty.");
    }
    double max = exp.getCompounds()[0].rt;
    double min = exp.getCompounds()[0].rt;
    for (Size i = 0; i < exp.getCompounds().size(); i++)
    {
      if (exp.getCompounds()[i].rt < min) min = exp.getCompounds()[i].rt;
      if (exp.getCompounds()[i].rt > max) max = exp.getCompounds()[i].rt;
    }
    return std::make_pair(min,max);
  }

  OpenSwath::LightTargetedExperiment
  OpenSwathHelper::sampleExperiment(
    const OpenSwath::LightTargetedExperiment & exp,
    Size bins,
    Size peptides_per_bin,
    unsigned int seed,
    bool sort_by_intensity,
    double top_fraction,
    const std::unordered_set<std::string> & priority_peptides)
  {
    OPENMS_PRECONDITION(bins >= 1, "bins must be >= 1");
    OPENMS_PRECONDITION(peptides_per_bin >= 1, "peptides_per_bin must be >= 1");
    OPENMS_PRECONDITION(!sort_by_intensity || (top_fraction > 0.0 && top_fraction <= 1.0),
                        "top_fraction must be in (0,1] when sort_by_intensity is true");

    // 0) initial candidate selection: exclude decoys 
    std::vector<OpenSwath::LightCompound> candidates;
    std::vector<OpenSwath::LightCompound> priority_candidates;

    std::unordered_set<String> good_ids;
    for (auto & tr : exp.getTransitions())
    {
      if (!tr.decoy)
        good_ids.insert(tr.getPeptideRef());
    }
    
    // Separate compounds into priority and regular candidates
    for (auto & cmp : exp.getCompounds())
    {
      if (good_ids.count(cmp.id))
      {
        // Check if this peptide sequence is in the priority list
        if (!priority_peptides.empty() && priority_peptides.count(cmp.sequence))
        {
          priority_candidates.push_back(cmp);
        }
        else
        {
          candidates.push_back(cmp);
        }
      }
    }

    // 1) optionally sort by library intensities and trim to top fraction
    if (sort_by_intensity && top_fraction > 0.0 && top_fraction <= 1.0)
    {
      // sum intensities per peptide across all transitions
      std::unordered_map<String, double> intensity_sum;
      for (auto & tr : exp.getTransitions())
      {
        if (!tr.decoy)
        {
          intensity_sum[tr.getPeptideRef()] += tr.library_intensity;
        }
      }

      // sort regular candidates by descending sum
      std::sort(candidates.begin(), candidates.end(), [&](auto & a, auto & b) {
        return intensity_sum[a.id] > intensity_sum[b.id];
      });

      // trim to top N%
      Size max_keep = std::max<Size>(1, static_cast<Size>(candidates.size() * top_fraction));
      candidates.resize(std::min(max_keep, candidates.size()));
      
      // Also sort priority candidates by intensity (but don't trim them)
      std::sort(priority_candidates.begin(), priority_candidates.end(), [&](auto & a, auto & b) {
        return intensity_sum[a.id] > intensity_sum[b.id];
      });
    }

    // Combine all available candidates for sampling
    std::vector<OpenSwath::LightCompound> all_candidates;
    all_candidates.reserve(priority_candidates.size() + candidates.size());
    all_candidates.insert(all_candidates.end(), priority_candidates.begin(), priority_candidates.end());
    all_candidates.insert(all_candidates.end(), candidates.begin(), candidates.end());

    if (all_candidates.size() < 3)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Insufficient candidates for sampling: " + String(all_candidates.size()) +
                                         " found, minimum 3 required for meaningful iRT calibration.");
    }

    // 2) estimate RT range using all available candidates
    double rt_min = std::numeric_limits<double>::max();
    double rt_max = std::numeric_limits<double>::lowest();
    for (auto & cmp : all_candidates)
    {
      rt_min = std::min(rt_min, cmp.rt);
      rt_max = std::max(rt_max, cmp.rt);
    }
    double bin_width = (rt_max - rt_min) / static_cast<double>(bins);

    // 3) sample priority peptides first, then fill remaining quota with uniform sampling
    std::vector<OpenSwath::LightCompound> picked;
    std::unordered_set<std::string> picked_sequences; // Track which sequences we've already picked
    
    // Initialize OpenMS random shuffler
    Math::RandomShuffler rshuffler;
    rshuffler.seed(seed == 0 ? std::random_device{}() : seed);
    
    // First pass: sample priority peptides
    if (!priority_candidates.empty())
    {
      OPENMS_LOG_INFO << "Sampling " << priority_candidates.size() 
                      << " priority peptides from the input experiment" << std::endl;
      
      for (Size b = 0; b < bins; ++b)
      {
        double lo = rt_min + b * bin_width;
        double hi = (b + 1 == bins ? rt_max : lo + bin_width);
        std::vector<OpenSwath::LightCompound> bucket;
        
        for (auto & cmp : priority_candidates)
        {
          if (cmp.rt >= lo && cmp.rt < hi && picked_sequences.find(cmp.sequence) == picked_sequences.end())
            bucket.push_back(cmp);
        }
        
        if (!bucket.empty())
        {
          rshuffler.portable_random_shuffle(bucket.begin(), bucket.end());
          Size take = std::min(peptides_per_bin, bucket.size());
          for (Size i = 0; i < take; ++i)
          {
            picked.push_back(bucket[i]);
            picked_sequences.insert(bucket[i].sequence);
          }
        }
      }
      
      OPENMS_LOG_INFO << "Successfully sampled " << picked.size() 
                      << " priority peptides" << std::endl;
    }
    
    // Second pass: fill remaining quota with regular sampling
    Size total_quota = bins * peptides_per_bin;
    if (picked.size() < total_quota && !candidates.empty())
    {
      OPENMS_LOG_INFO << "Filling remaining quota (" << (total_quota - picked.size()) 
                      << " peptides) from regular candidates" << std::endl;
      
      for (Size b = 0; b < bins; ++b)
      {
        double lo = rt_min + b * bin_width;
        double hi = (b + 1 == bins ? rt_max : lo + bin_width);
        
        // Count how many priority peptides we already have in this bin
        Size priority_in_bin = 0;
        for (auto & p : picked)
        {
          if (p.rt >= lo && p.rt < hi)
            priority_in_bin++;
        }
        
        // Calculate how many more we need from this bin
        Size needed = (peptides_per_bin > priority_in_bin) ? (peptides_per_bin - priority_in_bin) : 0;
        
        if (needed > 0)
        {
          std::vector<OpenSwath::LightCompound> bucket;
          for (auto & cmp : candidates)
          {
            if (cmp.rt >= lo && cmp.rt < hi && picked_sequences.find(cmp.sequence) == picked_sequences.end())
              bucket.push_back(cmp);
          }
          
          if (!bucket.empty())
          {
            rshuffler.portable_random_shuffle(bucket.begin(), bucket.end());
            Size take = std::min(needed, bucket.size());
            for (Size i = 0; i < take; ++i)
            {
              picked.push_back(bucket[i]);
              picked_sequences.insert(bucket[i].sequence);
            }
          }
        }
      }
    }

    // 4) assemble output experiment
    OpenSwath::LightTargetedExperiment out_exp;
    out_exp.compounds = picked;

    // copy matching transitions, excluding decoys if requested
    std::unordered_set<String> pep_ids;
    for (auto & cmp : picked)
      pep_ids.insert(cmp.id);
    for (auto & tr : exp.getTransitions())
    {
      if (pep_ids.count(tr.getPeptideRef()) && (!tr.decoy))
        out_exp.transitions.push_back(tr);
    }

    // copy associated proteins
    std::unordered_set<String> prot_ids;
    for (auto & cmp : picked)
      for (auto & pid : cmp.protein_refs)
        prot_ids.insert(pid);
    for (auto & prot : exp.getProteins())
    {
      if (prot_ids.count(prot.id))
        out_exp.proteins.push_back(prot);
    }

    return out_exp;
  }

  std::map<std::string, double> OpenSwathHelper::simpleFindBestFeature(
      const OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map,
      bool useQualCutoff, double qualCutoff)
  {
    std::map<std::string, double> result;
    for (const auto & trgroup_it : transition_group_map)
    {
      if (trgroup_it.second.getFeatures().empty() ) {continue;}

      // Find the feature with the highest score
      auto bestf = trgroup_it.second.getBestFeature();

      // Skip if we did not find a feature or do not exceed a certain quality
      if (useQualCutoff && bestf.getOverallQuality() < qualCutoff )
      {
        continue;
      }

      // If we have a found a best feature, add it to the vector
      String pepref = trgroup_it.second.getTransitions()[0].getPeptideRef();
      result[ pepref ] = bestf.getRT();
    }
    return result;
  }
}
