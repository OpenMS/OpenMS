// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

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
