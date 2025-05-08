// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  /**
    @brief A helper class that is used by several OpenSWATH tools
  */
  class OPENMS_DLLAPI OpenSwathHelper
  {

public:

    /**
      @brief Compute unique precursor identifier

      Uses transition_group_id and isotope number to compute a unique precursor
      id of the form "groupID_Precursor_ix" where x is the isotope number, e.g.
      the monoisotopic precursor would become "groupID_Precursor_i0".

      @param[in] transition_group_id Unique id of the transition group (peptide/compound)
      @param[in] isotope Precursor isotope number

      @return Unique precursor identifier
    */
    static String computePrecursorId(const String& transition_group_id, int isotope)
    {
      return transition_group_id + "_Precursor_i" + String(isotope);
    }

    /**
      @brief Compute transition group id

      Uses the unique precursor identifier to compute the transition group id
      (peptide/compound identifier), reversing the operation performed by
      computePrecursorId().

      @param[in] precursor_id Precursor identifier as computed by computePrecursorId()

      @return Original transition group id
    */
    static String computeTransitionGroupId(const String& precursor_id)
    {
      std::vector<String> substrings;
      precursor_id.split("_", substrings);

      if (substrings.size() == 3) return substrings[0];
      else if (substrings.size() > 3)
      {
        String r;
        for (Size k = 0; k < substrings.size() - 2; k++) r += substrings[k] + "_";
        return r.prefix(r.size() - 1);
      }
      return "";
    }

    /**
      @brief Select transitions between lower and upper and write them into the new TargetedExperiment

      Version for the OpenMS TargetedExperiment

      @param[in] targeted_exp Transition list for selection
      @param[out] selected_transitions Selected transitions for SWATH window
      @param[in] min_upper_edge_dist Distance in Th to the upper edge
      @param[in] lower Lower edge of SWATH window (in Th)
      @param[in] upper Upper edge of SWATH window (in Th)
    */
    static void selectSwathTransitions(const OpenMS::TargetedExperiment& targeted_exp,
                                       OpenMS::TargetedExperiment& selected_transitions,
                                       double min_upper_edge_dist,
                                       double lower, double upper);

    /**
      @brief Select transitions between lower and upper and write them into the new TargetedExperiment

      Version for the LightTargetedExperiment

      @param[in] targeted_exp Transition list for selection
      @param[out] selected_transitions Selected transitions for SWATH window
      @param[in] min_upper_edge_dist Distance in Th to the upper edge
      @param[in] lower Lower edge of SWATH window (in Th)
      @param[in] upper Upper edge of SWATH window (in Th)
    */
    static void selectSwathTransitions(const OpenSwath::LightTargetedExperiment& targeted_exp,
                                       OpenSwath::LightTargetedExperiment& selected_transitions,
                                       double min_upper_edge_dist,
                                       double lower, double upper);
    /**
     @brief Match transitions with their "best" window across m/z and ion mobility, save results in a vector.

     @param[in] transition_exp Transition list for selection
     @param[out] tr_win_map Mapping from transition (index) to the best matching entry in @p swath_maps
     @param[in] min_upper_edge_dist Distance in Th to the upper edge
     @param[in] swath_maps vector of SwathMap objects defining mz and im bounds
    */
    static void selectSwathTransitionsPasef(const OpenSwath::LightTargetedExperiment& transition_exp, std::vector<int>& tr_win_map,
		                     double min_upper_edge_dist, const std::vector< OpenSwath::SwathMap > & swath_maps);

    /**
      @brief Get the lower / upper offset for this SWATH map and do some sanity checks

     
      Sanity check for the whole map:
       - all scans need to have exactly one precursor
       - all scans need to have the same MS levels (otherwise extracting an XIC
         from them makes no sense)
       - all scans need to have the same precursor isolation window (otherwise
         extracting an XIC from them makes no sense)

      @param[in] swath_map Input SWATH map to check
      @param[out] lower Lower edge of SWATH window (in Th)
      @param[out] upper Upper edge of SWATH window (in Th)
      @param[out] center Center of SWATH window (in Th)

      @throw IllegalArgument exception if the sanity checks fail.
    */
    static void checkSwathMap(const OpenMS::PeakMap& swath_map,
                              double& lower, double& upper, double& center);

    /**
      @brief Check the map and select transition in one function

      Computes lower and upper offset for the SWATH map and performs some
      sanity checks (see checkSwathMap()). Then selects transitions.

      @param[in] exp Input SWATH map to check
      @param[in] targeted_exp Transition list for selection
      @param[out] selected_transitions Selected transitions for SWATH window
      @param[in] min_upper_edge_dist Distance in Th to the upper edge
    */
    template <class TargetedExperimentT>
    static bool checkSwathMapAndSelectTransitions(const OpenMS::PeakMap& exp,
                                                  const TargetedExperimentT& targeted_exp,
                                                  TargetedExperimentT& selected_transitions,
                                                  double min_upper_edge_dist)
    {
      if (exp.empty() || exp[0].getPrecursors().empty())
      {
        std::cerr << "WARNING: File " << exp.getLoadedFilePath()
                  << " does not have any experiments or any precursors. Is it a SWATH map? "
                  << "I will move to the next map."
                  << std::endl;
        return false;
      }
      double upper, lower, center;
      OpenSwathHelper::checkSwathMap(exp, lower, upper, center);
      OpenSwathHelper::selectSwathTransitions(targeted_exp, selected_transitions, min_upper_edge_dist, lower, upper);
      if (selected_transitions.getTransitions().size() == 0)
      {
        std::cerr << "WARNING: For File " << exp.getLoadedFilePath()
                  << " no transition were within the precursor window of " << lower << " to " << upper
                  << std::endl;
        return false;
      }
      return true;

    }

    /**
      @brief Computes the min and max retention time value
      
      Estimate the retention time span of a targeted experiment by returning
      the min/max values in retention time as a pair.

      @return A std::pair that contains (min,max)

    */
    static std::pair<double,double> estimateRTRange(const OpenSwath::LightTargetedExperiment & exp);

    /**
      @brief Returns the feature with the highest score for each transition group.
      
      Simple method to extract the best feature for each transition group (e.g.
      for RT alignment). A quality cutoff can be used to skip some low-quality
      features altogether.

      @param[in] transition_group_map Input data containing the picked and scored map
      @param useQualCutoff Whether to apply a quality cutoff to the data
      @param qualCutoff What quality cutoff should be applied (all data above the cutoff will be kept)

      @return Result of the best scoring peaks (stored as map of peptide id and RT)

    */
    static std::map<std::string, double> simpleFindBestFeature(const OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
                                                               bool useQualCutoff = false,
                                                               double qualCutoff = 0.0);
  };

} // namespace OpenMS

