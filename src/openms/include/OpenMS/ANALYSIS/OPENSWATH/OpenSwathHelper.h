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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>

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
      @brief Get the lower / upper offset for this SWATH map and do some sanity checks

     
      Sanity check for the whole map:
       - all scans need to have exactly one precursor
       - all scans need to have the same MS levels (otherwise extracting an XIC
         from them makes no sense)
       - all scans need to have the same precursor isolation window (otherwise
         extracting an XIC from them makes no sense)

      @param[in] swath_map Input SWATH map to check
      @param[in] lower Lower edge of SWATH window (in Th)
      @param[in] upper Upper edge of SWATH window (in Th)

      @throw throws IllegalArgument exception if the sanity checks fail.
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
      if (exp.size() == 0 || exp[0].getPrecursors().size() == 0)
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

