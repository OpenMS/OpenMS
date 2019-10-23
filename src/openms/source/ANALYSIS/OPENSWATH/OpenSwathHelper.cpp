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
        transition_exp_used.addTransition(tr);
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
