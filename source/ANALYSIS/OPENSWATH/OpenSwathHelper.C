// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

  void OpenSwathHelper::checkSwathMap(const OpenMS::MSExperiment<Peak1D>& swath_map,
                                      double& lower, double& upper)
  {
    if (swath_map.size() == 0 || swath_map[0].getPrecursors().size() == 0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Swath map has no Spectra");
    }
    const std::vector<Precursor> prec = swath_map[0].getPrecursors();
    lower = prec[0].getIsolationWindowLowerOffset();
    upper = prec[0].getIsolationWindowUpperOffset();
    UInt expected_mslevel = swath_map[0].getMSLevel();

    for (Size index = 0; index < swath_map.size(); index++)
    {
      const std::vector<Precursor> prec = swath_map[index].getPrecursors();
      if (prec.size() != 1)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Scan " + String(index) + " does not have exactly one precursor.");
      }
      if (swath_map[index].getMSLevel() != expected_mslevel)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Scan " + String(index) + " if of a different MS level than the first scan.");
      }
      if (
        fabs(lower - prec[0].getIsolationWindowLowerOffset()) > 0.1 ||
        fabs(upper - prec[0].getIsolationWindowUpperOffset()) > 0.1
        )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Scan " + String(index) + " has a different precursor isolation window than the first scan.");
      }

    }
  }

  void OpenSwathHelper::selectSwathTransitions(const OpenSwath::LightTargetedExperiment& targeted_exp,
                                               OpenSwath::LightTargetedExperiment& transition_exp_used, double min_upper_edge_dist,
                                               double lower, double upper)
  {
    std::set<std::string> matching_peptides;
    for (Size i = 0; i < targeted_exp.transitions.size(); i++)
    {
      const OpenSwath::LightTransition& tr = targeted_exp.transitions[i];
      if (lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < upper &&
          std::fabs(upper - tr.getPrecursorMZ()) >= min_upper_edge_dist)
      {
        transition_exp_used.transitions.push_back(tr);
        matching_peptides.insert(tr.getPeptideRef());
      }
    }
    std::set<std::string> matching_proteins;
    for (Size i = 0; i < targeted_exp.peptides.size(); i++)
    {
      if (matching_peptides.find(targeted_exp.peptides[i].id) != matching_peptides.end())
      {
        transition_exp_used.peptides.push_back( targeted_exp.peptides[i] );
        matching_proteins.insert(targeted_exp.peptides[i].protein_ref);
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

}
