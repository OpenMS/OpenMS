// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2019.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDScoreGetterSetter.h>

using namespace std;

namespace OpenMS
{

  /** @ingroup getScoresFunctions
      * @brief For protein groups. Groups are target if at least one protein is target
      * Therefore it needs an unordered set of decoy accessions to evaluate that.
    */
  void IDScoreGetterSetter::getScores_(
      ScoreToTgtDecLabelPairs &scores_labels,
      const std::vector<ProteinIdentification::ProteinGroup> &grps,
      const std::unordered_set<string> &decoy_accs)
  {
    for (const auto &grp : grps)
    {
      double score = grp.probability;
      bool target = false;
      for (const auto &acc : grp.accessions)
      {
        // In groups you usually want to check if at least one member is a real target
        if (decoy_accs.find(acc) == decoy_accs.end())
        {
          target = true;
          break;
        }
      }
      scores_labels.emplace_back(score, target);
    }
  }


  /** @ingroup setScoresFunctions
  * @brief For protein groups. Unaffected by keep_decoy_proteins. Always keeps all for now @todo.
  * score_type and higher_better unused since ProteinGroups do not carry that information.
  * You have to assume that groups will always have the same scores as the ProteinHits
  */
  void IDScoreGetterSetter::setScores_(const map<double, double> &scores_to_FDR,
                                      vector <ProteinIdentification::ProteinGroup> &grps,
                                      const string & /*score_type*/,
                                      bool /*higher_better*/)
  {
    for (auto &grp : grps)
    {
      grp.probability = (scores_to_FDR.lower_bound(grp.probability)->second);
    }
  }
} // namespace std
