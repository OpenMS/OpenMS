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
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h> // for "OPENMS_PRECONDITION"
#include <OpenMS/FILTERING/ID/IDFilter.h>

using namespace std;

// #define DEBUG_ID_CONSENSUS
// #undef  DEBUG_ID_CONSENSUS

namespace OpenMS
{
  ConsensusIDAlgorithm::ConsensusIDAlgorithm() :
    DefaultParamHandler("ConsensusIDAlgorithm")
  {
    defaults_.setValue("filter:considered_hits", 0, "The number of top hits in each ID run that are considered for consensus scoring ('0' for all hits).");
    defaults_.setMinInt("filter:considered_hits", 0);

    defaults_.setValue("filter:min_support", 0.0, "For each peptide hit from an ID run, the fraction of other ID runs that must support that hit (otherwise it is removed).");
    defaults_.setMinFloat("filter:min_support", 0.0);
    defaults_.setMaxFloat("filter:min_support", 1.0);
    defaults_.setValue("filter:count_empty", "false", "Count empty ID runs (i.e. those containing no peptide hit for the current spectrum) when calculating 'min_support'?");
    defaults_.setValidStrings("filter:count_empty", ListUtils::create<String>("true,false"));

    defaultsToParam_();
  }


  ConsensusIDAlgorithm::~ConsensusIDAlgorithm()
  {
  }


  void ConsensusIDAlgorithm::updateMembers_()
  {
    considered_hits_ = param_.getValue("filter:considered_hits");
    min_support_ = param_.getValue("filter:min_support");
    count_empty_ = (param_.getValue("filter:count_empty") == "true");
  }


  void ConsensusIDAlgorithm::apply(vector<PeptideIdentification>& ids,
                                   Size number_of_runs)
  {
    // abort if no IDs present
    if (ids.empty())
    {
      return;
    }

    number_of_runs_ = (number_of_runs != 0) ? number_of_runs : ids.size();

    // prepare data here, so that it doesn't have to happen in each algorithm:
    for (vector<PeptideIdentification>::iterator pep_it = ids.begin(); 
         pep_it != ids.end(); ++pep_it)
    {
      pep_it->sort();
      if ((considered_hits_ > 0) &&
          (pep_it->getHits().size() > considered_hits_))
      {
        pep_it->getHits().resize(considered_hits_);
      }
    }
    // make sure there are no duplicated hits (by sequence):
    IDFilter::removeDuplicatePeptideHits(ids, true);

    SequenceGrouping results;
    apply_(ids, results); // actual (subclass-specific) processing

    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(score_type);
    ids[0].setHigherScoreBetter(higher_better);
    for (SequenceGrouping::iterator res_it = results.begin(); 
         res_it != results.end(); ++res_it)
    {
      OPENMS_PRECONDITION(!res_it->second.second.empty(),
                          "Consensus score for peptide required");
      PeptideHit hit;

      if (res_it->second.second.size() == 2)
      {
        // filter by "support" value:
        double support = res_it->second.second[1];
        if (support < min_support_) continue;
        hit.setMetaValue("consensus_support", support);
      }
      
      hit.setSequence(res_it->first);
      hit.setCharge(res_it->second.first);
      hit.setScore(res_it->second.second[0]);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      OPENMS_LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
    ids[0].assignRanks();
  }


  void ConsensusIDAlgorithm::compareChargeStates_(Int& recorded_charge, 
                                                  Int new_charge,
                                                  const AASequence& peptide)
  {
    if (recorded_charge == 0) // update recorded charge
    {
      recorded_charge = new_charge;
    }
    else if ((new_charge != 0) && (recorded_charge != new_charge))
    { // maybe TODO: calculate correct charge from prec. m/z and peptide mass?
      String msg = "Conflicting charge states found for peptide '" +
        peptide.toString() + "': " + String(recorded_charge) + ", " + 
        String(new_charge);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                    msg, String(new_charge));
    }
  }

} // namespace OpenMS
