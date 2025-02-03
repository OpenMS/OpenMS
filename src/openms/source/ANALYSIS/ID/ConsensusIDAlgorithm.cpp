// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h> // for "OPENMS_PRECONDITION"
#include <OpenMS/PROCESSING/ID/IDFilter.h>

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
    defaults_.setValidStrings("filter:count_empty", {"true","false"});

    defaults_.setValue("filter:keep_old_scores", "false", "if set, keeps the original scores as user params");
    defaults_.setValidStrings("filter:keep_old_scores", {"true","false"});

    defaultsToParam_();
  }


  ConsensusIDAlgorithm::~ConsensusIDAlgorithm() = default;


  void ConsensusIDAlgorithm::updateMembers_()
  {
    considered_hits_ = param_.getValue("filter:considered_hits");
    min_support_ = param_.getValue("filter:min_support");
    count_empty_ = (param_.getValue("filter:count_empty") == "true");
    keep_old_scores_ = (param_.getValue("filter:keep_old_scores") == "true");
  }


  void ConsensusIDAlgorithm::apply(vector<PeptideIdentification>& ids,
                                   const map<String, String>& se_info,
                                   Size number_of_runs)
  {
    // abort if no IDs present
    if (ids.empty())
    {
      return;
    }

    number_of_runs_ = (number_of_runs != 0) ? number_of_runs : ids.size();

    // prepare data here, so that it doesn't have to happen in each algorithm:
    for (PeptideIdentification& pep : ids)
    {
      pep.sort();
      if ((considered_hits_ > 0) &&
          (pep.getHits().size() > considered_hits_))
      {
        pep.getHits().resize(considered_hits_);
      }
    }
    // make sure there are no duplicated hits (by sequence):
    IDFilter::removeDuplicatePeptideHits(ids, true);

    SequenceGrouping results;
    apply_(ids, se_info, results); // actual (subclass-specific) processing

    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(score_type);
    ids[0].setHigherScoreBetter(higher_better);
    for (SequenceGrouping::iterator res_it = results.begin(); 
         res_it != results.end(); ++res_it)
    {
      // filter by "support" value:
      if (res_it->second.support < min_support_) continue;
      PeptideHit hit;
      hit.setMetaValue("consensus_support", res_it->second.support);
      if (!res_it->second.target_decoy.empty())
        hit.setMetaValue("target_decoy", res_it->second.target_decoy);
      hit.setSequence(res_it->first);
      hit.setCharge(res_it->second.charge);
      hit.setScore(res_it->second.final_score);
      for (auto& ev : res_it->second.evidence)
      {
        hit.addPeptideEvidence(ev);
      }

      if (keep_old_scores_)
      {
        for (Size s = 0; s < res_it->second.scores.size(); ++s)
        {
          //TODO add SE name
          hit.setMetaValue(res_it->second.types[s]+"_score", res_it->second.scores[s]);
        }
      }
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      OPENMS_LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
    ids[0].assignRanks();
  }

  void ConsensusIDAlgorithm::apply(vector<PeptideIdentification>& ids,
                                   Size number_of_runs)
  {
    const auto empty = map<String,String>();
    apply(ids, empty, number_of_runs);
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
