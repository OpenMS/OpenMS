// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/ANALYSIS/ID/ConsensusID.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <cmath>
#include <numeric> // for "accumulate"

using namespace std;

#define DEBUG_ID_CONSENSUS
#undef  DEBUG_ID_CONSENSUS

namespace OpenMS
{
  ConsensusID::ConsensusID() :
    DefaultParamHandler("ConsensusID")
  {
    defaults_.setValue("algorithm", "PEPMatrix", "Algorithm used for the consensus scoring.\n"
                       "PEPMatrix -- calculates a consensus score based on posterior error probabilities (PEPs) and peptide sequence similarities. This algorithm uses a substitution matrix to score the similarity of sequences not listed by all search engines. Make sure that the scores for all peptide IDs are PEPs!\n"
                       "PEPIons -- calculates a consensus score based on posterior error probabilities (PEPs) and fragment ion similarities. Make sure that the scores for all peptide IDs are PEPs!\n"
                       "best -- uses the best score of any search engine as the consensus score of each peptide ID. Make sure that all peptide IDs use the same score type!\n"
                       "average -- uses the average score of all search engines as the consensus score of each peptide ID. Make sure that all peptide IDs use the same score type!\n"
                       "ranks -- calculates a consensus score based on the ranks of peptide IDs in results of the different search engines. The final score is in the range (0, 1], with 1 being the best score. The input peptide IDs do not need to have the same score type.\n"
);
    defaults_.setValidStrings("algorithm", ListUtils::create<String>("PEPMatrix,PEPIons,best,average,ranks"));
    defaults_.setValue("considered_hits", 10, "The number of top hits that are used for the consensus scoring ('0' for all hits).");
    defaults_.setMinInt("considered_hits", 0);

    defaults_.setValue("PEPMatrix:matrix", "identity", "Substitution matrix to use for alignment-based similarity scoring");
    defaults_.setValidStrings("PEPMatrix:matrix", ListUtils::create<String>("identity,PAM30MS"));
    defaults_.setValue("PEPMatrix:penalty", 5, "Alignment gap penalty (the same value will be used for gap opening and extension)");
    defaults_.setMinInt("PEPMatrix:penalty", 1);

    defaults_.setValue("PEPIons:mass_tolerance", 0.5, "Maximum difference between fragment masses (in Da) for fragments to be considered 'shared' between peptides .");
    defaults_.setMinFloat("PEPIons:mass_tolerance", 0.0);    
    // is the "min_shared" parameter necessary/useful?
    defaults_.setValue("PEPIons:min_shared", 2, "The minimal number of 'shared' fragments (between two suggested peptides) that is necessary to evaluate the similarity based on shared peak count (SPC).");
    defaults_.setMinInt("PEPIons:min_shared", 1);

    defaults_.setValue("ranks:number_of_runs", 0, "The number of runs used as input. This information is used in the 'ranks' algorithm to compute the new scores. If not given, the number of input identifications is used.");
    defaults_.setMinInt("ranks:number_of_runs", 0);

    defaultsToParam_();

    ::seqan::resize(rows(alignment_), 2);
  }


  void ConsensusID::updateMembers_()
  {
    // alignment scoring using SeqAn/similarity matrices:
    String matrix = param_.getValue("PEPMatrix:matrix");
    int penalty = param_.getValue("PEPMatrix:penalty");
    scoring_method_ = SeqAnScore(-penalty, -penalty);
    if (matrix == "identity")
    {
      ::seqan::setDefaultScoreMatrix(scoring_method_, 
                                     ::seqan::AdaptedIdentity());
    }
    else if (matrix == "PAM30MS")
    {
      ::seqan::setDefaultScoreMatrix(scoring_method_, ::seqan::PAM30MS());
    }
    else
    {
      String msg = "Matrix '" + matrix + "' is not known! Valid choices are: "
        "'identity', 'PAM30MS'.";
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       msg);
    }

    // similarity scoring based on shared peak count:
    mass_tolerance_ = param_.getValue("PEPIons:mass_tolerance");
    min_shared_ = param_.getValue("PEPIons:min_shared");

    // other parameters:
    number_of_runs_ = param_.getValue("ranks:number_of_runs");
    considered_hits_ = param_.getValue("considered_hits");
    
    // new parameters may affect the similarity calculation, so clear cache:
    similarities_.clear();
  }


  void ConsensusID::apply(vector<PeptideIdentification>& ids)
  {
    // abort if no IDs present
    if (ids.empty())
    {
      return;
    }

    // @TODO: what if there's only one ID run?!

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

    String algorithm = param_.getValue("algorithm");

    if (algorithm == "PEPMatrix")
    {
      PEPMatrixOrIons_(ids, true);
    }
    else if (algorithm == "PEPIons")
    {
      PEPMatrixOrIons_(ids, false);
    }
    else if (algorithm == "ranks")
    {
      ranks_(ids);
    }
    else if (algorithm == "average")
    {
      average_(ids);
    }
    else if (algorithm == "best")
    {
      best_(ids);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Algorithm '" + algorithm + "' is not known! Valid choices are: " + ListUtils::concatenate(StringList(defaults_.getEntry("algorithm").valid_strings), ", "));
    }
    ids[0].assignRanks();

#ifdef DEBUG_ID_CONSENSUS
    const vector<PeptideHit>& hits2 = ids[0].getHits();
    for (Size i = 0; i < hits2.size(); ++i)
    {
      LOG_DEBUG << "  " << hits2[i].getSequence() << " "
                << hits2[i].getScore() << endl;
    }
#endif
  }


  void ConsensusID::compareChargeStates_(Int& recorded_charge, Int new_charge,
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
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
                                    msg, String(new_charge));
    }
  }


  void ConsensusID::groupHits_(vector<PeptideIdentification>& ids,
                               SequenceGrouping& grouping)
  {
    bool higher_better = ids[0].isHigherScoreBetter();
    set<String> score_types;

    for (vector<PeptideIdentification>::iterator pep_it = ids.begin();
         pep_it != ids.end(); ++pep_it)
    {
      if (pep_it->isHigherScoreBetter() != higher_better)
      {
        // scores with different orientations definitely aren't comparable:
        String hi_lo = higher_better ? "higher/lower" : "lower/higher";
        String msg = "Score types '" + ids[0].getScoreType() + "' and '" +
          pep_it->getScoreType() + "' have different orientations (" + hi_lo +
          " is better) and cannot be compared meaningfully.";
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                      msg, higher_better ? "false" : "true");
      }
      score_types.insert(pep_it->getScoreType());
      pep_it->sort();
      for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
           hit_it != pep_it->getHits().end(); ++hit_it)
      {
        const AASequence& seq = hit_it->getSequence();
        SequenceGrouping::iterator pos = grouping.find(seq);
        if (pos == grouping.end()) // new sequence
        {
          grouping[seq] = make_pair(hit_it->getCharge(), 
                                    vector<double>(1, hit_it->getScore()));
        }
        else // previously seen sequence
        { 
          compareChargeStates_(pos->second.first, hit_it->getCharge(),
                               pos->first);
          pos->second.second.push_back(hit_it->getScore());
        }
      }
    }
    if (score_types.size() > 1)
    {
      String types;
      types.concatenate(score_types.begin(), score_types.end(), "'/'");
      LOG_WARN << "Warning: Different score types for peptide hits found ('"
               << types << "'). If the scores are not comparable, "
               << "results will be meaningless." << endl;
    }
  }


  void ConsensusID::ranks_(vector<PeptideIdentification>& ids)
  {
    // The idea here is that each peptide hit (sequence) gets assigned a score
    // from each ID run, based on its rank in the list of search results.
    // The best hit of a run will receive score 0, the second best 1, etc. up to
    // a score of N = considered_hits - 1 for the last hit (if there are that
    // many). A hit that was not observed in an ID run receives a score of
    // N + 1 = considered_hits from that run. In the end the scores for each
    // hit (sequence) are averaged and normalized to the range from 0
    // (exclusive, worst) to 1 (inclusive, best).

    Size number_of_runs = (number_of_runs_ > 0) ? number_of_runs_ : ids.size();
    Size considered_hits = considered_hits_;
    bool set_considered_hits = (considered_hits == 0);

    for (vector<PeptideIdentification>::iterator pep_it = ids.begin();
         pep_it != ids.end(); ++pep_it)
    {
      pep_it->assignRanks();
      for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
           hit_it != pep_it->getHits().end(); ++hit_it)
      {
        // give each hit a score based on the search rank (counting from 0):
        hit_it->setScore(hit_it->getRank() - 1);
      }
      pep_it->setScoreType("rank");

      // if "considered_hits" wasn't set, we find the max. number of hits:
      if (set_considered_hits && (pep_it->getHits().size() > considered_hits))
      {
        considered_hits = pep_it->getHits().size();
      }
    }

    SequenceGrouping grouping;
    groupHits_(ids, grouping);

    ids.clear();
    ids.resize(1);
    // the original score type doesn't matter here, so don't include it:
    ids[0].setScoreType("Consensus_ranks");
    ids[0].setHigherScoreBetter(true);
    for (SequenceGrouping::iterator group_it = grouping.begin(); 
         group_it != grouping.end(); ++group_it)
    {
      PeptideHit hit;
      hit.setSequence(group_it->first);
      double sum_scores = accumulate(group_it->second.second.begin(),
                                     group_it->second.second.end(), 0.0);
      // add score contributions equivalent to "not found":
      sum_scores += ((number_of_runs - group_it->second.second.size()) *
                     considered_hits);
      // normalize to range 0-1:
      hit.setScore(1.0 - sum_scores / (considered_hits * number_of_runs));
      hit.setCharge(group_it->second.first);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
  }


  void ConsensusID::average_(vector<PeptideIdentification>& ids)
  {
    SequenceGrouping grouping;
    groupHits_(ids, grouping);
    
    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(String("Consensus_average (") + score_type + ")");
    ids[0].setHigherScoreBetter(higher_better);
    for (SequenceGrouping::iterator group_it = grouping.begin(); 
         group_it != grouping.end(); ++group_it)
    {
      PeptideHit hit;
      hit.setSequence(group_it->first);
      double sum_scores = accumulate(group_it->second.second.begin(),
                                     group_it->second.second.end(), 0.0);
      hit.setScore(sum_scores / group_it->second.second.size());
      hit.setCharge(group_it->second.first);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
  }


  void ConsensusID::best_(vector<PeptideIdentification>& ids)
  {
    SequenceGrouping grouping;
    groupHits_(ids, grouping);
    
    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(String("Consensus_best (") + score_type + ")");
    ids[0].setHigherScoreBetter(higher_better);
    for (SequenceGrouping::iterator group_it = grouping.begin(); 
         group_it != grouping.end(); ++group_it)
    {
      PeptideHit hit;
      hit.setSequence(group_it->first);
      if (higher_better)
      {
        hit.setScore(*max_element(group_it->second.second.begin(),
                                  group_it->second.second.end()));
      }
      else
      {
        hit.setScore(*min_element(group_it->second.second.begin(),
                                  group_it->second.second.end()));
      }
      hit.setCharge(group_it->second.first);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
  }


  double ConsensusID::getSimilarityMatrix_(AASequence seq1, AASequence seq2)
  {
    // here we cannot take modifications into account:
    String unmod_seq1 = seq1.toUnmodifiedString();
    String unmod_seq2 = seq2.toUnmodifiedString();
    if (unmod_seq1 == unmod_seq2) return 1.0;
    // order of sequences matters for cache look-up:
    if (unmod_seq1 > unmod_seq2) swap(unmod_seq1, unmod_seq2);
    seq1 = AASequence::fromString(unmod_seq1);
    seq2 = AASequence::fromString(unmod_seq2);
    pair<AASequence, AASequence> seq_pair = make_pair(seq1, seq2);
    SimilarityCache::iterator pos = similarities_.find(seq_pair);
    if (pos != similarities_.end()) return pos->second; // score found in cache
    
    // use SeqAn similarity scoring:
    SeqAnSequence seqan_seq1 = unmod_seq1.c_str();
    SeqAnSequence seqan_seq2 = unmod_seq2.c_str();
    // seq. 1 against itself:
    ::seqan::assignSource(row(alignment_, 0), seqan_seq1);
    ::seqan::assignSource(row(alignment_, 1), seqan_seq1);
    double score_self1 = globalAlignment(alignment_, scoring_method_,
                                         ::seqan::NeedlemanWunsch());
    // seq. 1 against seq. 2:
    ::seqan::assignSource(row(alignment_, 1), seqan_seq2);
    double score_sim = globalAlignment(alignment_, scoring_method_, 
                                       ::seqan::NeedlemanWunsch());
    // seq. 2 against itself:
    ::seqan::assignSource(row(alignment_, 0), seqan_seq2);
    double score_self2 = globalAlignment(alignment_, scoring_method_,
                                         ::seqan::NeedlemanWunsch());
    if (score_sim < 0)
    {
      score_sim = 0;
    }
    else
    {
      score_sim /= min(score_self1, score_self2); // normalize
    }
    similarities_[seq_pair] = score_sim; // cache the similarity score

    return score_sim;
  }


  double ConsensusID::getSimilarityIons_(AASequence seq1, AASequence seq2)
  {
    if (seq1 == seq2) return 1.0;
    // order of sequences matters for cache look-up:
    if (seq2 < seq1) swap(seq1, seq2); // "operator>" not defined
    pair<AASequence, AASequence> seq_pair = make_pair(seq1, seq2);
    SimilarityCache::iterator pos = similarities_.find(seq_pair);
    if (pos != similarities_.end()) return pos->second; // score found in cache

    // compare b and y ion series of seq. 1 and seq. 2:
    vector<double> ions1(2 * seq1.size()), ions2(2 * seq2.size());
    // b ions, seq. 1:
    ions1[0] = seq1.getPrefix(1).getMonoWeight(); // includes N-terminal mods
    // y ions, seq. 1:
    ions1[seq1.size()] = seq1.getSuffix(1).getMonoWeight(); // inc. C-term. mods
    for (Size i = 1; i < seq1.size(); ++i)
    {
      ions1[i] = ions1[i - 1] + seq1[i].getMonoWeight();
      ions1[seq1.size() + i] = (ions1[seq1.size() + i - 1] + 
                                seq1[seq1.size() - i - 1].getMonoWeight());
    }
    // b ions, seq. 2:
    ions2[0] = seq2.getPrefix(1).getMonoWeight(); // includes N-terminal mods
    // y ions, seq. 2:
    ions2[seq2.size()] = seq2.getSuffix(1).getMonoWeight(); // inc. C-term. mods
    for (Size i = 1; i < seq2.size(); ++i)
    {
      ions2[i] = ions2[i - 1] + seq2[i].getMonoWeight();
      ions2[seq2.size() + i] = (ions2[seq2.size() + i - 1] + 
                                seq2[seq2.size() - i - 1].getMonoWeight());
    }

    // now compare fragment masses from both sequences to find best matches
    // within the allowed tolerance; note that:
    // 1. we can be more efficient than comparing "all against all"
    // 2. an ion from seq. 2 may be the best match for two (similar) ions from
    // seq. 1 - then we want to count that ion only once, not twice 
    sort(ions1.begin(), ions1.end());
    sort(ions2.begin(), ions2.end());
    set<double> matches; // each best-matching ion counts only once
    vector<double>::iterator start = ions2.begin();
    // for each fragment in seq. 1...
    for (vector<double>::iterator it1 = ions1.begin(); it1 != ions1.end();
         ++it1)
    {
      // ...find fragments from seq. 2 that are within the mass tolerance:
      vector<double>::iterator lower = lower_bound(start, ions2.end(),
                                                   *it1 - mass_tolerance_);
      if (lower == ions2.end()) break; // all values are too low
      vector<double>::iterator upper = upper_bound(lower, ions2.end(),
                                                   *it1 + mass_tolerance_);
      double best_match = 0.0, best_diff = mass_tolerance_ + 1.0;
      // find ion from seq. 2 (*it2) that is closest to ion from seq. 1 (*it1):
      for (vector<double>::iterator it2 = lower; it2 != upper; ++it2)
      {
        double diff = fabs(*it1 - *it2);
        if (diff < best_diff)
        {
          best_diff = diff;
          best_match = *it2;
        }
      }
      if (best_diff <= mass_tolerance_) matches.insert(best_match);
      
      start = lower; // "*it1" is increasing, so lower bounds can't get lower
    }

    double score_sim = 0.0;
    if (matches.size() >= min_shared_)
    {
      score_sim = matches.size() / float(min(ions1.size(), ions2.size()));
    }
    similarities_[seq_pair] = score_sim; // cache the similarity score

    return score_sim;
  }


  void ConsensusID::PEPMatrixOrIons_(vector<PeptideIdentification>& ids,
                                     const bool matrix)
  {
    for (vector<PeptideIdentification>::iterator id = ids.begin();
         id != ids.end(); ++id)
    {
      if (id->getScoreType() != "Posterior Error Probability")
      {
        String msg = "Score type must be 'Posterior Error Probablity'";
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                      msg, id->getScoreType());
      }
    }

    // mapping: peptide sequence -> (charge, [consensus score, similarity])
    SequenceGrouping results;
    
    for (vector<PeptideIdentification>::iterator id1 = ids.begin();
         id1 != ids.end(); ++id1)
    {
      for (vector<PeptideHit>::iterator hit1 = id1->getHits().begin();
           hit1 != id1->getHits().end(); ++hit1)
      {
        // have we scored this sequence already? if yes, skip:
        SequenceGrouping::iterator pos = results.find(hit1->getSequence());
        if (pos != results.end())
        { 
          compareChargeStates_(pos->second.first, hit1->getCharge(),
                               pos->first);
          continue;
        }
        
        // similarity scores and PEPs of best matches for all ID runs:
        vector<pair<double, double> > best_matches;
        best_matches.reserve(ids.size() - 1);
        for (vector<PeptideIdentification>::iterator id2 = ids.begin();
             id2 != ids.end(); ++id2)
        {
          if (id1 == id2) continue;
          
          // similarity scores and PEPs of all matches in current ID run
          // (to get the best match, we look for highest similarity, breaking
          // ties by better PEP - so we need to transform PEP so higher scores
          // are better, same as similarity):
          vector<pair<double, double> > current_matches;
          current_matches.reserve(id2->getHits().size());
          for (vector<PeptideHit>::iterator hit2 = id2->getHits().begin();
               hit2 != id2->getHits().end(); ++hit2)
          {
            // hopefully the compiler should optimize this by moving the check
            // ("matrix ?") outside of all the for-loops; if not:
            // maybe @TODO - if necessary, optimize this by:
            // a) turning "matrix" into a template argument and creating two
            // instantiations of this function, or
            // b) making the "getSimilarity..." functions static and passing in
            // a function pointer
            double sim_score = matrix ? 
              getSimilarityMatrix_(hit1->getSequence(), hit2->getSequence()) :
              getSimilarityIons_(hit1->getSequence(), hit2->getSequence());
            // use "1 - PEP" so higher scores are better (for "max_element"):
            current_matches.push_back(make_pair(sim_score,
                                                1.0 - hit2->getScore()));
          }
          best_matches.push_back(*max_element(current_matches.begin(),
                                              current_matches.end()));
        }
        double score = hit1->getScore();
        double sum_sim = 1.0; // sum of similarity scores
        for (vector<pair<double, double> >::iterator it = best_matches.begin();
             it != best_matches.end(); ++it)
        {
          score += it->first * (1.0 - it->second); // undo "1 - PEP" transform
          sum_sim += it->first;
        }
        score /= (sum_sim * sum_sim);

        vector<double> scores(2, score);
        // normalize similarity score to range 0-1:
        scores[1] = (sum_sim - 1.0) / best_matches.size();

        results[hit1->getSequence()] = make_pair(hit1->getCharge(), scores);
      }
    }

    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(matrix ? "Consensus_PEPMatrix" : "Consensus_PEPIons");
    ids[0].setHigherScoreBetter(false);
    for (SequenceGrouping::iterator res_it = results.begin(); 
         res_it != results.end(); ++res_it)
    {
      PeptideHit hit;
      hit.setSequence(res_it->first);
      hit.setScore(res_it->second.second[0]);
      hit.setCharge(res_it->second.first);
      hit.setMetaValue("consensus_similarity", res_it->second.second[1]);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
  }

} // namespace OpenMS
