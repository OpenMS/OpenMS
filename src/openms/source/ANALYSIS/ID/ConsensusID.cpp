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
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen and others $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusID.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <cmath>
#include <map>
#include <numeric> // for "accumulate"

// Extend SeqAn by a user-define scoring matrix.
namespace seqan
{

  // We have to create a new specialization of the _ScoringMatrix class
  // for amino acids.  For this, we first create a new tag.
  struct PAM30MS {};

  // Then, we specialize the class _ScoringMatrix.
  template <>
  struct ScoringMatrixData_<int, AminoAcid, PAM30MS>
  {
    enum
    {
      VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
      TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };
    static inline int const* getData()
    {
      // DEFINE the PAM30MS matrix from Habermann et al. MCP. 2004.
      static int const _data[TAB_SIZE] =
      {
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -17,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, 1,
      };

      return _data;
    }

  };
} // namespace seqan

using namespace std;

#define DEBUG_ID_CONSENSUS
#undef  DEBUG_ID_CONSENSUS

namespace OpenMS
{
  ConsensusID::ConsensusID() :
    DefaultParamHandler("ConsensusID")
  {
    defaults_.setValue("algorithm", "PEPMatrix", "Algorithm used for the consensus scoring.\n"
                       "PEPMatrix -- calculates a consensus score based on posterior error probabilities (PEPs) and peptide similarity scoring. This algorithm uses the PAM30MS matrix to score the similarity of sequences not listed by all search engines. Make sure that the scores for all peptide IDs are PEPs!\n"
                       "PEPIons -- calculates a consensus score based on posterior error probabilities (PEPs) and fragment ion similarities. Make sure that the scores for all peptide IDs are PEPs!\n"
                       "best -- uses the best score of any search engine as the consensus score of each peptide ID. Make sure that all peptide IDs use the same score type!\n"
                       "average -- uses the average score of all search engines as the consensus score of each peptide ID. Make sure that all peptide IDs use the same score type!\n"
                       "rank -- calculates a consensus score based on the ranks of peptide IDs in results of the different search engines. The final score is in the range (0, 1], with 1 being the best score. The input peptide IDs do not need to have the same score type.\n"
);
    defaults_.setValidStrings("algorithm", ListUtils::create<String>("PEPMatrix,PEPIons,best,average,rank"));
    defaults_.setValue("considered_hits", 10, "The number of top hits that are used for the consensus scoring ('0' for all hits).");
    defaults_.setMinInt("considered_hits", 0);
    defaults_.setValue("PEPIons:MinNumberOfFragments", 2, "The minimal number of similar (between two suggested sequences) fragment ion masses that is necessary to evaluate the shared peak count similarity (SPC).");
    defaults_.setMinInt("PEPIons:MinNumberOfFragments", 0);
    defaults_.setValue("PEPMatrix:common", 1.1, "Similarity threshold to accept the best score. E.g. for a given spectrum: engine 1 -> pep 1 with score x1 and engine2 -> pep2 with score x2. The better score from {x1,x2} will be used if the degree of similarity between pep1 and pep2 >= common, Note that 0 <= degree of similarity <= 1. Values > 1 will disable this option.");
    defaults_.setMinFloat("PEPMatrix:common", 0);
    defaults_.setMaxFloat("PEPMatrix:common", 1.1);
    defaults_.setValue("PEPMatrix:penalty", 5, "The alignment gap penalty as a positive integer (the same penalty will be used for opening and extension)");
    defaults_.setMinInt("PEPMatrix:penalty", 1);
    defaults_.setValue("PEPIons:common", 1.1, "Similarity threshold to accept the best score. E.g. for a given spectrum: engine 1 -> pep 1 with score x1 and engine2 -> pep2 with score x2. The better score from {x1,x2} will be used if the degree of similarity between pep1 and pep2 >= common, Note that 0 <= degree of similarity <= 1. Values > 1 will disable this option.");
    defaults_.setMinFloat("PEPIons:common", 0);
    defaults_.setMaxFloat("PEPIons:common", 1.1);
    defaults_.setValue("rank:number_of_runs", 0, "The number of runs used as input. This information is used in the 'rank' algorithm to compute the new scores. If not given, the number of input identifications is used.");
    defaults_.setMinInt("rank:number_of_runs", 0);

    defaultsToParam_();
  }

  void ConsensusID::apply(vector<PeptideIdentification>& ids)
  {
    // abort if no IDs present
    if (ids.empty())
    {
      return;
    }

    // prepare data here, so that it doesn't have to happen in each algorithm:
    Size considered_hits = param_.getValue("considered_hits");
    for (vector<PeptideIdentification>::iterator pep_it = ids.begin(); 
         pep_it != ids.end(); ++pep_it)
    {
      pep_it->sort();
      if ((considered_hits > 0) && (pep_it->getHits().size() > considered_hits))
      {
        pep_it->getHits().resize(considered_hits);
      }
    }

    String algorithm = param_.getValue("algorithm");

    if (algorithm == "PEPMatrix")
    {
      PEPMatrix_(ids);
    }
    else if (algorithm == "PEPIons")
    {
      PEPIons_(ids);
    }
    else if (algorithm == "rank")
    {
      rank_(ids);
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
          " is better) and cannot be compared meaningfully";
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
          // compare charge states:
          Int charge = hit_it->getCharge();
          if (pos->second.first == 0)
          {
            pos->second.first = charge;
          }
          else if ((charge != 0) && (pos->second.first != charge))
          {
            String msg = "Conflicting charge states found for peptide '" +
              pos->first.toString() + "': " + String(pos->second.first) +
              ", " + String(charge);
            throw Exception::InvalidValue(__FILE__, __LINE__,
                                          __PRETTY_FUNCTION__, msg, 
                                          String(charge));
          }
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


  void ConsensusID::rank_(vector<PeptideIdentification>& ids)
  {
    // The idea here is that each peptide hit (sequence) gets assigned a score
    // from each ID run, based on its rank in the list of search results.
    // The best hit of a run will receive score 0, the second best 1, etc. up to
    // a score of N = considered_hits - 1 for the last hit (if there are that
    // many). A hit that was not observed in an ID run receives a score of
    // N + 1 = considered_hits from that run. In the end the scores for each
    // hit (sequence) are averaged and normalized to the range from 0
    // (exclusive, worst) to 1 (inclusive, best).

    Size number_of_runs = param_.getValue("rank:number_of_runs");
    if (number_of_runs == 0) number_of_runs = ids.size();
    Size considered_hits = param_.getValue("considered_hits");
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
    ids[0].setScoreType("Consensus_rank");
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


  void ConsensusID::PEPMatrix_(vector<PeptideIdentification>& ids)
  {
    Map<AASequence, vector<double> > scores;

    int penalty = (UInt)param_.getValue("PEPMatrix:penalty");
    double common = (double)param_.getValue("PEPMatrix:common");

    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();

    for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
    {
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - ID run" << endl;
#endif

      //check the score type
      if (id->getScoreType() != score_type)
      {
        cerr << "Warning: You are working with different types of scores: '" << score_type << "' and '" << id->getScoreType() << "'" << endl;
      }
      if (id->isHigherScoreBetter() != higher_better)
      {
        cerr << "Warning: The score of the identifications have disagreeing score orientation!" << endl;
      }
      if (higher_better)
      {
        cerr << "You need to calculate posterior error probabilities as input scores!" << endl;
      }

      //iterate over the hits
      for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end(); ++hit)
      {
        double a_score = (double)hit->getScore();
        double a_sim = 1;
        double NumberAnnots = 1;


        set<String> myset;
        for (vector<PeptideHit>::const_iterator t = id->getHits().begin(); t != id->getHits().end(); ++t)
        {
          if (myset.find(t->getMetaValue("scoring")) == myset.end() && hit->getMetaValue("scoring") != t->getMetaValue("scoring"))
          {
            double a = 0;
            double zz = 0;
            vector<double> z;
            z.push_back((double)hit->getScore());
            // find the same or most similar peptide sequence in lists from other search engines
            for (vector<PeptideHit>::const_iterator tt = id->getHits().begin(); tt != id->getHits().end(); ++tt)
            {
              if (!(hit->getMetaValue("scoring") != t->getMetaValue("scoring"))) exit(1);

              if (tt->getMetaValue("scoring") == t->getMetaValue("scoring"))
              {
                //use SEQAN similarity scoring
                typedef ::seqan::String< ::seqan::AminoAcid> TSequence;
                TSequence seq1 = tt->getSequence().toUnmodifiedString().c_str();
                TSequence seq2 = hit->getSequence().toUnmodifiedString().c_str();
/////////////////////////introduce scoring with PAM30MS
                typedef int TValue;
                typedef ::seqan::Score<TValue, ::seqan::ScoreMatrix< ::seqan::AminoAcid, ::seqan::Default> > TScoringScheme;
                TScoringScheme pam30msScoring(-penalty, -penalty);
                ::seqan::setDefaultScoreMatrix(pam30msScoring, ::seqan::PAM30MS());
/////////////////////////introduce scoring with PAM30MS

//You can also use normal mutation based matrices, such as BLOSUM or the normal PAM matrix
                ::seqan::Pam250 pam(-5, -5);
                //::seqan::Score<int, ::seqan::Pam<> > pam(30, -10, -10);
                ::seqan::Align<TSequence, ::seqan::ArrayGaps> align, self1, self2;
                ::seqan::resize(rows(align), 2);
                ::seqan::resize(rows(self1), 2);
                ::seqan::resize(rows(self2), 2);
                ::seqan::assignSource(row(align, 0), seq1);
                ::seqan::assignSource(row(align, 1), seq2);
                ::seqan::assignSource(row(self1, 0), seq1);
                ::seqan::assignSource(row(self1, 1), seq1);
                ::seqan::assignSource(row(self2, 0), seq2);
                ::seqan::assignSource(row(self2, 1), seq2);

                vector<double> temp;
                temp.push_back(globalAlignment(self1, pam30msScoring, ::seqan::NeedlemanWunsch()));
                temp.push_back(globalAlignment(self2, pam30msScoring, ::seqan::NeedlemanWunsch()));

                double c = (double)globalAlignment(align, pam30msScoring, ::seqan::NeedlemanWunsch());
                double b = *(min_element(temp.begin(), temp.end()));
                c /= b;
                if (c < 0)
                {
                  c = 0;
                }
                if (c > a)
                {
                  a = c;
                  if (a >= common)
                  {
                    z.push_back((double)tt->getScore());
                    zz = *(min_element(z.begin(), z.end()));
                  }
                  else
                  {
                    zz = (double)tt->getScore() * a;
                  }
                }
              }
            }
            if (a >= common)
            {
              a_sim += 1;
            }
            else
            {
              a_sim += a;
            }
            NumberAnnots += 1;
            a_score += zz;
            myset.insert(t->getMetaValue("scoring"));
          }
        }
        // the meta value similarity corresponds to the sum of the similarities.
        // Note: if similarity equals the number of search engines, the same peptide has been assigned by all engines
        //::std::cout <<hit->getSequence()<<" a_score="<<a_score<<" a_sim="<< a_sim <<::std::endl;
        vector<double> ScoreSim;
        //test
        ScoreSim.push_back(a_score / (a_sim * a_sim));
        ScoreSim.push_back(a_sim);
        ScoreSim.push_back(NumberAnnots);
        ScoreSim.push_back(hit->getCharge());
        scores.insert(make_pair(hit->getSequence(), ScoreSim));
      }
    }

    //Replace IDs by consensus
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(String("Consensus_PEPMatrix (") + score_type + ")");
    ids[0].setHigherScoreBetter(false);
    for (Map<AASequence, vector<double> >::const_iterator it = scores.begin(); it != scores.end(); ++it)
    {
      PeptideHit hit;
      hit.setSequence(it->first);
      hit.setScore(it->second[0]);
      hit.setMetaValue("similarity", it->second[1]);
      hit.setMetaValue("Number of annotations", it->second[2]);
      hit.setCharge(static_cast<Int>(it->second[3]));
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " " << hit.getScore() << endl;
#endif
    }
  }


  void ConsensusID::PEPIons_(vector<PeptideIdentification>& ids)
  {
    Map<AASequence, vector<double> > scores;

    UInt MinNumberOfFragments = (UInt)(param_.getValue("PEPIons:MinNumberOfFragments"));

    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();
    double common = (double)param_.getValue("PEPIons:common");

    for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
    {
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - ID run" << endl;
#endif

      //check the score type
      if (id->getScoreType() != score_type)
      {
        cerr << "Warning: You are working with different types of scores: '" << score_type << "' and '" << id->getScoreType() << "'" << endl;
      }
      if (id->isHigherScoreBetter() != higher_better)
      {
        cerr << "Warning: The score of the identifications have disagreeing score orientation!" << endl;
      }
      if (higher_better)
      {
        cerr << "You need to calculate posterior error probabilities as input scores!" << endl;
      }

      //iterate over the hits
      for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end(); ++hit)
      {
        //double a_score=(double)hit->getMetaValue("PEP");
        double a_score = (double)hit->getScore();
        double a_sim = 1;
        double NumberAnnots = 1;

        set<String> myset;
        for (vector<PeptideHit>::const_iterator t = id->getHits().begin(); t != id->getHits().end(); ++t)
        {
          if (myset.find(t->getMetaValue("scoring")) == myset.end() && hit->getMetaValue("scoring") != t->getMetaValue("scoring"))
          {
            double a = 0;
            UInt SumIonSeries = 2;
            double zz = 0;
            vector<double> z;
            z.push_back((double)hit->getScore());
            //find the same or most similar peptide sequence in lists from other search engines
            for (vector<PeptideHit>::const_iterator tt = id->getHits().begin(); tt != id->getHits().end(); ++tt)
            {
              PeptideHit k = *tt;
              if (hit->getMetaValue("scoring") != t->getMetaValue("scoring") && tt->getMetaValue("scoring") == t->getMetaValue("scoring"))
              {
                //use similarity of b and y ion series for scoring
                AASequence S1, S2;
                S1 = tt->getSequence();
                S2 = hit->getSequence();
                //compare b and y ion series of S1 and S2
                vector<double> Yions_S1, Yions_S2, Bions_S1, Bions_S2;
                for (UInt r = 1; r <= S1.size(); ++r)
                {
                  Yions_S1.push_back(S1.getPrefix(r).getMonoWeight());
                  Bions_S1.push_back(S1.getSuffix(r).getMonoWeight());
                }
                for (UInt r = 1; r <= S2.size(); ++r)
                {
                  Yions_S2.push_back(S2.getPrefix(r).getMonoWeight());
                  Bions_S2.push_back(S2.getSuffix(r).getMonoWeight());
                }
                UInt Bs = 0;
                UInt Ys = 0;
                for (UInt xx = 0; xx < Yions_S1.size(); ++xx)
                {
                  for (UInt yy = 0; yy < Yions_S2.size(); ++yy)
                  {
                    if (fabs(Yions_S1[xx] - Yions_S2[yy]) <= 1)
                    {
                      Ys += 1;
                    }
                  }
                }
                for (UInt xx = 0; xx < Bions_S1.size(); ++xx)
                {
                  for (UInt yy = 0; yy < Bions_S2.size(); ++yy)
                  {
                    if (fabs(Bions_S1[xx] - Bions_S2[yy]) <= 1)
                    {
                      Bs += 1;
                    }
                  }
                }
                vector<UInt> tmp;
                tmp.push_back(Ys);
                tmp.push_back(Bs);
                UInt sum_tmp;
                sum_tmp = Bs + Ys;
                //# matching ions/number of AASeqences(S1)
                double c, b;
                b = *(min_element(tmp.begin(), tmp.end()));
                c = b / S1.size();
                if (sum_tmp > SumIonSeries && sum_tmp > MinNumberOfFragments)
                {
                  SumIonSeries = sum_tmp;
                  a = c;
                  if (a >= common)
                  {
                    z.push_back((double)tt->getScore());
                    zz = *(min_element(z.begin(), z.end()));

                  }
                  else
                  {
                    zz = (double)tt->getScore() * a;
                  }
                }
              }
            }
            NumberAnnots += 1;
            a_score += zz;
            a_sim += a;
            myset.insert(t->getMetaValue("scoring"));
          }
        }
        //the meta value similarity corresponds to the sum of the similarities. Note that if similarity equals the number of search engines, the
        //same peptide has been assigned by all engines
        vector<double> ScoreSim;
        ScoreSim.push_back(a_score / (a_sim * a_sim));
        ScoreSim.push_back(a_sim);
        ScoreSim.push_back(NumberAnnots);
        ScoreSim.push_back(hit->getCharge());
        scores.insert(make_pair(hit->getSequence(), ScoreSim));
      }
    }

    //Replace IDs by consensus
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(String("Consensus_PEPIons (") + score_type + ")");
    ids[0].setHigherScoreBetter(false);
    for (Map<AASequence, vector<double> >::const_iterator it = scores.begin(); it != scores.end(); ++it)
    {
      PeptideHit hit;
      hit.setSequence(it->first);
      hit.setScore(it->second[0]);
      hit.setMetaValue("similarity", it->second[1]);
      hit.setMetaValue("Number of annotations", it->second[2]);
      hit.setCharge(static_cast<Int>(it->second[3]));
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " " << hit.getScore() << endl;
#endif
    }
  }

} // namespace OpenMS
