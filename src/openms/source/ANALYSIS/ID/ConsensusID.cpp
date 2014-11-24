// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>


#include <map>
#include <cmath>

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

/*

{
  6,-7,-4,-3,-6,-4,-2,-2,-7,-5,-6,-7,-5,-8,-2,0,-1,-13,-8,-2,-7,-6,0,-17,
  -7,8,-6,-10,-8,-2,-9,-9,-2,-5,-7,0,-4,-9,-4,-3,-6,-2,-10,-8,5,-1,0,-17,
  -4,-6,8,2,-11,-3,-2,-3,0,-5,-6,-1,-9,-9,-6,0,-2,-8,-4,-8,-4,-2,0,-17,
  -3,-10,2,8,-14,-2,2,-3,-4,-7,-10,-4,-11,-15,-8,-4,-5,-15,-11,-8,-7,-3,0,-17,
  -6,-8,-11,-14,10,-14,-14,-9,-7,-6,-11,-14,-13,-13,-8,-3,-8,-15,-4,-6,-11,-14,0,-17,
  -4,-2,-3,-2,-14,8,1,-7,1,-8,-7,-3,-4,-13,-3,-5,-5,-13,-12,-7,-3,4,0,-17,
  -2,-9,-2,2,-14,1,8,-4,-5,-5,-7,-4,-7,-14,-5,-4,-6,-17,-8,-6,-7,-2,0,-17,
  -2,-9,-3,-3,-9,-7,-4,6,-9,-11,-11,-7,-8,-9,-6,-2,-6,-15,-14,-5,-8,-7,0,-17,
  -7,-2,0,-4,-7,1,-5,-9,9,-9,-8,-6,-10,-6,-4,-6,-7,-7,-3,-6,-4,-3,0,-17,
  -5,-5,-5,-7,-6,-8,-5,-11,-9,8,5,-6,-1,-2,-8,-7,-2,-14,-6,2,-6,-7,0,-17,
  -6,-7,-6,-10,-11,-7,-7,-11,-8,5,5,-7,0,-3,-8,-8,-5,-10,-7,0,-7,-7,0,-17,
  -7,0,-1,-4,-14,-3,-4,-7,-6,-6,-7,7,-2,-14,-6,-4,-3,-12,-9,-9,5,4,0,-17,
  -5,-4,-9,-11,-13,-4,-7,-8,-10,-1,0,-2,11,-4,-8,-5,-4,-13,-11,-1,-3,-3,0,-17,
  -8,-9,-9,-15,-13,-13,-14,-9,-6,-2,-3,-14,-4,9,-10,-6,-9,-4,2,-8,-12,-14,0,-17,
  -2,-4,-6,-8,-8,-3,-5,-6,-4,-8,-8,-6,-8,-10,8,-2,-4,-14,-13,-6,-5,-5,0,-17,
  0,-3,0,-4,-3,-5,-4,-2,-6,-7,-8,-4,-5,-6,-2,6,0,-5,-7,-6,-4,-5,0,-17,
  -1,-6,-2,-5,-8,-5,-6,-6,-7,-2,-5,-3,-4,-9,-4,0,7,-13,-6,-3,-5,-4,0,-17,
  -13,-2,-8,-15,-15,-13,-17,-15,-7,-14,-10,-12,-13,-4,-14,-5,-13,13,-5,-15,-7,-13,0,-17,
  -8,-10,-4,-11,-4,-12,-8,-14,-3,-6,-7,-9,-11,2,-13,-7,-6,-5,10,-7,-10,-11,0,-17,
  -2,-8,-8,-8,-6,-7,-6,-5,-6,2,0,-9,-1,-8,-6,-6,-3,-15,-7,7,-9,-8,0,-17,
  -7,5,-4,-7,-11,-3,-7,-8,-4,-6,-7,5,-3,-12,-5,-4,-5,-7,-10,-9,5,1,0,-17,
  -6,-1,-2,-3,-14,4,-2,-7,-3,-7,-7,4,-3,-14,-5,-5,-4,-13,-11,-8,1,4,0,-17,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-17,
  -17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,1,
              };*/
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
                                                 "ranked -- reorders the hits according to a consensus score computed from the ranks in the input runs. The score is normalized to the interval (0,100). The PeptideIdentifications do not need to have the same score type.\n"
                                                 "average -- reorders the hits according to the average score of the input runs. Make sure to use PeptideIdentifications with the same score type only!\n"
                                                 "PEPMatrix -- calculates a consensus score based on posterior error probabilities and scoring matrices for siimilarity. This algorithm uses the PAM30MS matrix to score sequences not listed by all engines. Make sure to use PeptideIdentifications with score types converted to PEPs only!\n"
                                                 "PEPIons -- calculates a consensus score based on posterior error probabilities and fragment ion siimilarity. Make sure to use PeptideIdentifications with score types converted to PEPs only!\n"
                                                 "Minimum -- calculates a consensus score based on the minimal score. Make sure to use PeptideIdentifications with score types converted to PEPs only!\n");
    defaults_.setValidStrings("algorithm", ListUtils::create<String>("ranked,average,PEPMatrix,PEPIons,Minimum"));
    defaults_.setValue("considered_hits", 10, "The number of top hits that are used for the consensus scoring.");
    defaults_.setMinInt("considered_hits", 1);
    defaults_.setValue("PEPIons:MinNumberOfFragments", 2, "The minimal number of similar (between two suggested sequences) fragment ion masses that is necessary to evaluate the shared peak count similarity (SPC).");
    defaults_.setMinInt("PEPIons:MinNumberOfFragments", 0);
    defaults_.setValue("number_of_runs", 0, "The number of runs used as input. This information is used in 'Ranked' and 'Average' to compute the new scores. If not given, the number of input identifications is taken.");
    defaults_.setMinInt("number_of_runs", 0);
    defaults_.setValue("PEPMatrix:common", 1.1, "Similarity threshold to accept the best score. E.g. for a given spectrum: engine 1 -> pep 1 with score x1 and engine2 -> pep2 with score x2. The better score from {x1,x2} will be used if the degree of similarity between pep1 and pep2 >= common, Note that 0 <= degree of similarity <= 1. Values > 1 will disable this option.");
    defaults_.setMinFloat("PEPMatrix:common", 0);
    defaults_.setMaxFloat("PEPMatrix:common", 1.1);
    defaults_.setValue("PEPMatrix:penalty", 5, "Give the gap penalty (the same penalty will be used for opening and extension) as a positive integer");
    defaults_.setValue("PEPIons:common", 1.1, "Similarity threshold to accept the best score. E.g. for a given spectrum: engine 1 -> pep 1 with score x1 and engine2 -> pep2 with score x2. The better score from {x1,x2} will be used if the degree of similarity between pep1 and pep2 >= common, Note that 0 <= degree of similarity <= 1. Values > 1 will disable this option.");
    defaults_.setMinFloat("PEPIons:common", 0);
    defaults_.setMaxFloat("PEPIons:common", 1.1);

    defaultsToParam_();
  }

  void ConsensusID::apply(vector<PeptideIdentification>& ids)
  {
    //Abort if no IDs present
    if (ids.empty())
    {
      return;
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
    else if (algorithm == "ranked")
    {
      ranked_(ids);
    }
    else if (algorithm == "average")
    {
      average_(ids);
    }
    else if (algorithm == "Minimum")
    {
      Minimum_(ids);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Algorithm '" + algorithm + "' was used but is not known! Please fix to something valid: " + ListUtils::concatenate(StringList(defaults_.getEntry("algorithm").valid_strings), ", ") + "!");
    }
    ids[0].assignRanks();

#ifdef DEBUG_ID_CONSENSUS
    const vector<PeptideHit>& hits2 = ids[0].getHits();
    for (Size i = 0; i < hits2.size(); ++i)
    {
      cout << "  " << hits2[i].getSequence() << " " << hits2[i].getScore() << endl;
    }
#endif
  }

  void ConsensusID::ranked_(vector<PeptideIdentification>& ids)
  {
    map<AASequence, double> scores;
    UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
    UInt number_of_runs = (UInt)(param_.getValue("number_of_runs"));
    String score_type = ids[0].getScoreType();

    //iterate over the different ID runs
    for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
    {
#ifdef DEBUG_ID_CONSENSUS
      cout << " - ID run" << endl;
#endif
      //make sure that the ranks are present
      id->assignRanks();
      //iterate over the hits
      UInt hit_count = 1;
      for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end() && hit_count <= considered_hits; ++hit)
      {
        if (scores.find(hit->getSequence()) == scores.end())
        {
#ifdef DEBUG_ID_CONSENSUS
          cout << " - New hit: " << hit->getSequence() << " " << hit->getRank() << endl;
#endif
          scores.insert(make_pair(hit->getSequence(), double(considered_hits + 1 - hit->getRank())));
        }
        else
        {
#ifdef DEBUG_ID_CONSENSUS
          cout << " - Added hit: " << hit->getSequence() << " " << hit->getRank() << endl;
#endif
          scores[hit->getSequence()] += (considered_hits + 1 - hit->getRank());
        }
        ++hit_count;
      }
    }
    //divide score by the maximum possible score and multiply with 100 => %
    Size max_score;
    if (number_of_runs == 0)
    {
      max_score = ids.size() * considered_hits;
    }
    else
    {
      max_score = number_of_runs * considered_hits;
    }
    for (map<AASequence, double>::iterator it = scores.begin(); it != scores.end(); ++it)
    {
      it->second = (it->second * 100.0f / max_score);
    }

    // replace IDs by consensus
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(String("Consensus_ranked (") + score_type + ")");
    ids[0].setHigherScoreBetter(true);

    for (map<AASequence, double>::const_iterator it = scores.begin(); it != scores.end(); ++it)
    {
      PeptideHit hit;
      hit.setSequence(it->first);
      hit.setScore(it->second);
      ids[0].insertHit(hit);
    }

  }

  void ConsensusID::average_(vector<PeptideIdentification>& ids)
  {
    map<AASequence, double> scores;
    UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
    UInt number_of_runs = (UInt)(param_.getValue("number_of_runs"));

    //store the score type (to make sure only IDs of the same type are averaged)
    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();

    //iterate over the different ID runs
    for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
    {
#ifdef DEBUG_ID_CONSENSUS
      cout << " - ID run" << endl;
#endif

      //check the score type
      if (id->getScoreType() != score_type)
      {
        cerr << "Warning: You are averaging different types of scores: '" << score_type << "' and '" << id->getScoreType() << "'" << endl;
      }
      if (id->isHigherScoreBetter() != higher_better)
      {
        cerr << "Warning: The score of the identifications have disagreeing score orientation!" << endl;
      }

      //make sure that the ranks are present
      id->assignRanks();
      //iterate over the hits
      UInt hit_count = 1;
      for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end() && hit_count <= considered_hits; ++hit)
      {
        if (scores.find(hit->getSequence()) == scores.end()) //.end zeigt auf ein Element nach dem letzten
        {
#ifdef DEBUG_ID_CONSENSUS
          cout << " - New hit: " << hit->getSequence() << " " << hit->getScore() << endl;
#endif
          scores.insert(make_pair(hit->getSequence(), hit->getScore()));
        }
        else
        {
#ifdef DEBUG_ID_CONSENSUS
          cout << " - Summed up: " << hit->getSequence() << " " << hit->getScore() << endl;
#endif
          scores[hit->getSequence()] += hit->getScore();
        }
        ++hit_count;
      }
    }
    //normalize score by number of id runs
    for (map<AASequence, double>::iterator it = scores.begin(); it != scores.end(); ++it)
    {
      if (number_of_runs == 0)
      {
        it->second = (it->second / ids.size());
      }
      else
      {
        it->second = (it->second / number_of_runs);
      }
    }

    //Replace IDs by consensus
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(String("Consensus_averaged (") + score_type + ")");
    ids[0].setHigherScoreBetter(higher_better);
    for (map<AASequence, double>::const_iterator it = scores.begin(); it != scores.end(); ++it)
    {
      PeptideHit hit;
      hit.setSequence(it->first);
      hit.setScore(it->second);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      cout << " - Output hit: " << hit.getSequence() << " " << hit.getScore() << endl;
#endif

    }
  }

//////////////////////////////////////////PEPMatrix

  void ConsensusID::PEPMatrix_(vector<PeptideIdentification>& ids)
  {
    Map<AASequence, vector<double> > scores;

    UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
    int penalty = (UInt)param_.getValue("PEPMatrix:penalty");
    double common = (double)param_.getValue("PEPMatrix:common");

    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();

    for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
    {
#ifdef DEBUG_ID_CONSENSUS
      cout << " - ID run" << endl;
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

      //make sure that the ranks are present
      id->assignRanks();

      //iterate over the hits
      UInt hit_count = 1;
      for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end() && hit_count <= considered_hits; ++hit)
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
                typedef::seqan::String< ::seqan::AminoAcid > TSequence;
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
        ++hit_count;
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
      hit.setCharge(it->second[3]);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      cout << " - Output hit: " << hit.getSequence() << " " << hit.getScore() << endl;
#endif

    }
  }

///////////////////////////////////////////////PEPIons

  void ConsensusID::PEPIons_(vector<PeptideIdentification>& ids)
  {
    Map<AASequence, vector<double> > scores;

    UInt considered_hits = (UInt)(param_.getValue("considered_hits"));
    UInt MinNumberOfFragments = (UInt)(param_.getValue("PEPIons:MinNumberOfFragments"));

    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();
    double common = (double)param_.getValue("PEPIons:common");

    for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
    {
#ifdef DEBUG_ID_CONSENSUS
      cout << " - ID run" << endl;
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

      //make sure that the ranks are present
      id->assignRanks();

      //iterate over the hits
      UInt hit_count = 1;

      for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end() && hit_count <= considered_hits; ++hit)
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
        ++hit_count;
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
      hit.setCharge(it->second[3]);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      cout << " - Output hit: " << hit.getSequence() << " " << hit.getScore() << endl;
#endif

    }
  }

//////////////////////////////////////////////////////////////////////////////////Minimum
  void ConsensusID::Minimum_(vector<PeptideIdentification>& ids)
  {
    Map<AASequence, double> scores;

    UInt considered_hits = (UInt)(param_.getValue("considered_hits"));


    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();
    if (higher_better == true)
    {
      cerr << "Warning: The score orientation is not suitable to take the minimum as the best hit!" << endl;
    }

    for (vector<PeptideIdentification>::iterator id = ids.begin(); id != ids.end(); ++id)
    {
#ifdef DEBUG_ID_CONSENSUS
      cout << " - ID run" << endl;
#endif
      //make sure that the ranks are present
      id->assignRanks();

      //iterate over the hits
      UInt hit_count = 1;

      double a_score = 1;
      AASequence a_pep;
      for (vector<PeptideHit>::const_iterator hit = id->getHits().begin(); hit != id->getHits().end() && hit_count <= considered_hits; ++hit)
      {
        if (hit->getScore() < a_score)
        {
          a_score = hit->getScore();
          a_pep = hit->getSequence();
        }

      }

      scores.insert(make_pair(a_pep, a_score));
      ++hit_count;
    }

    //Replace IDs by consensus
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(String("Consensus_Minimum(") + score_type + ")");
    ids[0].setHigherScoreBetter(false);

    for (Map<AASequence, double>::const_iterator it = scores.begin(); it != scores.end(); ++it)
    {
      PeptideHit hit;
      hit.setSequence(it->first);
      hit.setScore(it->second);
      ids[0].insertHit(hit);
    }
#ifdef DEBUG_ID_CONSENSUS
    cout << " - Output hit: " << hit.getSequence() << " " << hit.getScore() << endl;
#endif

  }

  /* this is not used by anyone ... delete?? */
  void ConsensusID::mapIdentifications_(vector<PeptideIdentification>& sorted_ids, const vector<PeptideIdentification>& ids)
  {
    double mz_delta = 0.01;
    double rt_delta = 0.01;
    for (vector<PeptideIdentification>::const_iterator it1 = ids.begin(); it1 != ids.end(); ++it1)
    {
      double rt1(it1->getRT());
      double mz1(it1->getMZ());
      PeptideIdentification new_ids;
      for (vector<PeptideIdentification>::const_iterator it2 = it1 + 1; it2 != ids.end(); ++it2)
      {
        double rt2(it2->getRT());
        double mz2(it2->getMZ());
        if (fabs(rt1 - rt2) < rt_delta && fabs(mz1 - mz2) < mz_delta)
        {
          if (new_ids.empty())
          {
            new_ids = (*it1);
          }
          else
          {
            for (vector<PeptideHit>::const_iterator pit = it2->getHits().begin(); pit != it2->getHits().end(); ++pit)
            {
              new_ids.insertHit(*pit);
            }
          }
        }
      }
      sorted_ids.push_back(new_ids);
    }
    return;
  }

} // namespace OpenMS
