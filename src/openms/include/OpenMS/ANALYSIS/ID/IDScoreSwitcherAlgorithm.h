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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <vector>
#include <set>

namespace OpenMS
{

  class OPENMS_DLLAPI IDScoreSwitcherAlgorithm:
    public DefaultParamHandler
  {
  public:
    IDScoreSwitcherAlgorithm ();

    /// This is a rough hierarchy of possible score types in MS
    /// In an ideal case this should be reimplemented to follow
    /// ontology hierarchies as soon as e.g. MS-OBO is complete
    /// and we switched the Metavalues to CV terms.
    enum class ScoreType
    {
      RAW,
      RAW_EVAL,
      PP,
      PEP,
      FDR,
      QVAL,
    };

    /// Switches all main scores in all hits in @p id according to
    /// the settings in the param object of the switcher class
    template <typename IDType>
    void switchScores(IDType& id, Size& counter)
    {
      for (typename std::vector<typename IDType::HitType>::iterator hit_it = id.getHits().begin();
           hit_it != id.getHits().end(); ++hit_it, ++counter)
      {
        if (!hit_it->metaValueExists(new_score_))
        {
          std::stringstream msg;
          msg << "Meta value '" << new_score_ << "' not found for " << *hit_it;
          throw Exception::MissingInformation(__FILE__, __LINE__,
                                              OPENMS_PRETTY_FUNCTION, msg.str());
        }

        const String& old_score_meta = (old_score_.empty() ? id.getScoreType() :
                                 old_score_);
        const DataValue& dv = hit_it->getMetaValue(old_score_meta);
        if (!dv.isEmpty()) // meta value for old score already exists
        {
          if (fabs((double(dv) - hit_it->getScore()) * 2.0 /
                   (double(dv) + hit_it->getScore())) > tolerance_)
          {
            std::stringstream msg;
            msg << "Meta value '" << old_score_meta << "' already exists "
              << "with a conflicting value for " << *hit_it;
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          msg.str(), dv.toString());
          } // else: values match, nothing to do
        }
        else
        {
          hit_it->setMetaValue(old_score_meta, hit_it->getScore());
        }
        hit_it->setScore(hit_it->getMetaValue(new_score_));
      }
      id.setScoreType(new_score_type_);
      id.setHigherScoreBetter(higher_better_);
    }

    /// Looks at the first Hit of the given @p id and according to the given @p type ,
    /// deduces a fitting score and score direction to be switched to.
    /// Then tries to switch all hits.
    void switchToGeneralScoreType(std::vector<PeptideIdentification>& id, ScoreType type, Size& counter)
    {
      if (id.empty()) return;
      String t = findScoreType(id[0], type);
      if (t.empty())
      {
        String msg = "First encountered ID does not have the requested score type.";
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, msg);
      } 
      else if (t == id[0].getScoreType())
      {
        // we assume that all the other peptide ids
        // also already have the correct score set
        return;
      }

      if (t.hasSuffix("_score"))
      {
        new_score_type_ = t.chop(6);
      }
      else
      {
        new_score_type_ = t;
      }
      new_score_ = t;

      if (type != ScoreType::RAW && higher_better_ != type_to_better_[type])
      {
        OPENMS_LOG_WARN << "Requested non-raw score type does not match the expected score direction. Correcting!\n";
        higher_better_ = type_to_better_[type];
      }
      for (auto& i : id)
      {
        switchScores(i, counter);
      }
    }

    /// Looks at the first Hit of the given @p id and according to the given @p type ,
    /// deduces a fitting score and score direction to be switched to.
    /// Then tries to switch all hits.
    void switchToGeneralScoreType(ConsensusMap& cmap, ScoreType type, Size& counter, bool unassigned_peptides_too = true)
    {
      String new_type = "";
      for (const auto& f : cmap)
      {
        const auto& ids = f.getPeptideIdentifications();
        if (!ids.empty())
        {
          new_type = findScoreType(ids[0], type);
          if (new_type == ids[0].getScoreType())
          {
            return;
          }
          else
          {
            break;
          }
        }
      }

      if (new_type.empty())
      {
        String msg = "First encountered ID does not have the requested score type.";
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, msg);
      }

      if (new_type.hasSuffix("_score"))
      {
        new_score_type_ = new_type.chop(6);
      }
      else
      {
        new_score_type_ = new_type;
      }
      new_score_ = new_type;

      if (type != ScoreType::RAW && higher_better_ != type_to_better_[type])
      {
        OPENMS_LOG_WARN << "Requested non-raw score type does not match the expected score direction. Correcting!\n";
        higher_better_ = type_to_better_[type];
      }

      const auto switchScoresSingle = [&counter,this](PeptideIdentification& id){switchScores(id,counter);};
      cmap.applyFunctionOnPeptideIDs(switchScoresSingle, unassigned_peptides_too);
    }


    /// finds a certain score type in an ID and its metavalues if present, otherwise returns empty string
    template <typename IDType>
    String findScoreType(IDType& id, IDScoreSwitcherAlgorithm::ScoreType type)
    {
      const String& curr_score_type = id.getScoreType();
      const std::set<String>& possible_types = type_to_str_[type];
      if (possible_types.find(curr_score_type) != possible_types.end())
      {
        OPENMS_LOG_INFO << "Requested score type already set as main score: " + curr_score_type + "\n";
        return curr_score_type;
      }
      else
      {
        if (id.getHits().empty())
        {
          OPENMS_LOG_WARN << "Identification entry used to check for alternative score was empty.\n";
          return "";
        }
        const auto& hit = id.getHits()[0];
        for (const auto& poss_str : possible_types)
        {
          if (hit.metaValueExists(poss_str)) return poss_str;
          else if (hit.metaValueExists(poss_str + "_score")) return poss_str + "_score";
        }
        OPENMS_LOG_WARN << "Score of requested type not found in the UserParams of the checked ID object.\n";
        return "";
      }
    }

  private:
    void updateMembers_() override;

    /// relative tolerance for score comparisons:
    const double tolerance_ = 1e-6;

    /// will be set according to the algorithm parameters
    String new_score_, new_score_type_, old_score_;
    /// will be set according to the algorithm parameters
    bool higher_better_; // for the new scores, are higher ones better?

    /// a map from ScoreType to their names as used around OpenMS
    std::map<ScoreType, std::set<String>> type_to_str_ =
        {
            {ScoreType::RAW, {"XTandem", "OMSSA", "SEQUEST:xcorr", "Mascot", "mvh"}},
            //TODO find out reasonable raw scores for SES that provide evalues as main score or see below
            //TODO there is no test for spectraST idXML, so I dont know its score
            //TODO check if we should combine RAW and RAW_EVAL:
            // What if a SE does not have an e-value score (spectrast, OMSSA, crux/sequest, myrimatch),
            // then you need additional ifs/trys
            {ScoreType::RAW_EVAL, {"expect", "SpecEValue", "E-Value", "evalue", "MS:1002053", "MS:1002257"}},
            {ScoreType::PP, {"Posterior Probability"}},
            {ScoreType::PEP, {"Posterior Error Probability", "pep", "MS:1001493"}}, // TODO add CV terms
            {ScoreType::FDR, {"FDR", "fdr", "false discovery rate"}},
            {ScoreType::QVAL, {"q-value", "qvalue", "MS:1001491", "q-Value", "qval"}}
        };

    /// a map from ScoreType to their ordering
    std::map<ScoreType, bool> type_to_better_ =
        {
            {ScoreType::RAW, true}, //TODO this might actually not always be true
            {ScoreType::RAW_EVAL, false},
            {ScoreType::PP, true},
            {ScoreType::PEP, false},
            {ScoreType::FDR, false},
            {ScoreType::QVAL, false}
        };
  };
} // namespace OpenMS
