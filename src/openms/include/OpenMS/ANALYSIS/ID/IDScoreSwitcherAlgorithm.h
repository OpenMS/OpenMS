// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <algorithm>
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

    /// Checks if the given @p score_name is of ScoreType @p type
    bool isScoreType(const String& score_name, const ScoreType& type)
    {
      String chopped = score_name;
      if (chopped.hasSuffix("_score"))
      {
        chopped = chopped.chop(6);
      }
      const std::set<String>& possible_types = type_to_str_[type];
      return possible_types.find(chopped) != possible_types.end();
    }

    /// Gets a @p ScoreType enum from a given score name @p score_name
    static ScoreType getScoreType(String score_type)
    {
      if (score_type.hasSuffix("_score"))
      {
        score_type = score_type.chop(6);
      }
      score_type.toLower();     
      score_type.erase(std::remove_if(score_type.begin(), score_type.end(), 
                [](unsigned char c) { return c == '-' || c == '_' || c == ' '; }), 
                score_type.end());

      const std::map<String, ScoreType> s_to_type =
      {
        {"raw", ScoreType::RAW},
        {"rawevalue", ScoreType::RAW_EVAL},
        {"qvalue", ScoreType::QVAL},
        {"fdr", ScoreType::FDR},
        {"falsediscoveryrate", ScoreType::FDR},
        {"pep", ScoreType::PEP},
        {"posteriorerrorprobability", ScoreType::PEP},
        {"posteriorprobabilty", ScoreType::PP},
        {"pp", ScoreType::PP}
      };

      if (auto it = s_to_type.find(score_type); it != s_to_type.end())
      {
        return it->second;
      }
      else
      {
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, String("Unknown score type ") + score_type);
      }
    }

    /**
     * @brief Determines whether a higher score type is better given a ScoreType enum.
     * 
     * @param score_type The score type to check.
     * @return True if a higher score type is better, false otherwise.
     */
    bool isScoreTypeHigherBetter(ScoreType score_type)
    {
      return type_to_better_[score_type];
    }

    /*
      * @brief Gets a vector of all score names that are used in OpenMS.
      *
      * @return A vector of all score names that are used in OpenMS (e.g., "q-value", "ln(hyperscore)").
    */
    std::vector<String> getScoreTypeNames();

    /**
     * @brief Switches the main scores of all hits in an identification object based on the new scoring settings.
     *
     * This method iterates through all hits in the provided identification object and updates their main scores
     * according to the new scoring settings defined in the switcher class's parameter object. If the old and new
     * score types share the same name (e.g., "q-value"), the method safeguards the original scores by storing them
     * as meta values with a "~" appended to the old score type. This prevents overwriting the meta value of the new score.
     *
     * @tparam IDType The type of the identification object, which must support getHits(), getScoreType(),
     *                setScoreType(), and setHigherScoreBetter() methods, along with the ability to handle meta values.
     * @param[in,out] id An identification object containing hits whose scores are to be switched. The object will
     *                   be modified in place, with updated scores and score type.
     * @param[in,out] counter A reference to a Size variable that counts the number of hits processed.
     *
     * @throws Exception::MissingInformation If a required meta value (specified as the new score) is not found
     *                                       in any of the hits, indicating incomplete or incorrect score setup.
     *
     * @note The method assumes that the identification object's hits are properly initialized with all necessary
     *       meta values. It also relies on the tolerance_ value to determine significant differences between scores.     
     */ 
    template <typename IDType>
    void switchScores(IDType& id, Size& counter)
    {
      for (auto hit_it = id.getHits().begin();
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
          // TODO: find a better way to check if old score type is something different (even if it has same name)
          // This currently, is a workaround for e.g., having Percolator_qvalue as meta value and same q-value as main score (getScore()).
          // Note by jpfeuffer: The problem with this is, that this may add the old score to some of the hits if different, but not
          // all, in case one is by chance the same. I would be fine with this, if it was done in the beginning and checked
          // for every score.
          if (fabs((double(dv) - hit_it->getScore()) * 2.0 /
                   (double(dv) + hit_it->getScore())) > tolerance_)
          {          
            hit_it->setMetaValue(old_score_meta + "~", hit_it->getScore());
          }
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

    /**
     * @brief Switches the scoring type of identification objects to a general score type.
     *
     * This method iterates over a vector of identification objects and changes their scoring type
     * to a specified general score type. It first checks the score type of the first identification
     * object in the vector to determine the necessary conversion. If the first ID does not have the
     * requested score type, an exception is thrown. The method also adjusts the score direction
     * (higher_better_) based on the specified score type if it's different from the raw score.
     *
     * @tparam IDType The type of the identification objects contained in the vector. Must have
     *                getScoreType() and other relevant methods for score manipulation.
     * @param[in,out] id A vector of identification objects whose score types are to be switched.
     * @param[in] type The desired general score type to switch to. This could be an enum or similar
     *                 representing different scoring systems (e.g., RAW, LOG, etc.).
     * @param[in,out] counter A reference to a Size variable that may be used to count certain
     *                        operations or changes made by this method. The exact usage depends on
     *                        the implementation details and needs.
     *
     * @throws Exception::MissingInformation If the first identification object in the vector does not
     *                                       have the requested score type, indicating that the
     *                                       operation cannot proceed.
     *
     * @note The method assumes that if the first identification object has the correct score type,
     *       all subsequent objects in the vector also have the correct score type. This assumption
     *       might need validation depending on the use case.
     */    
    template<class IDType>
    void switchToGeneralScoreType(std::vector<IDType>& id, ScoreType type, Size& counter)
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

    
    /**
     * @brief Searches for a specified score type within an identification object and its meta values.
     *
     * This method attempts to find a given score type in the main score type of an identification object (`id`)
     * or within its hits' meta values. It first checks if the current main score type of `id` matches any of
     * the possible score types for the specified `type`. If not found, it iterates through the meta values of
     * the first hit in `id` looking for a match. If the score type or a related meta value is found, it is
     * returned as a `String`. Otherwise, an empty `String` is returned, indicating the score type is not present.
     *
     * @tparam IDType The type of the identification object, which must support getScoreType(), getHits(), and
     *                meta value operations.
     * @param[in] id The identification object to search for the score type. It is expected to have a main score
     *               type and possibly additional scores stored as meta values in its hits.
     * @param[in] type The `ScoreType` to search for, defined in `IDScoreSwitcherAlgorithm`. This type specifies
     *                 the score of interest.
     *
     * @return A String representing the found score type. If the score type is not found,
     *         an empty String is returned.
     *
     * @note This method logs an informational message if the requested score type is already set as the main score,
     *       a warning if the identification entry is empty, and another warning if the score type is not found in
     *       the UserParams of the checked ID object. 
     *       It only checks the first hit of the `id` for meta values.
     */    
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
          if (hit.metaValueExists(poss_str)) 
          {
            return poss_str;
          }
          else if (hit.metaValueExists(poss_str + "_score")) 
          {
            return poss_str + "_score";
          }
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
            //TODO introduce real meaningful score names for XTandem, Mascot etc. (e.g., hyperscore)
            {ScoreType::RAW, {"svm", "MS:1001492", "XTandem", "OMSSA", "SEQUEST:xcorr", "Mascot", "mvh", "hyperscore", "ln(hyperscore)"}},
            //TODO find out reasonable raw scores for SES that provide E-Values as main score or see below
            //TODO there is no test for spectraST idXML, so I don't know its score
            //TODO check if we should combine RAW and RAW_EVAL:
            // What if a SE does not have an e-value score (spectrast, OMSSA, crux/sequest, myrimatch),
            // then you need additional if's/try's
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
