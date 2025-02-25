// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

  /**
    @brief This class is used to switch identification scores within identification or consensus feature maps.

    This class provides functionality to switch the main scoring type used in peptide or protein identification data.
    It supports switching between different score types, such as raw scores, E-values, posterior probabilities,
    posterior error probabilities, FDR, and q-values. The class also handles the direction of the score (whether a higher
    score is better) and can store the original scores as meta values to prevent data loss.

    The score switching process is configurable through parameters that specify the score types,
    as well as the desired score direction and how old scores are annotated in the meta information.

    The class can operate on individual identification objects or
    ConsensusMaps, updating the main scores of all hits based on the specified criteria.

  */
  class OPENMS_DLLAPI IDScoreSwitcherAlgorithm :
    public DefaultParamHandler
  {
  public:
    /// Default constructor. Initializes the parameter handler with default values.
    IDScoreSwitcherAlgorithm();

    /**
      @brief This is a rough hierarchy of possible score types in MS.

      In an ideal case, this should be reimplemented to follow
      ontology hierarchies as soon as e.g. MS-OBO is complete
      and we switched the Metavalues to CV terms.
    */
    enum class ScoreType
    {
      RAW,      ///< Raw score, e.g., search engine specific scores like hyperscore.
      RAW_EVAL, ///< Raw score with E-value, e.g., search engine specific scores like expect score.
      PP,       ///< Posterior probability.
      PEP,      ///< Posterior error probability.
      FDR,      ///< False discovery rate.
      QVAL,     ///< Q-value.
    };

    /**
      @brief Checks if the given score name corresponds to a specific score type.

      This method determines if a given score name, typically derived from an identification object or meta value,
      matches a specified ScoreType. It performs a case-insensitive comparison and optionally removes the "_score"
      suffix if present.

      @param score_name The name of the score to check.
      @param type The ScoreType to compare against.
      @return True if the score name matches the given ScoreType, false otherwise.
    */
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

    /**
      @brief Converts a string representation of a score type to a ScoreType enum.

      This static method attempts to map a given string, representing a score type, to the corresponding
      ScoreType enum value. It handles various common representations of score types, including those
      with or without the "_score" suffix, and ignores case and special characters like '-', '_', and ' '.

      @param score_type The string representation of the score type.
      @return The corresponding ScoreType enum value.
      @throws Exception::MissingInformation If the provided score_type string does not match any known
                                            score type.
    */
    static ScoreType toScoreTypeEnum(String score_type)
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
                                            OPENMS_PRETTY_FUNCTION, String("Unknown score type '") + score_type + "'.");
      }
    }

    /**
      @brief Determines whether a higher score type is better given a ScoreType enum.

      @param score_type The score type to check.
      @return True if a higher score type is better, false otherwise.
    */
    bool isScoreTypeHigherBetter(ScoreType score_type)
    {
      return type_to_better_[score_type];
    }

    /**
      @brief Gets a vector of all score names that are used in OpenMS.

      @return A vector of all score names that are used in OpenMS (e.g., "q-value", "ln(hyperscore)").
    */
    std::vector<String> getScoreNames();

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

      if (higher_better_ != type_to_better_[type])
      {
        OPENMS_LOG_WARN << "Requested score type does not match the expected score direction. Correcting!\n";
        higher_better_ = type_to_better_[type];
      }
      for (auto& i : id)
      {
        switchScores(i, counter);
      }
    }

    /**
      @brief Switches the score type of a ConsensusMap to a general score type.

      Looks at the first Hit of the given ConsensusMap and according to the given score type,
      deduces a fitting score and score direction to be switched to.
      Then tries to switch all hits.

      @param cmap The ConsensusMap containing peptide identifications whose scores need to be switched.
      @param type The desired general score type to switch to.
      @param counter A reference to a counter that will be incremented for each peptide identification processed.
      @param unassigned_peptides_too A boolean flag indicating whether to include unassigned peptides in the score switching process. Default is true.
      @throws Exception::MissingInformation If the first encountered ID does not have the requested score type.
    */
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

      if (higher_better_ != type_to_better_[type])
      {
        OPENMS_LOG_WARN << "Requested score type does not match the expected score direction. Correcting!\n";
        higher_better_ = type_to_better_[type];
      }

      const auto switchScoresSingle = [&counter,this](PeptideIdentification& id){switchScores(id,counter);};
      cmap.applyFunctionOnPeptideIDs(switchScoresSingle, unassigned_peptides_too);
    }

  /**
   @brief Determines the score type and orientation of the main score for a set of peptide identifications.

   This static method inspects a vector of PeptideIdentification objects to determine the overall score type and
   whether a higher score is considered better. It uses the first PeptideIdentification in the vector to make this
   determination, assuming that all identifications in the vector have the same score type and orientation.

   @param pep_ids The vector of PeptideIdentification objects to inspect.
   @param name Output parameter to store the determined overall score type.
   @param higher_better Output parameter to store whether a higher score is considered better.
   @param score_type Output parameter to store the determined score type.
   @note This method assumes that all PeptideIdentification objects in the input vector have the same score type and orientation.
  */
  void determineScoreNameOrientationAndType(
    const std::vector<PeptideIdentification>& pep_ids, 
    String& name, 
    bool& higher_better,
    ScoreType& score_type)
  {
    //TODO check all pep IDs? this assumes equality
    if (!pep_ids.empty())
    {
      name = pep_ids[0].getScoreType(); // The name of the score. Typically a name like "XTandem" or "Percolator_qvalue"
      higher_better = pep_ids[0].isHigherScoreBetter();
      
      // look up the score category ("RAW", "PEP", "q-value", etc.) for the given score name
      for (auto& [scoretype, names] : type_to_str_)
      {
        if (names.find(name) != names.end())
        {
          score_type = scoretype;
          OPENMS_LOG_INFO << "Found score type " << name << " to be of type " 
            << static_cast<std::underlying_type<ScoreType>::type>(scoretype) << std::endl;
          return;
        }
      }
    }
  }

  /**
   @brief Determines the score type and orientation of the main score in a ConsensusMap.

   This static method inspects a ConsensusMap to determine the overall score type and whether a higher score is
   considered better. It iterates through the ConsensusMap's features and uses the first PeptideIdentification found
   to determine the score type and orientation. If no assigned peptide identifications are found, it optionally
   considers unassigned peptide identifications.

   @param cmap The ConsensusMap to inspect.
   @param name Output parameter to store the determined overall score type.
   @param higher_better Output parameter to store whether a higher score is considered better.
   @param score_type Output parameter to store the determined score type.
   @param include_unassigned If true, unassigned peptide identifications are considered if no assigned ones are found. Default is true.
  */
  void determineScoreNameOrientationAndType(const ConsensusMap& cmap, 
    String& name,
    bool& higher_better,
    ScoreType& score_type,
    bool include_unassigned = true)
  {
    name = "";
    higher_better = true;

    // TODO: check all pep IDs? this assumes equality to first encountered
    for (const auto& cf : cmap)
    {
      const auto& pep_ids = cf.getPeptideIdentifications();
      if (!pep_ids.empty())
      {
        name = pep_ids[0].getScoreType();
        higher_better = pep_ids[0].isHigherScoreBetter();

        // look up the score category ("RAW", "PEP", "q-value", etc.) for the given score name
        for (auto& [scoretype, names] : type_to_str_)
        {
          if (names.find(name) != names.end())
          {
            score_type = scoretype;
            return;
          }
        }
      }
    }

    if (name.empty() && include_unassigned)
    {
      for (const auto& id : cmap.getUnassignedPeptideIdentifications())
      {
        name = id.getScoreType();
        higher_better = id.isHigherScoreBetter();

         // look up the score category ("RAW", "PEP", "q-value", etc.) for the given score name
        for (auto& [scoretype, names] : type_to_str_)
        {
          if (names.find(name) != names.end())
          {
            score_type = scoretype;
            return;
          }
        }        
        return;
      }
    }    
  }

    /**
     * @brief Switches the scores of peptide identifications in a ConsensusMap.
     *
     * This function iterates over all peptide identifications in the given ConsensusMap
     * and switches their scores using the switchScores function. It also increments the
     * provided counter for each peptide identification processed. Score names are
     * taken from the algorithm's parameters. If the requested score is already set as the
     * main score, the function returns without making any changes.
     *
     * @param cmap The ConsensusMap containing peptide identifications whose scores need to be switched.
     * @param counter A reference to a counter that will be incremented for each peptide identification processed.
     * @param unassigned_peptides_too A boolean flag indicating whether to include unassigned peptides in the score switching process. Default is true.
     */
    void switchScores(ConsensusMap& cmap, Size& counter, bool unassigned_peptides_too = true)
    {
      for (const auto& f : cmap)
      {
        const auto& ids = f.getPeptideIdentifications();
        if (!ids.empty())
        {
          if (new_score_ == ids[0].getScoreType()) // correct score or category already set
          {
            return;
          }
          else
          {
            break;
          }
        }
      }      
      const auto switchScoresSingle = [&counter,this](PeptideIdentification& id){switchScores(id,counter);};
      cmap.applyFunctionOnPeptideIDs(switchScoresSingle, unassigned_peptides_too);
    }
    
    /**
     * @brief Switches the scores of peptide identifications.
     *
     * This function iterates over all peptide identifications
     * and switches their scores using the switchScores function. It also increments the
     * provided counter for each peptide identification processed. Score names are
     * taken from the algorithm's parameters. If the requested score is already set as the
     * main score, the function returns without making any changes.
     *
     * @param pep_ids The peptide identifications whose scores need to be switched.
     * @param counter A reference to a counter that will be incremented for each peptide identification processed.     
     */
    void switchScores(std::vector<PeptideIdentification>& pep_ids, Size& counter)
    {
      if (pep_ids.empty()) return;

      if (new_score_ == pep_ids[0].getScoreType()) // correct score already set
      {
        return;
      }
      
      for (auto& id : pep_ids)
      {
        switchScores(id, counter);
      }
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

  /**
   * @brief Structure holding score switching information for IDScoreSwitcherAlgorithm.
   *
   * This structure contains both the original and requested score details, including
   * score names, their orientation (whether higher scores are better), and score types
   * before and after the switch. It also includes a flag to indicate if the main score
   * has been switched. Used to switch back to the original score if needed.   
   */
  struct IDSwitchResult
  {
    // the score name, orientation and type used before the switch
    String original_score_name; /// The name of the original score used before the switch.
    bool original_score_higher_better = true; /// whether a higher original score is better
    IDScoreSwitcherAlgorithm::ScoreType original_score_type = IDScoreSwitcherAlgorithm::ScoreType::RAW; /// the type of the original score
    // the score name, orientation and type used after the switch
    bool requested_score_higher_better = original_score_higher_better; /// whether a higher requested score is better
    IDScoreSwitcherAlgorithm::ScoreType requested_score_type = original_score_type; /// the type of the requested score
    String requested_score_name; // the search engine score name (e.g. "X!Tandem_score" or score category (e.g. "PEP")
    // wheter the main score was switched
    bool score_switched = false; /// flag indicating whether the main score was switched
  };

  /**
   * @brief Switches the score type of a ConsensusMap to the requested score type.
   *
   * This static method updates the scores within the provided ConsensusMap to the specified
   * score type. It determines the original score properties, checks if a switch is necessary
   * based on the requested score type, and performs the switch if required.
   *
   * @param cmap The ConsensusMap object whose score types are to be switched.
   * @param requested_score_type_as_string The desired score type as a string (e.g., "RAW", "PEP", "q-value").
   * @param include_unassigned Optional flag indicating whether to include unassigned IDs in the score switch. Defaults to true.
   *
   * @return An IDSwitchResult structure containing information about the score switch operation, including the original and requested score names, types, and whether a switch was performed.
   */
  static IDSwitchResult switchToScoreType(ConsensusMap& cmap, String requested_score_type_as_string, bool include_unassigned = true)
  {
    IDSwitchResult result;
    // fill in the original score name, orientation and type
    IDScoreSwitcherAlgorithm().determineScoreNameOrientationAndType(cmap, 
      result.original_score_name, 
      result.original_score_higher_better, 
      result.original_score_type,
      include_unassigned);

    // initalize with the assumption that the main score is the requested score
    result.requested_score_name = result.original_score_name; // the search engine score name (e.g. "X!Tandem_score" or score category (e.g. "PEP")
    result.requested_score_type = result.original_score_type;
    result.requested_score_higher_better = result.original_score_higher_better;          

    // no score type specified -> use main score
    if (requested_score_type_as_string.empty())
    {
      OPENMS_LOG_DEBUG << "No score type specified. Using main score." << std::endl;
      return result;
    }

    // ENUM for requested score type (e.g. "RAW", "PEP", "q-value")
    result.requested_score_type = IDScoreSwitcherAlgorithm::toScoreTypeEnum(requested_score_type_as_string);
    if (result.requested_score_type != result.original_score_type) // switch needed because we change type?
    { // user requests a different score type than the main score
      result.requested_score_higher_better = IDScoreSwitcherAlgorithm().isScoreTypeHigherBetter(result.requested_score_type);
      IDScoreSwitcherAlgorithm idsa;
      auto param = idsa.getDefaults();
      param.setValue("new_score", result.requested_score_name);
      param.setValue("new_score_orientation", result.requested_score_higher_better ? "higher_better" : "lower_better");
      param.setValue("proteins", "false");
      param.setValue("old_score", ""); // use default name generated for old score
      idsa.setParameters(param);            

      Size counter = 0; 
      idsa.switchToGeneralScoreType(cmap, result.requested_score_type, counter, include_unassigned);
      OPENMS_LOG_DEBUG << "Switched scores for " << counter << " IDs." << std::endl;
      result.score_switched = true;
    }

    // update after potential switch and read out actual score name
    IDScoreSwitcherAlgorithm().determineScoreNameOrientationAndType(cmap, 
      result.requested_score_name, 
      result.requested_score_higher_better, 
      result.requested_score_type,
      include_unassigned);

    return result;
  }

  /**
   * @brief Switches the score type of peptide identifications to the requested type.
   *
   * This static function modifies the provided vector of PeptideIdentification objects by switching
   * their main score to the specified type. If no score type is requested, the original main score
   * is retained. The function determines the original score's name, orientation, and type, and updates
   * these attributes based on the requested score type. If a different score type is requested,
   * it performs the switch and updates the relevant score information.
   *
   * @param pep_ids A vector of PeptideIdentification objects to be processed.
   * @param requested_score_type_as_string The desired score type as a string (e.g., "RAW", "PEP", "q-value").
   * @return IDSwitchResult A struct containing details about the original and requested score types,
   *                        whether a switch was performed, and the number of IDs updated.
   */
  static IDSwitchResult switchToScoreType(std::vector<PeptideIdentification>& pep_ids, String requested_score_type_as_string)
  {
    IDSwitchResult result;
    // fill in the original score name, orientation and type
    IDScoreSwitcherAlgorithm().determineScoreNameOrientationAndType(pep_ids, 
      result.original_score_name, 
      result.original_score_higher_better, 
      result.original_score_type
      );

    // initalize with the assumption that the main score is the requested score
    result.requested_score_name = result.original_score_name; // the search engine score name (e.g. "X!Tandem_score" or score category (e.g. "PEP")
    result.requested_score_type = result.original_score_type;
    result.requested_score_higher_better = result.original_score_higher_better;

    // no score type specified -> use main score
    if (requested_score_type_as_string.empty())
    {
      OPENMS_LOG_DEBUG << "No score type specified. Using main score." << std::endl;
      return result;
    }

    // ENUM for requested score type (e.g. "RAW", "PEP", "q-value")
    result.requested_score_type = IDScoreSwitcherAlgorithm::toScoreTypeEnum(requested_score_type_as_string);
    if (result.requested_score_type != result.original_score_type) // switch needed because we change type?
    { // user requests a different score type than the main score
      result.requested_score_higher_better = IDScoreSwitcherAlgorithm().isScoreTypeHigherBetter(result.requested_score_type);
      IDScoreSwitcherAlgorithm idsa;
      auto param = idsa.getDefaults();
      param.setValue("new_score", result.requested_score_name);
      param.setValue("new_score_orientation", result.requested_score_higher_better ? "higher_better" : "lower_better");
      param.setValue("proteins", "false");
      param.setValue("old_score", ""); // use default name generated for old score
      idsa.setParameters(param);            
      Size counter = 0;       
      idsa.switchToGeneralScoreType(pep_ids, result.requested_score_type, counter);
      OPENMS_LOG_DEBUG << "Switched scores for " << counter << " IDs." << std::endl;

      result.score_switched = true;
    }

    // update after potential switch and read out actual score name
    IDScoreSwitcherAlgorithm().determineScoreNameOrientationAndType(pep_ids, 
      result.requested_score_name, 
      result.requested_score_higher_better, 
      result.requested_score_type
      );

    return result;
  }

  /**
   * @brief Reverts the score type of a ConsensusMap to its original type based on the provided IDSwitchResult.
   *
   * This function checks if the scores have been switched and, if so, it switches them back to the original score type.
   * It updates the ConsensusMap by switching the scores, optionally including unassigned PSMs.
   *
   * @param cmap The ConsensusMap object whose scores will be modified.
   * @param isr The IDSwitchResult containing information about the score switch.
   * @param include_unassigned A boolean flag indicating whether to include unassigned PSMs in the score switching process. Defaults to true.
   */
  static void switchBackScoreType(ConsensusMap& cmap, IDSwitchResult isr, bool include_unassigned = true)
  {
    if (isr.score_switched)
    {
      // switch back to original score
      IDScoreSwitcherAlgorithm idsa;
      auto param = idsa.getDefaults();
      param.setValue("new_score", isr.original_score_name);
      param.setValue("new_score_orientation", isr.original_score_higher_better ? "higher_better" : "lower_better");
      param.setValue("proteins", "false");
      param.setValue("old_score", ""); // use default name generated for old score
      idsa.setParameters(param);
      Size counter = 0;
      idsa.switchScores(cmap, counter, include_unassigned);
      OPENMS_LOG_DEBUG << "Switched scores back for " << counter << " PSMs." << std::endl;
    }
  }

  /**
   * @brief Reverts the scoring type of peptide identifications to their original scores.
   *
   * This function checks if the scores have been switched. If so, it restores the original scoring parameters
   * using the provided IDSwitchResult. It updates the peptide identifications accordingly and logs the number
   * of PSMs (Peptide-Spectrum Matches) that were reverted.
   *
   * @param pep_ids A vector of PeptideIdentification objects to be updated.
   * @param isr An IDSwitchResult object containing information about the score switch state and original score details.
   */
  static void switchBackScoreType(std::vector<PeptideIdentification>& pep_ids, IDSwitchResult isr)
  {
    if (isr.score_switched)
    {
      // switch back to original score
      IDScoreSwitcherAlgorithm idsa;
      auto param = idsa.getDefaults();
      param.setValue("new_score", isr.original_score_name);
      param.setValue("new_score_orientation", isr.original_score_higher_better ? "higher_better" : "lower_better");
      param.setValue("proteins", "false");
      param.setValue("old_score", ""); // use default name generated for old score
      idsa.setParameters(param);
      Size counter = 0;
      idsa.switchScores(pep_ids, counter);
      OPENMS_LOG_DEBUG << "Switched scores back for " << counter << " PSMs." << std::endl;
    }
  }

  private:

    void updateMembers_() override; ///< documented in base class

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
