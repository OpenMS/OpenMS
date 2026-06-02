// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/OpenMSConfig.h>

#include <vector>
#include <set>
#include <map>
#include <algorithm>

namespace OpenMS
{

  /**
    @brief Utility class for score type handling in identification and quantification workflows.

    This class provides centralized handling of score types used in peptide/protein
    identification, quantification, and PTM localization. It defines the hierarchy of
    score types and provides utility methods for score type conversion, comparison, and lookup.

    The identification score type hierarchy is:
    - RAW: Raw search engine scores (e.g., XTandem hyperscore, Mascot score)
    - RAW_EVAL: E-value based scores (e.g., expect score)
    - PP: Posterior probability
    - PEP: Posterior error probability
    - FDR: False discovery rate
    - QVAL: Q-value

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI Scores
  {
  public:
    /**
      @brief Hierarchy of possible score types in MS identification.

      In an ideal case, this should be reimplemented to follow
      ontology hierarchies as soon as e.g. MS-OBO is complete
      and we switched the Metavalues to CV terms.
    */
    enum class IDType
    {
      RAW,      ///< Raw score, e.g., search engine specific scores like hyperscore.
      RAW_EVAL, ///< Raw score with E-value, e.g., search engine specific scores like expect score.
      PP,       ///< Posterior probability.
      PEP,      ///< Posterior error probability.
      FDR,      ///< False discovery rate.
      QVAL,     ///< Q-value.
    };

    /**
      @brief Structure to hold score detection results for any IDType.
    */
    struct ScoreSearchResult
    {
      bool is_main_score_type = false;  ///< True if the main score is already of the requested score type
      String score_name;                ///< Name of score to use (main score name if is_main_score_type=true, meta value name if found in meta values, empty if not found anywhere)
    };

    /**
      @brief Checks if the given score name corresponds to a specific ID score type.

      This method determines if a given score name, typically derived from an identification
      object or meta value, matches a specified IDType. It performs a case-insensitive comparison
      and optionally removes the "_score" suffix if present.

      @param[in] score_name The name of the score to check.
      @param[in] type The IDType to compare against.
      @return True if the score name matches the given IDType, false otherwise.
    */
    static bool isScoreType(const String& score_name, IDType type);

    /**
      @brief Converts a string representation of an ID score type to an IDType enum.

      This method attempts to map a given string, representing a score type, to the corresponding
      IDType enum value. It handles various common representations of score types, including those
      with or without the "_score" suffix, and ignores case and special characters like '-', '_', and ' '.

      @param[in] score_type The string representation of the score type.
      @return The corresponding IDType enum value.
      @throws Exception::MissingInformation If the provided score_type string does not match any known
                                            score type.
    */
    static IDType parseIDType(const String& score_type);

    /**
      @brief Determines whether a higher score is better for the given ID score type.

      @param[in] type The ID score type to check.
      @return True if a higher score is better, false otherwise.
    */
    static bool isHigherBetter(IDType type);

    /**
      @brief Gets a vector of all ID score names that are used in OpenMS.

      @return A vector of all ID score names (e.g., "q-value", "ln(hyperscore)").
    */
    static std::vector<String> getAllIDScoreNames();

    /**
      @brief Gets the set of known names for a specific ID score type.

      @param[in] type The ID score type.
      @return A set of strings representing known names for this score type.
    */
    static const std::set<String>& getIDNamesForType(IDType type);

    /**
      @brief Finds the ID score type for a given score name.

      Searches through all known score names to find a matching type.

      @param[in] name The score name to look up.
      @param[out] type Output parameter for the found score type.
      @return True if a matching type was found, false otherwise.
    */
    static bool findIDTypeByName(const String& name, IDType& type);

    /**
      @brief Normalizes a score name by removing the "_score" suffix if present.

      This is useful when checking if a score name matches a known score type,
      as OpenMS conventions allow both "q-value" and "q-value_score" forms.

      @param[in] score_name The score name to normalize.
      @return The normalized score name (without "_score" suffix).
    */
    static String normalizeScoreName(const String& score_name);

    /**
      @brief Checks if a score name is a known score type (after normalization).

      This method normalizes the score name and checks if it matches any known
      score type in the registry. Unlike isScoreType(), this doesn't require
      specifying which IDType to check - it checks all of them.

      @param[in] score_name The score name to check.
      @return True if the normalized name matches any known score type.
    */
    static bool isKnownScoreType(const String& score_name);

    /**
       @brief Searches for a general score type (e.g. PEP, QVAL) in an identification data structure.

       Returns a ScoreSearchResult for any identification type,
       checking if the main score of an identification object is already of the requested score type,
       and if not, searches for scores of that type in the meta values of the first hit.

       @tparam IdentificationType The type of the identification object (e.g., PeptideIdentification, ProteinIdentification)
       @param[in] id The identification object to analyze for scores
       @param[in] score_type The IDType to search for (e.g., IDType::PEP, IDType::QVAL, etc.)
       @return ScoreSearchResult containing whether main score is of the requested type and its name.

       @note This method only checks the first hit for meta values, similar to other methods in this class.
       @note Returns empty score_name if no score of the requested type is found.
     */
    template <typename IdentificationType>
    static ScoreSearchResult findScoreType(const IdentificationType& id, IDType score_type)
    {
      ScoreSearchResult result;

      // First check if main score is already of the requested score type using existing infrastructure
      const String& main_score_type = id.getScoreType();
      result.is_main_score_type = isScoreType(main_score_type, score_type);

      if (result.is_main_score_type)
      {
        // Main score is of the requested type, so return the main score name
        result.score_name = main_score_type;
      }
      else if (!id.getHits().empty())
      {
        // Main score is not of the requested type, look for it in meta values
        const auto& first_hit = id.getHits()[0];
        const std::set<String>& score_types = getIDNamesForType(score_type);

        // Search for scores of the requested type in meta values using the existing score type collection
        for (const String& score_name : score_types)
        {
          if (first_hit.metaValueExists(score_name))
          {
            result.score_name = score_name;
            break;
          }
          // Also check for "_score" suffix variant
          String score_name_with_suffix = score_name + "_score";
          if (first_hit.metaValueExists(score_name_with_suffix))
          {
            result.score_name = score_name_with_suffix;
            break;
          }
        }
      }
      // If neither main score nor meta values contain the requested type, score_name remains empty

      return result;
    }

  private:
    /// Holds the static score type lookup maps (thread-safe via C++11 function-local static)
    struct Maps_
    {
      std::map<IDType, std::set<String>> type_to_str;
      std::map<IDType, bool> type_to_better;
    };

    /// Returns the singleton Maps_ instance (thread-safe initialization guaranteed by C++11)
    static const Maps_& getMaps_();
  };

} // namespace OpenMS
