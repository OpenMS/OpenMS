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
    @brief Utility class for score type handling in identification workflows.

    This class provides centralized handling of score types used in peptide and protein
    identification. It defines the hierarchy of score types (raw scores, E-values,
    posterior probabilities, etc.) and provides utility methods for score type
    conversion, comparison, and lookup.

    The score type hierarchy is:
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
    enum class Type
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

      This method determines if a given score name, typically derived from an identification
      object or meta value, matches a specified Type. It performs a case-insensitive comparison
      and optionally removes the "_score" suffix if present.

      @param[in] score_name The name of the score to check.
      @param[in] type The Type to compare against.
      @return True if the score name matches the given Type, false otherwise.
    */
    static bool isScoreType(const String& score_name, Type type);

    /**
      @brief Converts a string representation of a score type to a Type enum.

      This method attempts to map a given string, representing a score type, to the corresponding
      Type enum value. It handles various common representations of score types, including those
      with or without the "_score" suffix, and ignores case and special characters like '-', '_', and ' '.

      @param[in] score_type The string representation of the score type.
      @return The corresponding Type enum value.
      @throws Exception::MissingInformation If the provided score_type string does not match any known
                                            score type.
    */
    static Type toType(const String& score_type);

    /**
      @brief Determines whether a higher score is better for the given score type.

      @param[in] type The score type to check.
      @return True if a higher score is better, false otherwise.
    */
    static bool isHigherBetter(Type type);

    /**
      @brief Gets a vector of all score names that are used in OpenMS.

      @return A vector of all score names (e.g., "q-value", "ln(hyperscore)").
    */
    static std::vector<String> getScoreNames();

    /**
      @brief Gets the set of known names for a specific score type.

      @param[in] type The score type.
      @return A set of strings representing known names for this score type.
    */
    static const std::set<String>& getNamesForType(Type type);

    /**
      @brief Finds the score type for a given score name.

      Searches through all known score names to find a matching type.

      @param[in] name The score name to look up.
      @param[out] type Output parameter for the found score type.
      @return True if a matching type was found, false otherwise.
    */
    static bool findTypeByName(const String& name, Type& type);

  private:
    /// Initialize static maps (called once)
    static void initializeMaps_();

    /// a map from Type to their names as used around OpenMS
    static std::map<Type, std::set<String>> type_to_str_;

    /// a map from Type to their ordering (higher better or not)
    static std::map<Type, bool> type_to_better_;

    /// flag to track initialization
    static bool maps_initialized_;
  };

} // namespace OpenMS
