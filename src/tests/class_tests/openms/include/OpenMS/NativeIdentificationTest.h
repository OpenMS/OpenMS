// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

// Checks that code working on the identification data of feature and consensus maps does what it
// does on their peptide identifications. Header-only, like TestFileValidation.h.

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>

#include <string>
#include <type_traits>
#include <vector>

namespace OpenMS::Internal::ClassTest
{
  /// Move the peptide identifications of a map into its identification data (feature links)
  inline void toNative(FeatureMap& map)
  {
    IdentificationDataConverter::importFeatureIDs(map);
  }
  inline void toNative(ConsensusMap& map)
  {
    IdentificationDataConverter::importConsensusIDs(map);
  }
  /// Move the identification data of a map back into peptide identifications
  inline void toLegacy(FeatureMap& map)
  {
    IdentificationDataConverter::exportFeatureIDs(map);
  }
  inline void toLegacy(ConsensusMap& map)
  {
    IdentificationDataConverter::exportConsensusIDs(map);
  }

  namespace NativeIdentificationDetail
  {
    inline std::string metaDifference(const MetaInfoInterface& a, const MetaInfoInterface& b)
    {
      std::vector<std::string> keys_a, keys_b;
      a.getKeys(keys_a);
      b.getKeys(keys_b);
      std::string difference;
      for (const auto& key : keys_a)
        if (! b.metaValueExists(key)) difference += " missing meta value " + key;
        else if (a.getMetaValue(key) != b.getMetaValue(key))
          difference += " meta value " + key + ": " + a.getMetaValue(key).toString() + " vs. " + b.getMetaValue(key).toString();
      for (const auto& key : keys_b)
        if (! a.metaValueExists(key)) difference += " extra meta value " + key;
      return difference;
    }

    inline std::string peptideDifference(const PeptideIdentificationList& a, const PeptideIdentificationList& b)
    {
      if (a.size() != b.size()) return " " + std::to_string(a.size()) + " vs. " + std::to_string(b.size()) + " peptide identifications";
      for (Size i = 0; i < a.size(); ++i)
      {
        if (a[i] == b[i]) continue;
        std::string difference = " peptide identification " + std::to_string(i) + ":" + metaDifference(a[i], b[i]);
        if (a[i].getIdentifier() != b[i].getIdentifier()) difference += " identifier " + a[i].getIdentifier() + " vs. " + b[i].getIdentifier();
        if (a[i].getScoreType() != b[i].getScoreType()) difference += " score type";
        if (a[i].getRT() != b[i].getRT() || a[i].getMZ() != b[i].getMZ()) difference += " position";
        if (a[i].getHits().size() != b[i].getHits().size()) return difference + " number of hits";
        for (Size h = 0; h < a[i].getHits().size(); ++h)
        {
          const auto& p = a[i].getHits()[h];
          const auto& q = b[i].getHits()[h];
          if (p == q) continue;
          difference += " hit " + std::to_string(h) + ":" + metaDifference(p, q);
          if (p.getSequence() != q.getSequence()) difference += " sequence " + p.getSequence().toString() + " vs. " + q.getSequence().toString();
          if (p.getScore() != q.getScore()) difference += " score " + std::to_string(p.getScore()) + " vs. " + std::to_string(q.getScore());
          if (p.getRank() != q.getRank()) difference += " rank";
          if (p.getPeptideEvidences() != q.getPeptideEvidences()) difference += " evidence";
          break;
        }
        return difference;
      }
      return {};
    }

    template<class Map>
    bool hasPeptideIdentifications(const Map& map)
    {
      if (! map.getUnassignedPeptideIdentifications().empty() || ! map.getProteinIdentifications().empty()) return true;
      const auto any = [](const auto& self, const auto& feature) -> bool {
        if (! feature.getPeptideIdentifications().empty()) return true;
        if constexpr (std::is_same_v<std::remove_cvref_t<decltype(feature)>, Feature>)
          for (const auto& subordinate : feature.getSubordinates())
            if (self(self, subordinate)) return true;
        return false;
      };
      for (const auto& feature : map)
        if (any(any, feature)) return true;
      return false;
    }
  } // namespace NativeIdentificationDetail

  /// The first difference between the identifications (and then anything else) of two maps; empty if they are equal
  template<class Map>
  std::string mapDifference(const Map& expected, const Map& actual)
  {
    using namespace NativeIdentificationDetail;
    const auto& runs_a = expected.getProteinIdentifications();
    const auto& runs_b = actual.getProteinIdentifications();
    if (runs_a.size() != runs_b.size()) return std::to_string(runs_a.size()) + " vs. " + std::to_string(runs_b.size()) + " protein identification runs";
    for (Size i = 0; i < runs_a.size(); ++i)
    {
      if (runs_a[i] == runs_b[i]) continue;
      std::string difference = "protein identification run " + std::to_string(i) + ":" + metaDifference(runs_a[i], runs_b[i]);
      if (runs_a[i].getIdentifier() != runs_b[i].getIdentifier())
        difference += " identifier " + runs_a[i].getIdentifier() + " vs. " + runs_b[i].getIdentifier();
      if (runs_a[i].getSearchParameters() != runs_b[i].getSearchParameters()) difference += " search parameters";
      if (runs_a[i].getHits() != runs_b[i].getHits()) difference += " protein hits";
      return difference;
    }
    auto difference = peptideDifference(expected.getUnassignedPeptideIdentifications(), actual.getUnassignedPeptideIdentifications());
    if (! difference.empty()) return "unassigned:" + difference;
    if (expected.size() != actual.size()) return std::to_string(expected.size()) + " vs. " + std::to_string(actual.size()) + " features";
    for (Size i = 0; i < expected.size(); ++i)
    {
      difference = peptideDifference(expected[i].getPeptideIdentifications(), actual[i].getPeptideIdentifications());
      if (! difference.empty()) return "feature " + std::to_string(i) + ":" + difference;
    }
    if (! (expected.getIdentificationData() == actual.getIdentificationData())) return "identification data";
    if (! (expected == actual)) return "other than identifications";
    return {};
  }

  /**
    @brief Give the peptide identifications of a hand-built map what identification data needs: their search run and a score type

    The search run is the first protein identification run of @p map (a new one named "search" if it has none, or if it
    has no identifier); every peptide identification gets its identifier, and the score type @p score_type if it has none.
  */
  template<class Map>
  void addSearchRun(Map& map, const std::string& score_type = "score")
  {
    if (map.getProteinIdentifications().empty()) map.getProteinIdentifications().emplace_back();
    auto& run = map.getProteinIdentifications()[0];
    if (run.getIdentifier().empty()) run.setIdentifier("search");
    const std::string identifier = run.getIdentifier();
    map.applyFunctionOnPeptideIDs([&](PeptideIdentification& id) {
      id.setIdentifier(identifier);
      if (id.getScoreType().empty()) id.setScoreType(score_type);
    });
  }

  /**
    @brief Compare an operation on the peptide identifications of a map with the same operation on its identification data

    Runs @p operation(map, false) on a copy of @p input and @p operation(map, true) on a copy whose peptide identifications
    were moved into identification data; the latter must not produce peptide identifications. Its result is moved back
    into peptide identifications and compared with the former. @p operation returns the resulting map (of any type).

    @return The first difference (see mapDifference()), or an empty string if the results are equal
  */
  template<class Map, class Operation>
  std::string nativeDifference(const Map& input, Operation operation)
  {
    Map native_input = input;
    toNative(native_input);
    auto expected = operation(Map(input), false);
    auto actual = operation(std::move(native_input), true);
    if (NativeIdentificationDetail::hasPeptideIdentifications(actual)) return "the operation on identification data produced peptide identifications";
    toLegacy(actual);
    return mapDifference(expected, actual);
  }
} // namespace OpenMS::Internal::ClassTest
