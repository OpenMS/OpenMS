// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

// Builds and reads the identifications of feature and consensus maps in tests, as peptide and protein identifications
// apart from the map (MapPeptideIdentifications). Header-only, like TestFileValidation.h.

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/METADATA/ID/MapPeptideIdentifications.h>

#include <string>
#include <vector>

namespace OpenMS::Internal::ClassTest
{
  /// The identifications of a map as peptide and protein identifications (see IdentificationDataConverter::exportMapIdentifications())
  template<class Map>
  MapPeptideIdentifications identificationsOf(const Map& map)
  {
    return IdentificationDataConverter::exportMapIdentifications(map);
  }
  /// The peptide identifications of the feature at @p path of a map: its index, then the indices of its subordinates
  template<class Map>
  PeptideIdentificationList peptidesOf(const Map& map, const std::vector<Size>& path)
  {
    return identificationsOf(map).at(path);
  }
  /// The unassigned peptide identifications of a map
  template<class Map>
  PeptideIdentificationList unassignedOf(const Map& map)
  {
    return identificationsOf(map).unassigned;
  }
  /// The protein identifications of a map
  template<class Map>
  std::vector<ProteinIdentification> proteinsOf(const Map& map)
  {
    return identificationsOf(map).proteins;
  }
  /// Give a map the identifications @p identifications (see IdentificationDataConverter::importMapIdentifications())
  template<class Map>
  void setIdentifications(Map& map, const MapPeptideIdentifications& identifications)
  {
    IdentificationDataConverter::importMapIdentifications(map, identifications);
  }

  /**
    @brief Give hand-built peptide identifications what identification data needs: their search run and a score type

    The search run is the first protein identification run of @p identifications (a new one named "search" if there is
    none, or if it has no identifier); every peptide identification gets its identifier, and the score type
    @p score_type if it has none.
  */
  inline void addSearchRun(MapPeptideIdentifications& identifications, const std::string& score_type = "score")
  {
    if (identifications.proteins.empty()) identifications.proteins.emplace_back();
    auto& run = identifications.proteins[0];
    if (run.getIdentifier().empty()) run.setIdentifier("search");
    const std::string identifier = run.getIdentifier();
    const auto complete = [&](PeptideIdentificationList& peptides) {
      for (auto& id : peptides)
      {
        id.setIdentifier(identifier);
        if (id.getScoreType().empty()) id.setScoreType(score_type);
      }
    };
    for (auto& [path, peptides] : identifications.features) complete(peptides);
    complete(identifications.unassigned);
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

  } // namespace NativeIdentificationDetail

  /// The first difference between two sets of identifications of a map; empty if they are equal
  inline std::string identificationDifference(const MapPeptideIdentifications& expected, const MapPeptideIdentifications& actual)
  {
    using namespace NativeIdentificationDetail;
    const auto& runs_a = expected.proteins;
    const auto& runs_b = actual.proteins;
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
    auto difference = peptideDifference(expected.unassigned, actual.unassigned);
    if (! difference.empty()) return "unassigned:" + difference;
    if (expected.features.size() != actual.features.size())
      return std::to_string(expected.features.size()) + " vs. " + std::to_string(actual.features.size()) + " features with peptide identifications";
    for (auto a = expected.features.begin(), b = actual.features.begin(); a != expected.features.end(); ++a, ++b)
    {
      if (a->first != b->first) return "feature " + std::to_string(a->first[0]) + " vs. " + std::to_string(b->first[0]) + " with peptide identifications";
      difference = peptideDifference(a->second, b->second);
      if (! difference.empty()) return "feature " + std::to_string(a->first[0]) + ":" + difference;
    }
    return {};
  }

  /// The first difference between the identifications (and then anything else) of two maps; empty if they are equal
  template<class Map>
  std::string mapDifference(const Map& expected, const Map& actual)
  {
    auto difference = identificationDifference(identificationsOf(expected), identificationsOf(actual));
    if (! difference.empty()) return difference;
    if (expected.size() != actual.size()) return std::to_string(expected.size()) + " vs. " + std::to_string(actual.size()) + " features";
    if (! (expected.getIdentificationData() == actual.getIdentificationData())) return "identification data";
    if (! (expected == actual)) return "other than identifications";
    return {};
  }
} // namespace OpenMS::Internal::ClassTest
