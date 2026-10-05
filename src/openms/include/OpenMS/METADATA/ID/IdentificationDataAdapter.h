// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

namespace OpenMS
{
class FeatureMap;
class ConsensusMap;

/**
  @brief Explicit conversions between owning identifications and established peptide/protein values.

  Imports reject heterogeneous primary PSM score definitions; normalize them first. Export
  rejects unrepresentable information by default. ALLOW returns a loss report instead of
  silently discarding information. No inference freshness is inferred from current matches.
  @ingroup Metadata
*/
class OPENMS_DLLAPI IdentificationDataAdapter
{
public:
  enum class LossPolicy
  {
    STRICT,
    ALLOW
  };
  enum class MissingLinkPolicy
  {
    REJECT,
    PRUNE
  };
  struct OPENMS_DLLAPI ExportOptions
  {
    LossPolicy loss_policy = LossPolicy::STRICT;
    bool include_inference = true;
    std::optional<std::string> inference_result;
  };
  using QueryReference = IdentificationData::QueryReference;
  struct OPENMS_DLLAPI ImportResult
  {
    IdentificationData data;
    std::vector<QueryReference> queries; ///< In the original peptide-identification order.
  };
  struct OPENMS_DLLAPI LegacyResult
  {
    std::vector<ProteinIdentification> proteins;
    PeptideIdentificationList peptides;
    std::vector<QueryReference> queries; ///< Parallel to peptides.
    std::vector<std::string> losses;
  };
  struct OPENMS_DLLAPI FeatureAssociation
  {
    bool unassigned = false;
    std::vector<Size> feature_path; ///< Top-level index followed by subordinate indices.
    QueryReference query;
    std::vector<IdentificationData::MatchId> matches;
  };
  struct OPENMS_DLLAPI FeatureImportResult
  {
    IdentificationData data;
    std::vector<FeatureAssociation> associations;
  };

  /// Import with exact observation mappings, including identifications with no hits.
  static ImportResult importLegacy(const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides);
  /// Convenience import without the optional legacy-order mapping.
  static IdentificationData fromLegacy(const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides);
  /// Strict export with automatic selection only when an inference result is unambiguous.
  static LegacyResult toLegacy(const IdentificationData& data);
  static LegacyResult toLegacy(const IdentificationData& data, const ExportOptions& options);
  /// Explicit chemical materialization using the run's retained custom modification definitions.
  static PeptideHit materializePeptide(const IdentificationData::Run& run, const IdentificationData::Match& match, IdentificationData::ScoreId score);
  /// Collect live associations while retaining map values outside the identification model.
  static FeatureImportResult fromFeatureMap(const FeatureMap& map);
  static FeatureImportResult fromConsensusMap(const ConsensusMap& map);
  /// Prune deleted candidate links, or throw without changing any associations.
  static Size reconcileAssociations(const IdentificationData& data, std::vector<FeatureAssociation>& associations, MissingLinkPolicy policy);
  /// Replace identification annotations only; retain measured features, channels and intensities.
  static std::vector<std::string> applyToFeatureMap(const IdentificationData& data,
                                                    const std::vector<FeatureAssociation>& associations,
                                                    FeatureMap& map,
                                                    const ExportOptions& options,
                                                    MissingLinkPolicy policy);
  static std::vector<std::string> applyToConsensusMap(const IdentificationData& data,
                                                      const std::vector<FeatureAssociation>& associations,
                                                      ConsensusMap& map,
                                                      const ExportOptions& options,
                                                      MissingLinkPolicy policy);
};
} // namespace OpenMS
