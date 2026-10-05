// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <set>

namespace OpenMS
{
/**
  @brief Protein inference over explicitly selected owning runs without editing their search evidence.

  The bridge materializes temporary peptide values, normalizes explicitly declared PEP/PP
  scores to posterior probabilities and calls BasicProteinInferenceAlgorithm once across
  all selected runs. Exact ordered candidate memberships and resulting assignments are
  recorded independently of the original matches. Repeated peptidoforms must have
    consistent qualified parent mappings across inputs; conflicting mappings are rejected.
    This is an in-memory algorithm.
  @ingroup Analysis_ID
*/
class OPENMS_DLLAPI IdentificationDataInference
{
public:
  enum class ProbabilityType
  {
    POSTERIOR_ERROR_PROBABILITY,
    POSTERIOR_PROBABILITY
  };
  struct OPENMS_DLLAPI Input
  {
    std::string run_uuid;
    IdentificationData::ScoreId score;
    ProbabilityType probability = ProbabilityType::POSTERIOR_ERROR_PROBABILITY;
  };
  /// Returns a result; the caller decides whether to attach or replace existing inference.
  static IdentificationData::InferenceResult
  infer(const IdentificationData& data, const std::vector<Input>& inputs, const std::string& identifier, const Param& parameters);
  static IdentificationData::InferenceResult infer(const IdentificationData& data, const std::vector<Input>& inputs, const std::string& identifier);
  /// Keep selected qualified proteins, drop any group that loses a member and preserve empty assignments.
  static void retainProteins(IdentificationData::InferenceResult& result, const std::set<IdentificationData::QualifiedAccession>& retained);
};
} // namespace OpenMS
