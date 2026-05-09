// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/config.h>

#include <optional>
#include <vector>

namespace OpenMS
{
  /**
    @brief Inference level for OpenSWATH-derived confidence inference.

    See: 
    - Rosenberger, G., Bludau, I., Schmitt, U. et al. Statistical control of peptide and protein error rates in large-scale targeted data-independent acquisition analyses. Nat Methods 14, 921–927 (2017). https://doi.org/10.1038/nmeth.4398
    - Rosenberger, G., Liu, Y., Röst, H. et al. Inference and quantification of peptidoforms in large sample cohorts by SWATH-MS. Nat Biotechnol 35, 781–788 (2017). https://doi.org/10.1038/nbt.3908

    @ingroup TargetedQuantitation
  */
  enum class InferenceLevel
  {
    Peptidoform,
    Peptide,
    Protein,
    Gene
  };

  /**
    @brief Context for peptide-, protein-, and gene-level inference.

    @ingroup TargetedQuantitation
  */
  enum class InferenceContext
  {
    Global,
    ExperimentWide,
    RunSpecific
  };

  /// String conversion helper for inference levels.
  inline String toString(const InferenceLevel level)
  {
    switch (level)
    {
      case InferenceLevel::Peptidoform: return "peptidoform";
      case InferenceLevel::Peptide: return "peptide";
      case InferenceLevel::Protein: return "protein";
      case InferenceLevel::Gene: return "gene";
    }
    return "";
  }

  /// String conversion helper for inference contexts.
  inline String toString(const InferenceContext context)
  {
    switch (context)
    {
      case InferenceContext::Global: return "global";
      case InferenceContext::ExperimentWide: return "experiment-wide";
      case InferenceContext::RunSpecific: return "run-specific";
    }
    return "";
  }

  /**
    @brief Peak-group and precursor evidence for peptidoform inference.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI IPFPrecursorRow
  {
    Int64 feature_id = -1;
    double ms2_peakgroup_pep = 1.0;
    std::optional<double> ms1_precursor_pep;
    std::optional<double> ms2_precursor_pep;
  };

  /**
    @brief Precursor-layer posterior error probability used by IPF.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI IPFPrecursorProbabilityRow
  {
    Int64 feature_id = -1;
    double precursor_peakgroup_pep = 1.0;
  };

  /**
    @brief Transition-layer evidence for peptidoform inference.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI IPFTransitionRow
  {
    Int64 feature_id = -1;
    Int64 transition_id = -1;
    Int64 peptide_id = -1;
    double pep = 1.0;
    Int32 bmask = 0;
    Int32 num_peptidoforms = 0;
    std::optional<Int64> alignment_group_id;
  };

  /**
    @brief Alignment-group membership for cross-run IPF signal propagation.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI IPFAlignmentRow
  {
    Int64 alignment_group_id = -1;
    Int64 feature_id = -1;
  };

  /**
    @brief Final peptidoform inference result row.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI IPFResultRow
  {
    Int64 feature_id = -1;
    Int64 peptide_id = -1;
    double precursor_peakgroup_pep = 1.0;
    double qvalue = 1.0;
    double pep = 1.0;
  };

  /**
    @brief Compact level-context input row independent of file format.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI LevelContextInputRow
  {
    std::optional<Int64> run_id;
    String group_id;
    Int64 entity_id = -1;
    bool decoy = false;
    double score = 0.0;
    InferenceContext context = InferenceContext::RunSpecific;
  };

  /**
    @brief Final level-context result row.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI LevelContextResultRow
  {
    std::optional<Int64> run_id;
    Int64 entity_id = -1;
    InferenceContext context = InferenceContext::RunSpecific;
    double score = 0.0;
    double pvalue = 1.0;
    double qvalue = 1.0;
    double pep = 1.0;
  };
} // namespace OpenMS
