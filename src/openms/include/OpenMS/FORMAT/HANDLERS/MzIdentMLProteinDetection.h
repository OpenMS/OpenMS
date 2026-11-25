// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <string>

namespace OpenMS::Internal
{
  /**
    @brief Mapping between OpenMS protein inference results and the mzIdentML protein detection model

    Shared by the mzIdentML writer (MzIdentMLHandler) and reader (MzIdentMLDOMHandler), so that a store/load round
    trip restores the same score types and inference engine.

    The mzIdentML 1.2/1.3 model (specification section "Protein grouping") and its OpenMS counterpart:
    - @c ProteinAmbiguityGroup (PAG): a protein group, the unit of a protein identification. In OpenMS: an
      indistinguishable protein group (ProteinIdentification::getIndistinguishableProteins()); its proteins are
      written as "leading protein"s, @p accessions[0] as "group representative". The group score is a PAG cvParam
      (see groupScoreAccession()), "protein group passes threshold" says whether the group is reported as identified.
    - @c ProteinDetectionHypothesis (PDH): one protein of a group, with its score (see proteinScoreAccession()) and
      the peptide evidence (PSMs) that supports it.
    - @c ProteinDetectionProtocol: the inference engine (as @c AnalysisSoftware), the protein score type
      (SCORE_TYPE_PARAM, HIGHER_SCORE_BETTER_PARAM) and the cutoff applied to group scores, if one was recorded
      (Constants::UserParam::PROTEIN_GROUP_SCORE_THRESHOLD).
    - Cluster: groups that share peptides, given by the same "cluster identifier" on their PAGs. In OpenMS: a general
      protein group (ProteinIdentification::getProteinGroups()). There is no CV term for a cluster score, so it is
      stored in the userParam CLUSTER_SCORE_PARAM on the PAGs of the cluster.

    @note Do not use this class. It is only needed by the mzIdentML handlers.
  */
  class OPENMS_DLLAPI MzIdentMLProteinDetection
  {
  public:
    /// Name of the userParam on a @c ProteinAmbiguityGroup that holds the score of its cluster (general protein group)
    static const std::string CLUSTER_SCORE_PARAM;

    /// Name of the userParam in @c ProteinDetectionProtocol/AnalysisParams with the OpenMS protein score type (ProteinIdentification::getScoreType()).
    /// Protein and group scores without a PSI-MS term are written as userParams named after the score type.
    static const std::string SCORE_TYPE_PARAM;

    /// Name of the userParam in @c ProteinDetectionProtocol/AnalysisParams that says whether higher protein scores are better
    static const std::string HIGHER_SCORE_BETTER_PARAM;

    /**
      @brief Accession of the PSI-MS term for a protein (@c ProteinDetectionHypothesis) score

      @param[in] score_type Protein score type of the ProteinIdentification
      @param[in] inference_engine Inference engine that computed the score (ProteinIdentification::getInferenceEngine())
      @return The accession, or an empty string if there is no term for this score type (it is written as userParam then)
    */
    static std::string proteinScoreAccession(const std::string& score_type, const std::string& inference_engine);

    /// Accession of the PSI-MS term for a protein group (@c ProteinAmbiguityGroup) score of type @p score_type, or an empty string if there is none
    static std::string groupScoreAccession(const std::string& score_type);

    /// Accession of the PSI-MS term for a cutoff on protein group scores of type @p score_type (in @c ProteinDetectionProtocol/Threshold)
    static std::string groupThresholdAccession(const std::string& score_type);

    /**
      @brief OpenMS score type for a PSI-MS protein or protein group score term written by this mapping (or an equivalent one)

      @return "Posterior Probability" or "q-value", or an empty string for any other term
    */
    static std::string scoreTypeOfAccession(const std::string& accession);

    /// Name of the PSI-MS software term for an OpenMS inference engine name (ProteinIdentification::getInferenceEngine()), or an empty string if there is none
    static std::string softwareTermName(const std::string& inference_engine);

    /// OpenMS inference engine name for a software name read from mzIdentML (inverse of softwareTermName(); other names are returned unchanged)
    static std::string inferenceEngineOfSoftware(const std::string& software_name);
  };
} // namespace OpenMS::Internal
