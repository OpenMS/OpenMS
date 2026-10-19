// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <map>
#include <string>
#include <vector>

namespace OpenMS
{
  class PeptideHit;

  /**
    @brief Measures the false discovery proportion (FDP) of identifications with an entrapment database.

    An entrapment database adds sequences to the searched targets that cannot be in the sample (for example
    shuffled copies of the targets, see DecoyDatabase -entrapment_ratio). The entrapment sequences are searched as
    targets, so the target-decoy FDR does not know about them. A target PSM that matches only entrapment
    sequences is false by construction, which gives an estimate of the true FDP among the accepted PSMs that
    is independent of the decoys:

      FDP = (N_E / (N_T + N_E)) * (1 + T/E)

    where N_E is the number of accepted entrapment-only PSMs, N_T + N_E all accepted (non-decoy) PSMs and T/E the
    ratio of searchable target to entrapment peptides (1 for a 1:1 shuffled entrapment database). This is the
    "combined" estimator of Wen et al., Nat. Methods (2025), which is an upper-bound-like estimate of the FDP.

    A PSM counts as entrapment only when @em all its proteins are entrapment proteins: a peptide that also occurs in
    a target protein is a target. Decoy PSMs (target_decoy = decoy) are ignored.

    Results can be stratified by a meta value of the peptide hits (for example a modification class), which shows
    whether the nominal FDR holds within each class.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI EntrapmentEvaluation
  {
  public:
    /// how a PSM is counted
    enum class Kind
    {
      TARGET,     ///< matches at least one target (non-entrapment) protein
      ENTRAPMENT, ///< matches only entrapment proteins
      DECOY       ///< decoy PSM (ignored)
    };

    /// result for one class at one q-value threshold
    struct OPENMS_DLLAPI Result
    {
      Size accepted_psms = 0;            ///< accepted non-decoy PSMs (targets + entrapment)
      Size entrapment_psms = 0;          ///< accepted entrapment-only PSMs
      Size accepted_peptides = 0;        ///< distinct accepted peptide sequences (with modifications)
      Size entrapment_peptides = 0;      ///< distinct accepted entrapment-only peptide sequences
      double te_ratio = 1.0;             ///< T/E used in the estimate

      /// combined FDP estimate at the PSM level (0 if nothing accepted)
      double psmFDP() const;
      /// combined FDP estimate at the peptide level (0 if nothing accepted)
      double peptideFDP() const;
    };

    /**
      @brief Classify a peptide hit by its protein accessions

      @param hit The hit; its peptide evidences carry the accessions
      @param entrapment_prefix Accession prefix of entrapment proteins (e.g. "ENTRAP_")
      @return DECOY if the hit is annotated as decoy, ENTRAPMENT if it has accessions and all start with
              @p entrapment_prefix, TARGET otherwise
    */
    static Kind classify(const PeptideHit& hit, const std::string& entrapment_prefix);

    /**
      @brief Estimate T/E, the ratio of distinct searchable target to entrapment peptides, from a database

      Proteins whose accession starts with @p decoy_string (as prefix or suffix, see @p decoy_is_prefix) are skipped.
      Peptides are counted after in-silico digestion with @p enzyme, within the length limits; a peptide that occurs
      in a target protein counts as target only.

      @return T/E, or 1 if the database has no entrapment peptides
    */
    static double estimateTERatio(const std::vector<FASTAFile::FASTAEntry>& db,
                                  const std::string& entrapment_prefix,
                                  const std::string& decoy_string,
                                  bool decoy_is_prefix,
                                  const std::string& enzyme = "Trypsin",
                                  Size missed_cleavages = 2,
                                  Size min_length = 7,
                                  Size max_length = 40);

    /**
      @brief Return the q-value of a hit

      Uses the hit's score when the score type of @p pep_id is a q-value (its name contains "q-value"), otherwise
      the meta values "q-value" or "MS:1001491" of the hit.

      @throws Exception::MissingInformation if no q-value is available
    */
    static double qValue(const PeptideIdentification& pep_id, const PeptideHit& hit);

    /**
      @brief Evaluate the top hit of each spectrum at a q-value threshold

      @param pep_ids Identifications (only the first hit of each is used; hits are expected to be sorted)
      @param q_threshold Accept hits with q <= this value
      @param entrapment_prefix Accession prefix of entrapment proteins
      @param te_ratio T/E (see estimateTERatio())
      @param class_meta_value If not empty, results are additionally reported per value of this hit meta value
             (hits without it are in class "NA")
      @return Results keyed by class; key "all" holds the overall result
    */
    static std::map<std::string, Result> evaluate(const PeptideIdentificationList& pep_ids,
                                                  double q_threshold,
                                                  const std::string& entrapment_prefix,
                                                  double te_ratio,
                                                  const std::string& class_meta_value = "");
  };
} // namespace OpenMS
