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
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{
  class PeptideHit;

  /**
    @brief Measures the false discovery proportion (FDP) of identifications with an entrapment database.

    An entrapment database adds sequences to the searched targets that cannot be in the sample (for example
    shuffled copies of the targets, see DecoyDatabase -entrapment_ratio). The entrapment sequences are searched as
    targets, so the target-decoy FDR does not know about them. A discovery that matches only entrapment sequences is
    false by construction, which gives estimates of the true FDP among the accepted discoveries that are independent
    of the decoys. With N_T accepted original-target and N_E accepted entrapment discoveries, N = N_T + N_E, and r the
    ratio of searchable entrapment to target peptides (T/E = 1/r), following Wen et al., "Assessment of false discovery
    rate control in tandem mass spectrometry analysis using entrapment", Nat. Methods (2025):

      - lower bound:  FDP = N_E / N. A lower bound on the FDP; it can only show that FDR control fails.
      - combined:     FDP = N_E * (1 + T/E) / N. An estimated upper bound; below the nominal FDR it is evidence of
                      control.
      - paired:       FDP = (N_E + N_{E>=s>T} + 2 * N_{E>T>=s}) / N, for r = 1 with every entrapment peptide paired to
                      the target peptide it was made from. N_{E>=s>T} counts accepted entrapment peptides whose paired
                      target is not accepted, N_{E>T>=s} those whose paired target is accepted but ranks below or level
                      with them. Also an estimated upper bound, and a tighter one than the combined estimate. As in the
                      paper's FDRBench, it is computed at the peptide level only.

    Discoveries are counted at two levels: the top PSM of each spectrum, and peptides (distinct unmodified sequences,
    I and L taken as one residue). The paper evaluates peptides rather than PSMs. A peptide is accepted at its
    peptide-level q-value (meta value "peptide q-value", as written by FalseDiscoveryRate -peptide true) when its
    PSMs carry one; otherwise at the best q-value of its PSMs, which is then not a peptide-level FDR.

    A discovery counts as entrapment only when @em all its proteins are entrapment proteins: a peptide that also occurs
    in a target protein is a target. Given the target peptides of the database (targetPeptides()), a PSM whose sequence
    equals a target peptide with I and L taken as one residue also counts as target: I and L have the same mass, so the
    two cannot be told apart. Decoy PSMs (target_decoy = decoy) are ignored.

    Results can be stratified by a meta value of the peptide hits (for example a modification class), which shows
    whether the nominal FDR holds within each class. A peptide takes the class of its best PSM.

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

    /// what is counted as one discovery
    enum class Level
    {
      PSM,     ///< the top hit of each spectrum
      PEPTIDE  ///< distinct unmodified peptide sequences (I = L)
    };

    /// result for one class at one q-value threshold and level
    struct OPENMS_DLLAPI Result
    {
      Size accepted = 0;                       ///< N: accepted non-decoy discoveries (targets + entrapment)
      Size entrapment = 0;                     ///< N_E: accepted entrapment-only discoveries
      Size entrapment_target_not_accepted = 0; ///< N_{E>=s>T}: accepted entrapment peptides whose paired target is not accepted
      Size entrapment_above_target = 0;        ///< N_{E>T>=s}: accepted entrapment peptides whose paired target is accepted but does not rank above them
      bool paired = false;                     ///< whether the paired counts are available (peptide level, with peptide pairs)
      double te_ratio = 1.0;                   ///< T/E used in the combined estimate

      /// lower bound N_E / N (0 if nothing accepted)
      double lowerBoundFDP() const;
      /// combined estimate N_E * (1 + T/E) / N (0 if nothing accepted)
      double combinedFDP() const;
      /// paired estimate (N_E + N_{E>=s>T} + 2 * N_{E>T>=s}) / N (0 if nothing accepted, NaN if not available)
      double pairedFDP() const;
    };

    /// entrapment peptide -> the target peptide it was made from (both unmodified, I replaced by L)
    using PeptidePairs = std::unordered_map<std::string, std::string>;

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

      Proteins whose accession carries @p decoy_string are skipped: as a prefix, or (see @p decoy_is_prefix) as a suffix,
      optionally followed by "_<n>" as DecoyDatabase writes it with shuffle_decoy_ratio > 1.
      Peptides are counted after in-silico digestion with @p enzyme, within the length limits, with I and L taken as
      one residue; a peptide that occurs in a target protein counts as target only.

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
      @brief The distinct searchable peptides of the target proteins of a database, with every I replaced by L

      Digests like estimateTERatio(); decoy and entrapment proteins are skipped. Pass the result to evaluate().
    */
    static std::unordered_set<std::string> targetPeptides(const std::vector<FASTAFile::FASTAEntry>& db,
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
      @brief Pair every entrapment peptide with the target peptide it was made from

      Expects entrapment proteins as DecoyDatabase -entrapment_ratio 1 writes them: the accession is
      @p entrapment_prefix followed by the target's accession, and the sequence is the target's with the residues
      permuted within each peptide, so cleavage sites and peptide boundaries stay in place. Each target peptide (as
      digested like estimateTERatio()) is paired with the entrapment peptide at the same position. Pairs whose two
      peptides differ in composition (the entrapment protein was not made this way) are skipped and counted in
      @p mismatched; an entrapment peptide that equals its target peptide (I = L) is not paired. Decoy proteins
      (@p decoy_string, as for estimateTERatio()) are skipped.

      @return The pairs; empty if the database has no entrapment protein with a matching target (for example with
              entrapment_ratio > 1, whose accessions carry a "_<n>" suffix)
    */
    static PeptidePairs pairPeptides(const std::vector<FASTAFile::FASTAEntry>& db,
                                     const std::string& entrapment_prefix,
                                     const std::string& decoy_string,
                                     bool decoy_is_prefix,
                                     const std::string& enzyme = "Trypsin",
                                     Size missed_cleavages = 2,
                                     Size min_length = 7,
                                     Size max_length = 40,
                                     Size* mismatched = nullptr);

    /**
      @brief Evaluate the identifications at a q-value threshold

      @param pep_ids Identifications (only the first hit of each is used; hits are expected to be sorted)
      @param q_threshold Accept discoveries with q <= this value
      @param entrapment_prefix Accession prefix of entrapment proteins
      @param te_ratio T/E (see estimateTERatio())
      @param level Count PSMs or peptides
      @param class_meta_value If not empty, results are additionally reported per value of this hit meta value
             (hits without it are in class "NA")
      @param target_peptides If given (see targetPeptides()), an entrapment-only hit whose unmodified sequence, with I
             replaced by L, is one of these counts as target
      @param pairs If given (see pairPeptides()) and @p level is PEPTIDE, the paired counts are computed. The paired
             target of an entrapment peptide is the target peptide entry of @p pep_ids (all classes). An entrapment
             peptide without a pair counts as if its target ranked below it.
      @param score_meta_value If not empty, a hit meta value (higher is better) that ranks discoveries with equal
             q-values for the paired comparison; without it, or at equal scores, an entrapment peptide ranks level with
             its target, which counts against the FDP
      @return Results keyed by class; key "all" holds the overall result (a class named "all" is reported as "class:all")
    */
    static std::map<std::string, Result> evaluate(const PeptideIdentificationList& pep_ids,
                                                  double q_threshold,
                                                  const std::string& entrapment_prefix,
                                                  double te_ratio,
                                                  Level level = Level::PSM,
                                                  const std::string& class_meta_value = "",
                                                  const std::unordered_set<std::string>* target_peptides = nullptr,
                                                  const PeptidePairs* pairs = nullptr,
                                                  const std::string& score_meta_value = "");
  };
} // namespace OpenMS
