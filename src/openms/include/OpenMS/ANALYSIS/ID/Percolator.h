// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/PercolatorTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <memory>
#include <vector>

namespace OpenMS
{
  /**
    @brief In-process Percolator: semi-supervised target/decoy rescoring with q-values + PEPs.

    Wraps a stripped-down, vendored copy of Percolator (training + posterior
    estimation). Provides two entry points:

    - High-level rescore(std::vector<PeptideIdentification>&, ...): PSM-specific sugar
      that stamps scores back onto PeptideHits as meta values.

    - Low-level rescore(const RescoreInput&): domain-agnostic. Accepts any (feature
      matrix, target/decoy labels, optional grouping keys) triple. Suitable for PSM
      rescoring, transition rescoring, peak-group rescoring, or any other setting
      where a semi-supervised target/decoy classifier is meaningful. The internal
      vendored types assume "PSM-shaped" data; the wrapper synthesizes opaque row
      IDs to satisfy them, so callers never see Percolator-native types.

    @section preconditions Preconditions
    - At least ~100 decoys and enough discriminable targets to pass SanityCheck,
      else InvalidValue is thrown.
    - Features are continuous-valued doubles. Categorical features must be
      one-hot-encoded by the caller.
    - If the input contains rows that share structure (e.g., multiple transitions
      per precursor), cv_group_keys must partition them so related rows go into
      the same CV fold. Otherwise q-values will be optimistic.

    @section thread_safety Thread safety
    A single instance is NOT concurrent-safe; construct one per worker.

    @section reproducibility Reproducibility
    Reproducible given seed, thread count, and input ordering. Changing thread
    count can perturb results through FeatureMemoryPool allocation ordering.

    @see PSM rescoring example: PercolatorAdapter / ProSE.
    @see Transition rescoring example: OpenSwath layer (to be added when needed).
  */
  class OPENMS_DLLAPI Percolator : public DefaultParamHandler
  {
  public:
    Percolator();
    ~Percolator() override;

    Percolator(const Percolator&) = delete;
    Percolator& operator=(const Percolator&) = delete;

    /**
      @brief Rescore PSMs in place. Domain-specific convenience wrapper over rescore(RescoreInput).

      Each PeptideHit gets three new meta values:
      - "percolator_score" (SVM discriminant score)
      - "percolator_q_value" (PSM-level q-value)
      - "percolator_pep" (posterior error probability)

      @param peptide_ids Target + decoy PSMs, mixed. Mutated in place.
      @param feature_names Meta-value names on each PeptideHit to use as features.
                          Must be numeric. If empty, auto-discover from the first hit's
                          numeric meta values (excluding a blocklist of internal keys).
      @throws Exception::InvalidValue if sanity checks fail (too few decoys, no
              discriminative feature, etc.).
    */
    void rescore(std::vector<PeptideIdentification>& peptide_ids,
                 const StringList& feature_names = {});

    /**
      @brief Rescore a feature matrix domain-agnostically.

      @param input Feature matrix, target/decoy labels, optional CV grouping keys,
                   and feature names. See RescoreInput for row/column contract.
      @return Per-row SVM scores, PSM-level q-values, and PEPs, aligned 1:1 with
              input.features (no reordering).

      @note This method does not interpret the rows semantically. "PSM" terminology
            in the returned struct refers to the underlying Percolator algorithm,
            not to what your rows represent. For transition rescoring, the
            q_values/peps are transition-level.

      @note No peptide-level roll-up is performed. If you need "best row per
            unique entity" aggregation, do it above this call.

      @throws Exception::InvalidValue if sanity checks fail (too few decoys, no
              discriminative feature, malformed input dimensions).
    */
    RescoreOutput rescore(const RescoreInput& input);

    /**
      @brief Fill PIN-compatible optional fields on a RescoreInput.

      Populates `input.scan_numbers`, `input.spec_file_numbers`,
      `input.exp_masses`, and `input.calc_masses` from the given
      PeptideIdentifications, using the same derivation that
      PercolatorInfile::store would apply when writing a .pin file:

      - scan: parsed via SpectrumLookup::extractScanNumber from the
        PeptideIdentification's spectrum_reference (or `spectrum_id`
        meta value, or fallback to 1-based index).
      - spec_file: hashes `file_origin` + `id_merge_index` (same as the
        PIN SpecId prefix). Zero when single-file / unset.
      - exp_mass: `pid.getMZ()` (kept as m/z — Percolator doesn't convert
        to neutral for the sort hash).
      - calc_mass: from `hit.metaValueExists("CalcMass")` if present,
        else `hit.getSequence().getMonoWeight()`.

      Call this BEFORE rescore(input) if you need the in-process output to
      match running the external percolator binary through the .pin / .pout
      pipeline on the same inputs. Without it, rows get row_index as scan,
      which produces a different CV fold split and therefore different
      trained weights and final scores.

      The helper must be called with a PeptideIdentifications vector that
      parallels `input.features` (same ordering, one row per hit per pid).

      @param peptide_ids Source of PIN-equivalent metadata.
      @param flatten_hits If true, iterate all hits per PeptideIdentification
             (matches high-level rescore row ordering). If false, use only
             the first hit per pid.
      @param input Output: the four PIN-compat fields are written here.
    */
    static void fillPINCompatibleFields(
        const std::vector<PeptideIdentification>& peptide_ids,
        bool flatten_hits,
        RescoreInput& input);

  protected:
    void updateMembers_() override;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
}
