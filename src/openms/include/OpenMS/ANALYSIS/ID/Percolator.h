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

    Wraps a vendored subset of Percolator (training and posterior estimation).
    The public API is grouped into two sets of entry points.

    Combined training and scoring (original Percolator semantics; per-fold
    score normalization on the merged set):
    - rescore(std::vector<PeptideIdentification>&, ...) — overload for
      PSM-shaped data; writes scores back to each PeptideHit as meta values.
    - rescore(const RescoreInput&) — domain-agnostic. Accepts a feature
      matrix, target/decoy labels, optional CV grouping keys, and feature
      names. Applicable to PSM rescoring, transition rescoring, peak-group
      rescoring, or any other setting where a semi-supervised target/decoy
      classifier is appropriate. The vendored implementation requires
      PSM-shaped records internally; the wrapper attaches synthetic row
      identifiers to satisfy that requirement and these never surface in
      the public API.

    Separate training and scoring (for train-on-A / predict-on-B workflows;
    uses fold-averaged weights):
    - train(const RescoreInput&) → PercolatorModel — training only.
    - score(const RescoreInput&, const PercolatorModel&) → RescoreOutput —
      prediction only; computes raw weights × raw features + bias.
    - saveModel / loadModel — plain-text model persistence.

    rescore() and train()+score() produce scores on the same scale but are
    not bitwise equal on the same input; see score() for the rationale.

    @section preconditions Preconditions
    - At least ~100 decoys and enough discriminable targets to pass SanityCheck,
      else InvalidValue is thrown.
    - Features are continuous-valued doubles. Categorical features must be
      one-hot-encoded by the caller.
    - If the input contains rows that share structure (e.g., multiple transitions
      per precursor), cv_group_keys must partition them so related rows go into
      the same CV fold. Otherwise q-values will be optimistic.

    @section thread_safety Thread safety
    A single instance is not concurrent-safe; construct one instance per
    worker. The vendored Percolator code additionally relies on several
    process-wide statics (FeatureNames::numFeatures,
    SanityCheck::initDefaultDir*, etc.) that are reset at the start of
    each rescore() / train() / score() call. Concurrent calls across
    *different instances* therefore also race on these globals; serialize
    at the call site if parallelism is required.

    @section reproducibility Reproducibility
    Results are reproducible given the same seed, thread count, and input
    ordering. Changing the thread count can perturb results because of
    FeatureMemoryPool allocation ordering.

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

      @note Target-Decoy Competition (enabled by default via
            `post_processing_tdc`) deduplicates PSMs by (scan, expMass).
            Rows that lose competition are dropped from the internal score
            set and returned with `score = 0.0`, `q_value = 1.0`,
            `pep = 1.0`. To distinguish a TDC-eliminated row from a row
            that genuinely scored zero, either disable TDC or populate
            scan_numbers / exp_masses so that the deduplication key is
            informative.

      @throws Exception::InvalidValue if sanity checks fail (too few decoys, no
              discriminative feature, malformed input dimensions).
    */
    RescoreOutput rescore(const RescoreInput& input);

    /**
      @brief SVM weights trained in the last rescore()/train() call.

      A single fold-averaged weight vector in raw feature space, with the
      bias appended as the final element. Identical to
      PercolatorModel::weights, wrapped in a one-element outer vector to
      preserve the signature expected by older call sites that consumed
      per-fold weights. Intended for diagnostics and for writing out a
      Percolator-compatible .weights file.

      Updated by rescore() and train(). score() does not modify this
      buffer; after a train()+score() sequence the contents reflect the
      train() call. Loading a model via loadModel() and then calling
      score() leaves the buffer unchanged (empty if this instance has
      never trained).

      @return Outer size = 1 (averaged); inner size = num_features + 1.
              Empty until rescore()/train() has been called on this instance.
    */
    const std::vector<std::vector<double>>& getSvmWeights() const;

    /**
      @brief Pi0 (null fraction) from the last rescore()/score() call.

      The target-null fraction estimated by the vendored
      `PosteriorEstimator::estimatePi0` on the merged, post-TDC Scores set
      immediately before PEP calculation. Equivalent to the value reported
      on stderr by the external percolator binary ("New pi_0 estimate on
      final list..." / "Selecting pi_0=..."). Intended for parity testing
      and diagnostics; not required for production scoring.

      @return pi0 in [0, 1], or -1.0 if no rescore()/score() has run yet
              on this instance, or after a train() call (train() resets
              last_pi0 since it doesn't run scoring). Always 1.0 when the
              instance was configured with use_pi0=false.
    */
    double getPi0() const;

    /**
      @brief Train a Percolator model on feature rows. No scoring.

      Runs the full semi-supervised cross-validation training on @p input.
      All rows participate in training; subset sampling, if desired, must
      be performed by the caller. The returned model holds fold-averaged
      SVM weights in raw feature space and can be passed to score(),
      including on a different RescoreInput.

      Side effects: the weights are mirrored into getSvmWeights() and
      per-call scratch state in the wrapper is reset.

      @param input Training data. Target/decoy labels required; CV grouping
                   and PIN-compat fields used if supplied.
      @return PercolatorModel with feature_names, raw-space weights+bias,
              recorded normalizer type and seed.
      @throws Exception::InvalidValue on malformed input or training failure.
    */
    PercolatorModel train(const RescoreInput& input);

    /**
      @brief Score feature rows using a pre-trained model. No training.

      Applies @p model.weights to @p input.features unchanged: raw weights
      × raw features + bias. Callers must not pre-normalize the features.
      The normalization transform is already folded into the raw weights
      by `Normalizer::unnormalizeweight()` at training time; reapplying it
      would double-count the transform.

      The post-scoring pipeline (q-values, optional TDC, PEPs, SVM score
      rescaling) operates on the input's own target/decoy distribution.
      q-values and PEPs are therefore always evaluated against the
      scoring dataset, not the training dataset.

      The returned `RescoreOutput.scores` are not raw dot products.
      Internally
          raw_score(i) = Σ input.features[i][j] * model.weights[j]
                       + model.weights[n_features]
      `Scores::normalizeScores` then rescales the entire score vector so
      that the median decoy score maps to 0 and the score at `test_fdr`
      maps to 1:
          scores[i] = (raw_score(i) - fdrScore) / (fdrScore - medianDecoyScore)
      The rescore() CV-merge path rescales per fold instead of once
      globally, so scores from rescore(X) and from score(X, train(X))
      share a scale but are not bitwise equal. If raw dot-product values
      are required, compute them directly from model.weights.

      Target-Decoy Competition (enabled by default via `post_processing_tdc`)
      deduplicates rows by (scan, expMass). Rows that lose competition are
      returned with the defaults `score = 0.0`, `q_value = 1.0`,
      `pep = 1.0`. See the rescore(RescoreInput) note for mitigation.

      @param input Scoring data. Target/decoy labels required for q-value
                   and PEP computation. feature_names must be populated.
      @param model Model produced by train() (possibly from another
                   Percolator instance, or loaded from disk). feature_names
                   must be populated and match @p input.feature_names
                   positionally; weight count must equal n_features+1.
      @return Per-row SVM scores (rescaled per above), q-values, PEPs,
              aligned 1:1 with @p input.features.
      @throws Exception::InvalidValue on malformed input, model mismatch,
              empty feature_names on either side, or pathological score
              distributions that fail Percolator's sanity checks.
    */
    RescoreOutput score(const RescoreInput& input, const PercolatorModel& model);

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
        else `hit.getSequence().getMZ(hit.getCharge())` (m/z, matching
        the PercolatorInfile::store fallback).

      This helper must be called before rescore(input) when the
      in-process output is required to match running the external
      percolator binary through the .pin / .pout pipeline on the same
      inputs. Without it, the row index is used in place of the scan
      number, producing a different CV fold split and consequently
      different trained weights and final scores.

      The PeptideIdentifications vector passed in must parallel
      `input.features` exactly: same ordering, one row per hit per
      PeptideIdentification.

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

    /**
      @brief Serialize a PercolatorModel to a plain-text file.

      Writes a comment header line, then the header keys `format_version`,
      `normalizer`, `seed`, `n_features`, and `bias` (one per line, in
      `key: value` form), followed by one `feature_name<TAB>weight` data
      row per feature. The bias is stored in the header rather than as a
      data row, so feature names are opaque strings with no reserved
      values. The format is intended to be human-readable and
      diff-friendly; it is not interoperable with the external percolator
      binary's multi-column .weights format.

      Weights are written at `std::numeric_limits<double>::max_digits10`
      precision so that loadModel() round-trips losslessly.

      @throws Exception::UnableToCreateFile if @p filename cannot be opened.
      @throws Exception::InvalidValue if the model violates its internal
              invariant (weights.size() == feature_names.size() + 1).
    */
    static void saveModel(const PercolatorModel& model, const String& filename);

    /**
      @brief Deserialize a PercolatorModel written by saveModel().

      The reader is strict: unknown header keys, duplicate keys, missing
      required keys, invalid normalizer values, and a declared
      `n_features` that does not match the actual feature-row count are
      all rejected.

      @throws Exception::FileNotFound if @p filename does not exist.
      @throws Exception::ParseError on any format violation (unknown or
              duplicate header key, missing required field, unsupported
              format_version, invalid normalizer, or row/count mismatch).
    */
    static PercolatorModel loadModel(const String& filename);

  protected:
    void updateMembers_() override;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
}
