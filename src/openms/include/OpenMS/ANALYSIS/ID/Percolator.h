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
    estimation). Provides several entry points grouped into two flavors.

    Fused train+score (original semantics; CV-merge-per-fold scoring):
    - rescore(std::vector<PeptideIdentification>&, ...) — PSM-specific sugar
      that stamps scores back onto PeptideHits as meta values.
    - rescore(const RescoreInput&) — domain-agnostic. Accepts any (feature
      matrix, target/decoy labels, optional grouping keys) triple. Suitable
      for PSM rescoring, transition rescoring, peak-group rescoring, or any
      other setting where a semi-supervised target/decoy classifier is
      meaningful. The internal vendored types assume "PSM-shaped" data; the
      wrapper synthesizes opaque row IDs to satisfy them, so callers never
      see Percolator-native types.

    Split train/score (for train-on-A / predict-on-B workflows; uses averaged
    weights across folds):
    - train(const RescoreInput&) → PercolatorModel — training only.
    - score(const RescoreInput&, const PercolatorModel&) → RescoreOutput —
      predict-only; applies raw weights × raw features + bias.
    - saveModel / loadModel — plain-text model persistence.

    rescore() and train()+score() produce scores on the same scale but not
    bitwise equal on the same input; see score() for why.

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
    Additionally, the vendored Percolator code relies on several process-wide
    statics (FeatureNames::numFeatures, SanityCheck::initDefaultDir* etc.)
    that are reset at the start of each rescore()/train()/score() call.
    Concurrent calls across *different instances* therefore also race on
    these globals — serialize at the call site if you need parallelism.

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

      @note Target-Decoy Competition (enabled by default via
            `post_processing_tdc`) deduplicates PSMs by (scan, expMass) —
            losing rows are dropped from the internal score set and come
            back with `score = 0.0`, `q_value = 1.0`, `pep = 1.0`. Callers
            that need to tell "weeded" rows from "genuinely bad" rows
            should turn TDC off, or populate scan_numbers/exp_masses so
            the deduplication key is meaningful.

      @throws Exception::InvalidValue if sanity checks fail (too few decoys, no
              discriminative feature, malformed input dimensions).
    */
    RescoreOutput rescore(const RescoreInput& input);

    /**
      @brief SVM weights trained in the last rescore()/train() call.

      Single averaged weight vector in raw feature space (+1 bias at the
      end). Same content as PercolatorModel::weights wrapped in a 1-element
      outer vector for compatibility with older call sites that expected
      per-fold weights. Useful for diagnostics or to write out a pseudo
      .weights file.

      Updated by rescore() and train(). score() does NOT touch this buffer —
      after a train()+score() sequence the weights reflect the train() call.
      Loading a model via loadModel() + calling score() on it leaves the
      buffer unchanged (empty if this instance never trained).

      @return Outer size = 1 (averaged); inner size = num_features + 1.
              Empty until rescore()/train() has been called on this instance.
    */
    const std::vector<std::vector<double>>& getSvmWeights() const;

    /**
      @brief Pi0 (null fraction) from the last rescore()/score() call.

      The target-null fraction estimated by the vendored
      `PosteriorEstimator::estimatePi0` on the merged, post-TDC Scores set
      right before PEP calculation. Matches the value the external percolator
      binary reports on stderr ("New pi_0 estimate on final list..." /
      "Selecting pi_0=..."). Useful for parity tests and diagnostics; not
      needed for production scoring.

      @return pi0 in [0, 1], or -1.0 if no rescore()/score() has run yet.
              Always 1.0 when the instance was configured with use_pi0=false.
    */
    double getPi0() const;

    /**
      @brief Train a Percolator model on feature rows. No scoring.

      Runs the full semi-supervised cross-validation training on @p input
      (all rows used — no subset sampling; do that at the call site if
      desired). The returned model holds raw-space averaged SVM weights
      suitable for passing to score(), possibly on a different RescoreInput.

      Side effect: mirrors the weights into getSvmWeights(). Resets per-call
      scratch state in the wrapper.

      @param input Training data. Target/decoy labels required; CV grouping
                   and PIN-compat fields used if supplied.
      @return PercolatorModel with feature_names, raw-space weights+bias,
              recorded normalizer type and seed.
      @throws Exception::InvalidValue on malformed input or training failure.
    */
    PercolatorModel train(const RescoreInput& input);

    /**
      @brief Score feature rows using a pre-trained model. No training.

      Applies @p model.weights to @p input.features as-is: raw weights × raw
      features + bias. Callers MUST NOT pre-normalize features — the
      normalization transform is already baked into the raw weights by
      Normalizer::unnormalizeweight() at train time. Doing it again
      double-counts the transform.

      The post-scoring pipeline (q-values, optional TDC, PEPs, SVM score
      rescaling) runs on the input's own target/decoy distribution — q-values
      and PEPs are always dataset-local, they do not carry over from the
      training dataset. This is the scientifically correct behavior.

      The returned `RescoreOutput.scores` are NOT raw dot products. Internally
      we compute
          raw_score(i) = Σ input.features[i][j] * model.weights[j]
                       + model.weights[n_features]
      and then `Scores::normalizeScores` rescales globally so the median
      decoy score maps to 0 and the score at `test_fdr` maps to 1:
          scores[i] = (raw_score(i) - fdrScore) / (fdrScore - medianDecoyScore)
      (Note: rescore()'s CV-merge path rescales per fold rather than once
      globally, so scores from rescore(X) vs. score(X, train(X)) live on the
      same scale but are not bitwise equal.) If you need raw dot-product
      values, compute them from model.weights directly.

      Target-Decoy Competition (enabled by default via `post_processing_tdc`)
      deduplicates rows by (scan, expMass); losing rows come back with the
      defaults `score = 0.0`, `q_value = 1.0`, `pep = 1.0`. See the
      rescore(RescoreInput) note for mitigation.

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

    /**
      @brief Serialize a PercolatorModel to a plain-text file.

      Writes a comment header line, then header keys `format_version`,
      `normalizer`, `seed`, `n_features`, `bias` (one per line, format
      `key: value`), then one `feature_name<TAB>weight` data row per
      feature. The bias sits in the header — feature names are therefore
      opaque strings with no reserved values. Intended to be human-readable
      and diff-friendly; not interoperable with the external percolator
      binary's multi-column .weights format.

      Weights are written at `std::numeric_limits<double>::max_digits10`
      precision so loadModel() round-trips losslessly.

      @throws Exception::UnableToCreateFile if @p filename cannot be opened.
      @throws Exception::InvalidValue if the model violates its internal
              invariant (weights.size() == feature_names.size() + 1).
    */
    static void saveModel(const PercolatorModel& model, const String& filename);

    /**
      @brief Deserialize a PercolatorModel written by saveModel().

      Strict reader: unknown header keys, duplicate keys, missing required
      keys, invalid normalizer values, or a declared `n_features` that
      does not match the feature row count are all rejected.

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
