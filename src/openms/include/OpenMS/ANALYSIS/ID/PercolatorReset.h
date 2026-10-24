// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/TypeAliases.h>

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Percolator-RESET: semi-supervised rescoring with FDR control at the peptide (stem) level.

    Implements the single-decoy RESET procedure of Freestone, Käll, Noble and Keich (J Proteome Res 2025,
    "How to train a postprocessor ... while maintaining control of the FDR"; theory: Freestone, Noble and Keich,
    arXiv:2411.15771), as Percolator 3.07+ ('--reset-algorithm') and the authors' Python implementation do:

    1. PSM-level competition: per spectrum (Input::spectrum) the hit with the best search engine score.
    2. Stem reduction: per stem (unmodified sequence; Input::stems) the best of these PSMs by search engine score; it
       carries that PSM's features and label (label-blind: a sequence with target and decoy PSMs is one stem). Every
       stem is considered ("PSM-only" protocol of the paper: ProSE's decoys are reversed proteins, which give no fixed
       target-decoy peptide pairing; a pseudo-reverse pairing was tried and removed, review round 2).
    3. Coin: every decoy stem is a training decoy or an estimating decoy, by an independent Bernoulli(1/2) coin:
       bit 0 of a 64-bit hash of (salt, stem). Pseudo-targets = targets + estimating decoys.
    4. Options::iterations rounds (default 25; the paper uses 5, which do not converge on DIA pseudo-spectra): positive set = pseudo-targets that pass Options::train_fdr against the
       training decoys with SeqStep+ at c = 3/4 (decoys + 1 counted three times; without positives first without the
       +1, then with the threshold relaxed by 0.005 steps); negative set = all training decoys; class weights
       (C+, C-) from Options::c_values by 3-fold cross-validation maximizing the held-out pseudo-target discoveries
       (same SeqStep+), then a linear L2-SVM (Percolator's TRON solver) on all of the round's training rows rescores
       every stem. Only training-decoy labels enter (masking): estimating decoys are pseudo-targets throughout. The C
       grid uses all 9 pairs (the design handoff; the authors' Python restricts it to C- >= C+).
    5. Final ranking: training decoys are dropped; targets and estimating decoys are ranked by the last model and
       q-values follow from FDP = 2 (D_est + 1) / max(T, 1) (SeqStep+ at c = 2/3), monotonized from the bottom.

    The stem-level q-values are the procedure's FDR-controlled output. Everything at PSM level is heuristic (see
    Output::psm_q_values).

    Feature-symmetry audit (Options::audit, option 2 of the RESET handoff): before training, the training decoys are
    compared with an unlabelled null proxy, the pseudo-targets below the median search engine score of the training
    decoys (the training decoys are restricted to the same score range; the handoff's "lowest third of the targets"
    is mostly correct targets in runs with many identifications). Per feature
    the AUC (Mann-Whitney) and the two-sample Kolmogorov-Smirnov statistic are reported, and a cross-fitted logistic
    regression on all features gives a multivariate AUC. Features with |AUC - 0.5| > 0.05, and the set if the
    multivariate AUC exceeds 0.55, are flagged; flags are reported, nothing is dropped. The null proxy holds the
    estimating decoys, too (about a third of its nulls), which attenuates real asymmetries by about a third.
  */
  class OPENMS_DLLAPI PercolatorReset
  {
  public:
    /// Rows: one per PSM (all hits of all spectra), in any order.
    struct Input
    {
      std::vector<std::vector<double>> features; ///< feature rows (all of the same length)
      std::vector<std::string> feature_names;    ///< one name per feature column
      std::vector<bool> is_decoy;                ///< per row
      std::vector<double> search_scores;         ///< per row: search engine score, higher = better
      std::vector<Size> spectrum;                ///< per row: spectrum (competition unit) of the PSM
      std::vector<std::string> stems;            ///< per row: unmodified peptide sequence
    };

    struct Options
    {
      /// coin salt (never pick the best of several salts); unset = a fresh random salt, recorded in Output::salt
      std::optional<std::uint64_t> salt;
      int iterations = 25;              ///< training rounds (5 in the paper; they do not converge on DIA data)
      double train_fdr = 0.01;          ///< FDR for the positive set and for the cross-validation counts
      std::vector<double> c_values = {0.1, 1.0, 10.0}; ///< grid for C+ and C- (all pairs)
      std::string learner = "svm";      ///< "svm" (linear) or "mlp" (small MLP ensemble, experimental)
      /// "model" (default, user decision 24 Oct 2026): after each round a stem is represented by its PSM with the best
      /// score under the new model; "search": by its PSM with the best search engine score (the paper). With "model"
      /// masking holds, but the FDR also needs null targets and decoys to be exchangeable over their whole PSM
      /// collections, including their sizes (review r4: a null target with many more PSMs than a decoy gains from the
      /// maximum); entrapment validation decides.
      std::string representative = "model";
      bool audit = false;               ///< run the feature-symmetry audit
      std::vector<std::string> zero_features = {"enzInt"}; ///< features set to 0 (target-biased, as in the paper)
      int threads = 1;                  ///< threads for the cross-validation solves
      int mlp_members = 8;              ///< MLP ensemble size
      int mlp_hidden = 16;              ///< MLP hidden units (two hidden layers)
      int mlp_epochs = 30;              ///< MLP training epochs (the first round, or every round without warm starts)
      /// MLP warm starts: > 0 = the networks of each cross-validation task and of the final model continue from the
      /// previous round for this many epochs (a stem keeps its fold over the rounds); 0 = trained afresh every round
      int mlp_warm_epochs = 3;
      /// MLP cross-fitting: the stems (folds fixed by a hash of salt and sequence) are scored by the ensemble of the chosen
      /// C pair that did not train on their fold, in every round and at the end (no final fit); otherwise one final
      /// ensemble trained on all training stems scores everything, which ranks the training decoys it trained on too low
      bool mlp_crossfit = true;
      /// MLP cross-fitting: the fold ensembles' scores are normalized (as Percolator's folds: the fold's score at
      /// train_fdr against its training decoys -> 0, its median training decoy -> -1) only if every fold has at least this
      /// many training decoys (0 = always)
      int mlp_norm_min_decoys = 1000;
    };

    /// Role of a row's stem (decoy stems by their coin)
    enum class Role : int { TARGET = 0, ESTIMATING_DECOY = 1, TRAINING_DECOY = 2 };

    /// One considered or dropped stem
    struct Stem
    {
      std::string stem;           ///< unmodified sequence
      Size representative = 0;    ///< row of the stem's best PSM (by search engine score)
      bool decoy = false;         ///< label of the search-score representative; fixed (masking) even if
                                  ///< representative "model" later picks a PSM of the other label
      Role role = Role::TARGET;
      bool considered = false;    ///< ranked (always true: PSM-only protocol; kept for the stem table format)
      double search_score = 0.0;  ///< search engine score of the representative
      double score = 0.0;         ///< final model score of the representative
      double q_value = 1.0;       ///< RESET q-value (1 for training decoys and stems not considered)
    };

    struct AuditFeature
    {
      std::string name;
      double auc = 0.5;  ///< P(training decoy > null proxy) + P(=) / 2
      double ks = 0.0;   ///< Kolmogorov-Smirnov statistic
      bool flagged = false;
    };

    struct Output
    {
      /// per row: score of the final model (all rows, also those that lost a competition)
      std::vector<double> scores;
      /// per row, heuristic PSM-level q-value: per spectrum, the hit with the best search engine score is the PSM (the
      /// competition that selects the stems' PSMs, so a stem's representative is always reported); PSMs of
      /// training-decoy sequences (decoy rows whose sequence's coin is "training") are dropped, and the remaining PSMs (targets, estimating decoys) get
      /// q = min over lower ranks of 2 (D_est + 1) / max(T, 1). 1 for the other hits and for training decoys.
      /// Not FDR-controlled (multiplicity of PSMs per peptide).
      std::vector<double> psm_q_values;
      /// per row, heuristic PSM-level posterior error probability: 2 p / (1 - p) (at most 1), with p the
      /// estimating-decoy fraction of the ranked PSMs by isotonic (antitonic in score) regression. 1 where q is 1 by
      /// definition (other hits, training decoys).
      std::vector<double> psm_peps;
      /// per row: q-value of the row's stem (Stem::q_value)
      std::vector<double> stem_q_values;
      /// per row: role of the row's stem
      std::vector<Role> roles;
      /// all stems (considered or not)
      std::vector<Stem> stems;
      std::vector<AuditFeature> audit;
      double audit_multivariate_auc = 0.5;
      bool audit_flagged = false;
      double audit_proxy_ratio = 0.0; ///< null-proxy stems per training decoy in the audit (about 3 if all null)
      Size audit_decoys = 0, audit_proxy = 0; ///< audit group sizes (training decoys, null proxy)
      /// training summary per round: chosen C+, C-, held-out discoveries, positives
      std::vector<std::string> rounds;
      Size target_stems = 0, decoy_stems = 0, training_decoys = 0, estimating_decoys = 0;
      Size discoveries_1pct = 0;  ///< target stems with q <= 0.01
      std::vector<std::string> used_features; ///< features the learner used (non-constant, not zeroed)
      std::vector<double> weights;            ///< linear learner: raw-feature weights of the final model, bias last
      std::uint64_t salt = 0;                 ///< the effective salt and the learner (for reports)
      std::string learner;
      std::string learner_settings;           ///< learner "mlp": its settings (for reports); empty for "svm"
    };

    /// Runs RESET; throws Exception::InvalidValue on malformed input or without decoy stems.
    static Output run(const Input& input, const Options& options);

    /// SeqStep(+) q-values of a ranked list (best first): fdp_k = min(1, factor (D_k + plus_one) / max(T_k, 1)),
    /// q = minimum of fdp over ranks >= k. @p is_decoy per rank. Exposed for tests.
    static std::vector<double> seqStepQValues(const std::vector<bool>& is_decoy, double factor, double plus_one);

    /// The coin of a decoy stem: true = training decoy. Exposed for tests.
    static bool isTrainingDecoy(std::uint64_t salt, const std::string& stem);
  };
}
