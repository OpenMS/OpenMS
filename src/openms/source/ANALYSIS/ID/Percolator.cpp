// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Design note: Percolator's vendored internals (PSMDescription, ScoreHolder,
// DataSet) are PSM-shaped. This wrapper treats them as an implementation
// detail and exposes a domain-agnostic API. The PSM-specific fields (peptide
// sequence, ScanNr, ExpMass, CalcMass, proteinIds) are populated with
// synthetic placeholders for non-PSM callers — see populateVendoredFromInput_.
// If you touch the vendored code, make sure no branch in training/scoring
// parses those fields; they're opaque row IDs, nothing more.

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include "PercolatorImpl.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/METADATA/SpectrumLookup.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace OpenMS
{

namespace P = OpenMS::Internal::Percolator;

void Percolator::Impl::reset()
{
  // Second-pass PSMs are NOT owned by any DataSet — delete explicitly.
  // Their feature slots were already deallocated by scoreAndAddPSM; we only
  // need to free the PSMDescription objects.
  for (auto* psm : second_pass_psms)
  {
    if (psm) P::PSMDescription::deletePtr(psm);
  }
  second_pass_psms.clear();
  second_pass_pool.reset();

  // First-pass (training) PSM ownership lives with DataSet: DataSet::~DataSet()
  // calls PSMDescription::deletePtr on each PSM. Don't double-free.
  owned_psms.clear();
  synthesized_ids.clear();
  target_ds.reset();
  decoy_ds.reset();
  set_handler.reset();  // destructor cascades into DataSet destructors → PSM deletion
}

Percolator::Percolator()
  : DefaultParamHandler("Percolator"), impl_(std::make_unique<Impl>())
{
  defaults_.setValue("c_pos",            0.1,  "Positive regularization constant (C+).");
  defaults_.setMinFloat("c_pos",         0.0);
  defaults_.setValue("c_neg",            0.0,  "Negative regularization constant (C-). 0 = auto-select.");
  defaults_.setMinFloat("c_neg",         0.0);
  defaults_.setValue("test_fdr",         0.01, "Test FDR threshold for SVM direction selection.");
  defaults_.setMinFloat("test_fdr",      0.0);
  defaults_.setMaxFloat("test_fdr",      1.0);
  defaults_.setValue("train_fdr",        0.01, "Training FDR for SVM iterations.");
  defaults_.setMinFloat("train_fdr",     0.0);
  defaults_.setMaxFloat("train_fdr",     1.0);
  defaults_.setValue("num_iterations",   10,   "Number of SVM training iterations.");
  defaults_.setMinInt("num_iterations",  1);
  defaults_.setValue("initial_direction","",
                     "Feature name to use as the initial scoring direction. "
                     "Prefix with '-' if LOWER values are more target-like "
                     "(e.g. 'COMET:lnExpect' means higher is better; '-COMET:lnExpect' means lower is better). "
                     "Leave empty (default) to let Percolator auto-select by testing every feature "
                     "in both directions and picking the one that separates the most targets at train_fdr.");
  defaults_.setValue("pep_method",       "logistic_regression",
                     "PEP estimator: 'logistic_regression' or 'isotonic'.");
  defaults_.setValidStrings("pep_method", {"logistic_regression","isotonic"});
  defaults_.setValue("seed",             1,    "Random seed for CV splitting.");
  defaults_.setValue("target_decoy_metavalue", "target_decoy",
                     "Meta-value on PeptideHit indicating target/decoy ('target'/'decoy').");

  defaults_.setValue("normalizer", "stdv",
                     "Feature normalizer: 'stdv' (zero-mean, unit-stdv, default), "
                     "'uni' (min/max), or 'none' (no normalization).");
  defaults_.setValidStrings("normalizer", {"stdv","uni","none"});

  defaults_.setValue("train_best_positive", "true",
                     "During TRAINING, use only the best PSM per spectrum as a positive example. "
                     "Affects the SVM (weights change). Orthogonal to post_processing_tdc. "
                     "Default is 'true' because typical OpenMS inputs are multi-hit-per-spectrum "
                     "(Comet/MSGF+/Sage top-N, OSW peak groups). Inert on one-hit-per-spectrum "
                     "inputs. Note: the upstream percolator binary defaults this to off — "
                     "PercolatorAdapter preserves the CLI contract by mirroring its flag state "
                     "when writing this param.");
  defaults_.setValidStrings("train_best_positive", {"true","false"});

  defaults_.setValue("post_processing_tdc", "true",
                     "AFTER training, before PEP/q-value calculation, apply Target-Decoy "
                     "Competition: keep only the highest-scoring PSM per (scan, mass, charge) "
                     "and drop the rest from the score set. Does NOT touch the SVM — only the "
                     "score distribution used for FDR. Default is 'true' because without TDC on "
                     "multi-hit input, q-values count every candidate as an independent 'chance' "
                     "and drift optimistic. Inert on one-hit-per-spectrum inputs. Note: the "
                     "upstream percolator binary defaults this to off — PercolatorAdapter "
                     "preserves the CLI contract by mirroring its flag state when writing this "
                     "param.");
  defaults_.setValidStrings("post_processing_tdc", {"true","false"});

  defaults_.setValue("nested_xval_bins", 1,
                     "Number of nested CV bins for grid search of Cpos/Cneg. "
                     "Upstream default for cross validation is 3; 1 disables nested search.");
  defaults_.setMinInt("nested_xval_bins", 1);

  defaults_.setValue("subset_max_train", 100000,
                     "Cap on the training set size (0 = use all). Default 100000 bounds training "
                     "time on very large inputs (>1M PSMs) at negligible quality loss. Sampling "
                     "is reservoir-based, keyed by (scan, exp_mass) so all PSMs sharing a "
                     "spectrum travel together; the trained model is then applied to ALL input "
                     "rows via a second scoring pass (matches upstream Caller.cpp behavior). "
                     "Inert when the input size is at or below the cap. Note: upstream "
                     "percolator binary defaults this to 0 — PercolatorAdapter preserves the "
                     "CLI contract by passing its own int option (CLI default 0) through to "
                     "this param.");
  defaults_.setMinInt("subset_max_train", 0);

  defaults_.setValue("report_as_main_score", "none",
                     "Which Percolator output to stamp as the main PeptideHit score after rescoring. "
                     "'none' (default) leaves hit.getScore() and the PeptideIdentification score_type "
                     "untouched — callers who want the original search score preserved (or who run "
                     "their own post-processing like PercolatorAdapter) rely on this. "
                     "'q-value' sets main score to the PSM-level q-value and marks "
                     "higher_score_better=false. "
                     "'pep' sets main score to the PEP and marks higher_score_better=false. "
                     "'svm' sets main score to the SVM discriminant and marks higher_score_better=true. "
                     "All three Percolator outputs remain accessible as percolator_{score,q_value,pep} "
                     "meta values regardless of this setting.");
  defaults_.setValidStrings("report_as_main_score", {"none","q-value","pep","svm"});

  defaults_.setValue("use_pi0", "true",
                     "Enable pi0 correction when computing q-values and PEPs (default). "
                     "When true, pi0 is estimated per CV fold via bootstrap spline fit on "
                     "the target score distribution using decoys as the null. When false, "
                     "pi0 is fixed at 1.0 (pure target-decoy FDR, no bootstrap estimation). "
                     "Affects both Scores::calcQvals and Scores::calcPep. Upstream percolator "
                     "has no CLI flag for this; the in-library Scores(bool usePi0) "
                     "constructor is the only toggle, matched here.");
  defaults_.setValidStrings("use_pi0", {"true","false"});

  defaultsToParam_();
}

Percolator::~Percolator()
{
  if (impl_) impl_->reset();
}

void Percolator::updateMembers_()
{
  impl_->c_pos          = param_.getValue("c_pos");
  impl_->c_neg          = param_.getValue("c_neg");
  impl_->test_fdr       = param_.getValue("test_fdr");
  impl_->train_fdr      = param_.getValue("train_fdr");
  impl_->num_iterations = static_cast<int>(param_.getValue("num_iterations"));
  impl_->seed           = static_cast<int>(param_.getValue("seed"));
  impl_->pep_method     = param_.getValue("pep_method").toString();
  impl_->normalizer     = param_.getValue("normalizer").toString();
  impl_->train_best_positive = (param_.getValue("train_best_positive").toString() == "true");
  impl_->post_processing_tdc = (param_.getValue("post_processing_tdc").toString() == "true");
  impl_->nested_xval_bins    = static_cast<int>(param_.getValue("nested_xval_bins"));
  impl_->subset_max_train    = static_cast<int>(param_.getValue("subset_max_train"));
  impl_->initial_direction   = param_.getValue("initial_direction").toString();
  impl_->use_pi0             = (param_.getValue("use_pi0").toString() == "true");
  impl_->report_as_main_score = param_.getValue("report_as_main_score").toString();
}

const std::vector<std::vector<double>>& Percolator::getSvmWeights() const
{
  return impl_->svm_weights;
}

struct RowKeys
{
  int    scan;
  int    spec_file;
  double exp_mass;
  double calc_mass;
};

// Populate a single PSMDescription from one row of RescoreInput. The row's
// feature pointer comes from SetHandler's FeatureMemoryPool.
static P::PSMDescription* makePSMFromRow(
  const std::vector<double>& feats,
  const RowKeys& keys,
  size_t row_index,
  std::string& synth_id,
  P::FeatureMemoryPool& pool)
{
  auto* psm = new P::PSMDescription();
  // Synthesize an opaque row ID as the "peptide" string — never parsed.
  char buf[32];
  std::snprintf(buf, sizeof(buf), "row_%08zu", row_index);
  synth_id = buf;
  psm->setId(synth_id);
  psm->setPeptide(synth_id);
  psm->scan        = static_cast<unsigned int>(keys.scan);
  psm->specFileNr  = static_cast<unsigned int>(keys.spec_file);
  psm->expMass     = keys.exp_mass;
  psm->calcMass    = keys.calc_mass;
  psm->retentionTime_ = 0.0;

  // Allocate feature row from the pool.
  double* feat_ptr = pool.allocate();
  for (size_t j = 0; j < feats.size(); ++j)
  {
    feat_ptr[j] = feats[j];
  }
  psm->features = feat_ptr;
  return psm;
}

// Build the RowKeys tuple for row `i`, applying the same precedence and
// defaults used during PSM registration. Centralized so the reservoir-sampling
// pass and the PSM-building pass use identical (scan, exp_mass) keys.
static RowKeys makeRowKeys(const RescoreInput& input, size_t i)
{
  RowKeys keys;
  if (!input.scan_numbers.empty())       keys.scan = input.scan_numbers[i];
  else if (!input.cv_group_keys.empty()) keys.scan = input.cv_group_keys[i];
  else                                   keys.scan = static_cast<int>(i);
  keys.spec_file = input.spec_file_numbers.empty() ? 0 : input.spec_file_numbers[i];
  keys.exp_mass  = input.exp_masses.empty()        ? 0.0 : input.exp_masses[i];
  keys.calc_mass = input.calc_masses.empty()       ? 0.0 : input.calc_masses[i];
  return keys;
}

// Materialize the (scan, exp_mass) arrays and hand off to the vendored
// reservoir-sampling helper (thirdparty/percolator/InProcessSampler.cpp).
// See that file's header for provenance and re-sync notes.
static std::vector<size_t> selectTrainingSubset(
  const RescoreInput& input, size_t n_rows, int subset_max_train)
{
  const size_t max_train = (subset_max_train <= 0)
    ? 0u : static_cast<size_t>(subset_max_train);

  std::vector<int> scans(n_rows);
  std::vector<double> exp_masses(n_rows, 0.0);
  for (size_t i = 0; i < n_rows; ++i)
  {
    const RowKeys k = makeRowKeys(input, i);
    scans[i]      = k.scan;
    exp_masses[i] = k.exp_mass;
  }
  return P::sampleTrainingRowIndices(scans, exp_masses, max_train);
}

RescoreOutput Percolator::rescore(const RescoreInput& input)
{
  if (input.features.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "RescoreInput has empty features matrix", "");
  }
  if (input.features.size() != input.is_decoy.size())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "features and is_decoy length mismatch", "");
  }
  const size_t n_rows = input.features.size();
  const size_t n_feat = input.features.front().size();
  for (const auto& row : input.features)
  {
    if (row.size() != n_feat)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "feature rows must all have the same length", "");
    }
  }

  impl_->reset();
  P::PseudoRandom::setSeed(static_cast<unsigned long>(impl_->seed));

  // Validate optional PIN-compat fields up front (before sampling, which
  // consults input.scan_numbers / exp_masses).
  auto validate_opt = [&](const auto& v, const char* name)
  {
    if (!v.empty() && v.size() != n_rows)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        std::string(name) + " must be empty or match n_rows", "");
    }
  };
  validate_opt(input.cv_group_keys,      "cv_group_keys");
  validate_opt(input.scan_numbers,       "scan_numbers");
  validate_opt(input.spec_file_numbers,  "spec_file_numbers");
  validate_opt(input.exp_masses,         "exp_masses");
  validate_opt(input.calc_masses,        "calc_masses");

  // Reservoir-sample the training subset (upstream Caller.cpp pattern). When
  // subset_max_train is 0 or >= n_rows, training_rows is the identity [0..n).
  const std::vector<size_t> training_rows =
    selectTrainingSubset(input, n_rows, impl_->subset_max_train);
  const bool did_subsample = training_rows.size() < n_rows;

  // Set up SetHandler + feature pool. Pass 0 for maxPSMs because sampling is
  // already done above (the SetHandler cap only matters for its own readPSMs
  // path, which we don't use).
  impl_->set_handler = std::make_unique<P::SetHandler>(0);
  P::FeatureMemoryPool& pool = impl_->set_handler->getFeaturePool();
  pool.createPool(static_cast<size_t>(n_feat));

  // Register feature names (used only for logging).
  P::DataSet::resetFeatureNames();
  P::FeatureNames& fn = P::DataSet::getFeatureNames();
  for (size_t j = 0; j < n_feat; ++j)
  {
    const std::string nm = (j < input.feature_names.size())
      ? static_cast<std::string>(input.feature_names[j])
      : "feat_" + std::to_string(j);
    fn.insertFeature(nm);
  }
  fn.initFeatures();

  // Create target + decoy DataSets.
  impl_->target_ds = std::make_unique<P::DataSet>();
  impl_->decoy_ds  = std::make_unique<P::DataSet>();
  impl_->target_ds->setLabel(P::LabelType::TARGET);
  impl_->decoy_ds->setLabel(P::LabelType::DECOY);
  // DataSet::initFeatureTables is declared in DataSet.h but upstream never
  // defines it. Not needed — feature storage comes from FeatureMemoryPool.

  // Populate rows (training subset only — the second pass below scores the
  // full set using the trained model, mirroring upstream Caller.cpp flow).
  impl_->synthesized_ids.resize(n_rows);
  impl_->owned_psms.reserve(training_rows.size());
  for (size_t i : training_rows)
  {
    const RowKeys keys = makeRowKeys(input, i);
    auto* psm = makePSMFromRow(input.features[i], keys, i,
                               impl_->synthesized_ids[i], pool);
    impl_->owned_psms.push_back(psm);
    if (input.is_decoy[i])
    {
      impl_->decoy_ds->registerPsm(psm);
    }
    else
    {
      impl_->target_ds->registerPsm(psm);
    }
  }

  impl_->set_handler->push_back_dataset(impl_->target_ds.release());
  impl_->set_handler->push_back_dataset(impl_->decoy_ds.release());

  // Normalizer.
  int norm_type = P::Normalizer::STDV;
  if      (impl_->normalizer == "uni")  norm_type = P::Normalizer::UNI;
  else if (impl_->normalizer == "none") norm_type = P::Normalizer::NONORM;
  P::Normalizer::setType(norm_type);
  P::Normalizer* normalizer = P::Normalizer::getNormalizer();

  // SanityCheck: "default" fingerprint → generic PSM-level sanity.
  // Reset any prior static state (initial direction) so successive rescore()
  // calls on the same instance don't leak configuration across calls.
  P::SanityCheck::setInitDefaultDir(0);
  P::SanityCheck::setInitDefaultDirName(impl_->initial_direction);

  P::SanityCheck* sanity = P::SanityCheck::initialize("");
  if (!sanity)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Could not initialize Percolator SanityCheck", "");
  }
  // Resolve initial_direction name → feature index. In the subprocess path
  // this is called from SetHandler::readPSMs (which we bypass). Must run
  // after feature names are registered (above) and before preIterationSetup.
  try
  {
    sanity->checkAndSetDefaultDir();
  }
  catch (const std::exception& e)
  {
    delete sanity;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      std::string("initial_direction not resolved: ") + e.what(),
      impl_->initial_direction);
  }

  impl_->set_handler->normalizeFeatures(normalizer);

  // Build the full score set and populate from SetHandler.
  P::Scores all_scores(/*usePi0=*/impl_->use_pi0);
  all_scores.populateWithPSMs(*impl_->set_handler);

  // Cross-validation setup and training.
  P::CrossValidation cv(
    /*quickValidation=*/false,
    /*reportPerformanceEachIteration=*/false,
    /*testFdr=*/impl_->test_fdr,
    /*selectionFdr=*/impl_->train_fdr,
    /*initialSelectionFdr=*/impl_->train_fdr,
    /*selectedCpos=*/impl_->c_pos,
    /*selectedCneg=*/impl_->c_neg,
    /*niter=*/static_cast<unsigned int>(impl_->num_iterations),
    /*usePi0=*/impl_->use_pi0,
    /*nestedXvalBins=*/static_cast<unsigned int>(impl_->nested_xval_bins),
    /*trainBestPositive=*/impl_->train_best_positive,
    /*numThreads=*/1,
    /*skipNormalizeScores=*/false,
    /*decoyFractionTraining=*/1.0,
    /*numFolds=*/3);

  try
  {
    int first_positives = cv.preIterationSetup(
      all_scores, sanity, normalizer, impl_->set_handler->getFeaturePool());
    (void)first_positives;
    cv.train(normalizer);
    cv.postIterationProcessing(all_scores, sanity);

    // Capture average SVM weights per fold BEFORE the TDC weeding step
    // (weeding mutates all_scores and may invalidate per-fold weight meaning).
    // getAvgWeights returns un-normalized weights in raw feature space.
    std::vector<double> avg_weights;
    cv.getAvgWeights(avg_weights, normalizer);
    impl_->svm_weights.clear();
    if (!avg_weights.empty()) impl_->svm_weights.push_back(avg_weights);

    // Second pass: when we subsampled for training, score the FULL input set
    // using the trained weights. Mirrors upstream Caller.cpp (lines 1196-1236
    // in rel-3-08-01): reset allScores, rebuild via readAndScoreTab, then
    // postMergeStep + calcQvals + normalizeScores. This is what makes the
    // subset_max_train docstring claim "trained model applied to ALL PSMs" true.
    if (did_subsample)
    {
      // Fresh feature pool — the existing `pool` is in use by training PSMs.
      impl_->second_pass_pool = std::make_unique<P::FeatureMemoryPool>();
      impl_->second_pass_pool->createPool(static_cast<size_t>(n_feat));
      P::FeatureMemoryPool& full_pool = *impl_->second_pass_pool;

      // Rebuild all_scores from scratch on every row, scored via avg_weights.
      // The training-set scores/q-values/PEPs produced by postIterationProcessing
      // above are discarded — the full-set values supersede them.
      all_scores = P::Scores(impl_->use_pi0);
      impl_->second_pass_psms.reserve(n_rows);

      for (size_t i = 0; i < n_rows; ++i)
      {
        const RowKeys keys = makeRowKeys(input, i);
        auto* psm = makePSMFromRow(input.features[i], keys, i,
                                   impl_->synthesized_ids[i], full_pool);
        impl_->second_pass_psms.push_back(psm);

        P::ScoreHolder sh;
        sh.pPSM  = psm;
        sh.label = input.is_decoy[i] ? P::LabelType::DECOY : P::LabelType::TARGET;
        // scoreAndAddPSM computes sh.score = Σ feat[j]·w[j] + w[bias], pushes
        // the ScoreHolder into all_scores, and deallocates feat back to the
        // pool (one-way trip — features not needed again).
        all_scores.scoreAndAddPSM(sh, avg_weights, full_pool);
      }

      // Sort + count + pi0 estimation on the full set.
      all_scores.postMergeStep();
      // Upstream sequence: calcQvals on full set, then normalizeScores rescales
      // SVM scores so median-decoy maps to 0 and first-FDR-crossing maps to 1.
      all_scores.calcQvals(impl_->test_fdr);
      all_scores.normalizeScores(impl_->test_fdr, avg_weights);
    }

    // Optional Target-Decoy Competition: best PSM per spectrum by score.
    // Must happen after postIterationProcessing / second pass (so final scores
    // are set) but before calcPep (PEPs depend on the deduplicated score set).
    if (impl_->post_processing_tdc)
    {
      all_scores.weedOutRedundantTDC();
    }

    // Post-train PEP estimation: matches upstream Caller.cpp which calls
    // calcPep after postIterationProcessing. Without this, every PEP stays
    // at its default (uninitialized) value. Args match upstream defaults:
    // non-IRLS (logistic regression), non-interpolating, non-PAVA.
    if (impl_->pep_method == "isotonic")
    {
      all_scores.calcPep(/*spline=*/false, /*interp=*/false, /*pava=*/true);
    }
    else
    {
      all_scores.calcPep(/*spline=*/false, /*interp=*/false, /*pava=*/false);
    }
  }
  catch (const std::exception& e)
  {
    delete sanity;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      std::string("Percolator training failed: ") + e.what(), "");
  }
  delete sanity;

  // Extract outputs. Each ScoreHolder's pPSM->peptide is the synthetic row ID.
  RescoreOutput out;
  out.scores.assign(n_rows, 0.0);
  out.q_values.assign(n_rows, 1.0);
  out.peps.assign(n_rows, 1.0);

  size_t matched = 0;
  double score_min = std::numeric_limits<double>::infinity();
  double score_max = -std::numeric_limits<double>::infinity();
  for (const auto& sh : all_scores)
  {
    if (sh.pPSM == nullptr) continue;
    const std::string& id = sh.pPSM->getFullPeptide();
    // Parse "row_NNNNNNNN" → row index
    if (id.size() < 12 || id.compare(0, 4, "row_") != 0) continue;
    char* end = nullptr;
    const unsigned long row = std::strtoul(id.c_str() + 4, &end, 10);
    if (row >= n_rows) continue;
    out.scores[row]   = sh.score;
    out.q_values[row] = sh.q;
    out.peps[row]     = sh.pep;
    ++matched;
    if (sh.score < score_min) score_min = sh.score;
    if (sh.score > score_max) score_max = sh.score;
  }
  OPENMS_LOG_DEBUG << "[Percolator] extracted " << matched << " / " << n_rows
                   << " rows from all_scores (size=" << all_scores.size() << "); "
                   << "score range [" << score_min << ", " << score_max << "]"
                   << std::endl;

  return out;
}

namespace
{
  /// Numeric-feature auto-discovery from a PeptideHit's meta values.
  /// Excludes: the target/decoy meta, OpenMS internal scoring keys, and
  /// keys previously stamped by a Percolator run.
  StringList discoverFeatureNames_(const PeptideHit& hit, const String& td_meta)
  {
    static const std::vector<String> blocklist {
      "percolator_score", "percolator_q_value", "percolator_pep",
      "q-value", "q-value_score", "score", "RT", "MZ",
      "target_decoy", "isDecoy"
    };
    StringList out;
    std::vector<String> keys;
    hit.getKeys(keys);
    for (const String& k : keys)
    {
      if (k == td_meta) continue;
      bool skipped = false;
      for (const String& b : blocklist)
      {
        if (k == b) { skipped = true; break; }
      }
      if (skipped) continue;
      const DataValue& v = hit.getMetaValue(k);
      if (v.valueType() == DataValue::DOUBLE_VALUE ||
          v.valueType() == DataValue::INT_VALUE)
      {
        out.push_back(k);
      }
    }
    std::sort(out.begin(), out.end());
    return out;
  }
}

void Percolator::fillPINCompatibleFields(
    const std::vector<PeptideIdentification>& peptide_ids,
    bool flatten_hits,
    RescoreInput& input)
{
  if (peptide_ids.empty()) return;

  // Count rows to pre-size the vectors — matches the high-level rescore's
  // iteration pattern (one row per PeptideHit if flatten_hits, else one per
  // pid using its first hit).
  size_t n_rows = 0;
  for (const auto& pid : peptide_ids)
  {
    if (pid.getHits().empty()) continue;
    n_rows += flatten_hits ? pid.getHits().size() : 1;
  }
  input.scan_numbers.assign(n_rows, 0);
  input.spec_file_numbers.assign(n_rows, 0);
  input.exp_masses.assign(n_rows, 0.0);
  input.calc_masses.assign(n_rows, 0.0);

  // Build a lookup so spec_file numbers stay stable across invocations.
  std::unordered_map<std::string, int> spec_file_to_idx;

  // Spec-lookup regex is derived from the first pid's scan identifier (same
  // as PercolatorInfile::preparePin_).
  const String first_sid = PercolatorInfile::getScanIdentifier(peptide_ids.front(), 0);
  const boost::regex scan_regex(
    SpectrumLookup::getRegExFromNativeID(first_sid));

  size_t row = 0;
  size_t pid_index = 0;
  for (const auto& pid : peptide_ids)
  {
    ++pid_index;
    if (pid.getHits().empty()) continue;

    const String scan_identifier = PercolatorInfile::getScanIdentifier(pid, pid_index);
    const Int scan_number = SpectrumLookup::extractScanNumber(
        scan_identifier, scan_regex, /*no_error=*/true);

    const std::string file_key =
      static_cast<std::string>(pid.getMetaValue("file_origin", String())) +
      "|" +
      static_cast<std::string>(pid.getMetaValue("id_merge_index", String()));

    int spec_file = 0;
    auto it = spec_file_to_idx.find(file_key);
    if (it != spec_file_to_idx.end())
    {
      spec_file = it->second;
    }
    else
    {
      spec_file = static_cast<int>(spec_file_to_idx.size());
      spec_file_to_idx[file_key] = spec_file;
    }

    const double exp_mass = pid.hasMZ() ? pid.getMZ() : 0.0;

    const size_t hit_count = flatten_hits ? pid.getHits().size() : 1u;
    for (size_t h = 0; h < hit_count; ++h)
    {
      const PeptideHit& hit = pid.getHits()[h];
      double calc_mass = 0.0;
      if (hit.metaValueExists("CalcMass"))
      {
        calc_mass = hit.getMetaValue("CalcMass");
      }
      else
      {
        try { calc_mass = hit.getSequence().getMonoWeight(); }
        catch (...) { calc_mass = 0.0; }
      }

      input.scan_numbers[row]      = scan_number;
      input.spec_file_numbers[row] = spec_file;
      input.exp_masses[row]        = exp_mass;
      input.calc_masses[row]       = calc_mass;
      ++row;
    }
  }
}

void Percolator::rescore(std::vector<PeptideIdentification>& peptide_ids,
                         const StringList& feature_names)
{
  if (peptide_ids.empty()) return;

  const String td_meta = param_.getValue("target_decoy_metavalue").toString();

  // Auto-discover features from first hit if not provided.
  StringList features = feature_names;
  if (features.empty())
  {
    const PeptideIdentification& first = peptide_ids.front();
    if (first.getHits().empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "first PeptideIdentification has no hits", "");
    }
    features = discoverFeatureNames_(first.getHits().front(), td_meta);
    if (features.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "no numeric features found on first PeptideHit", "");
    }
  }

  // Build RescoreInput. Each hit contributes one row.
  RescoreInput ri;
  ri.feature_names = features;
  struct HitLoc { size_t pid_i; size_t hit_i; };
  std::vector<HitLoc> hit_locs;
  ri.features.reserve(peptide_ids.size());
  ri.is_decoy.reserve(peptide_ids.size());
  ri.cv_group_keys.reserve(peptide_ids.size());

  std::hash<std::string> str_hash;
  for (size_t i = 0; i < peptide_ids.size(); ++i)
  {
    auto& hits = peptide_ids[i].getHits();
    for (size_t j = 0; j < hits.size(); ++j)
    {
      const PeptideHit& hit = hits[j];
      if (!hit.metaValueExists(td_meta))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "PeptideHit missing target/decoy meta value", td_meta);
      }

      std::vector<double> row;
      row.reserve(features.size());
      for (const String& f : features)
      {
        if (!hit.metaValueExists(f))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "PeptideHit missing feature meta value", f);
        }
        row.push_back(static_cast<double>(hit.getMetaValue(f)));
      }
      ri.features.push_back(std::move(row));
      ri.is_decoy.push_back(hit.getMetaValue(td_meta).toString() == "decoy");

      // CV group key: hash (identifier, rt).
      std::string spec_key = peptide_ids[i].getIdentifier() + "|" +
                             std::to_string(peptide_ids[i].getRT());
      ri.cv_group_keys.push_back(static_cast<int>(str_hash(spec_key) & 0x7fffffff));

      hit_locs.push_back({i, j});
    }
  }

  RescoreOutput ro = rescore(ri);  // low-level

  // Stamp back onto the hits.
  for (size_t row = 0; row < ro.scores.size(); ++row)
  {
    const auto& loc = hit_locs[row];
    PeptideHit& hit = peptide_ids[loc.pid_i].getHits()[loc.hit_i];
    hit.setMetaValue("percolator_score",   ro.scores[row]);
    hit.setMetaValue("percolator_q_value", ro.q_values[row]);
    hit.setMetaValue("percolator_pep",     ro.peps[row]);
  }

  // Optional: promote one of the three Percolator outputs to the main hit score
  // and update PeptideIdentification metadata accordingly.
  const std::string& as_main = impl_->report_as_main_score;
  if (as_main != "none")
  {
    // Track which pids we've updated so score_type/higher_score_better gets set
    // exactly once even when a pid has multiple hits.
    std::unordered_set<size_t> pids_touched;
    pids_touched.reserve(peptide_ids.size());

    for (size_t row = 0; row < ro.scores.size(); ++row)
    {
      const auto& loc = hit_locs[row];
      PeptideIdentification& pid = peptide_ids[loc.pid_i];
      PeptideHit& hit = pid.getHits()[loc.hit_i];

      if (as_main == "q-value")   hit.setScore(ro.q_values[row]);
      else if (as_main == "pep")  hit.setScore(ro.peps[row]);
      else /* "svm" */            hit.setScore(ro.scores[row]);

      if (pids_touched.insert(loc.pid_i).second)
      {
        if (as_main == "svm")
        {
          pid.setScoreType("svm");
          pid.setHigherScoreBetter(true);
        }
        else if (as_main == "pep")
        {
          pid.setScoreType("Posterior Error Probability");
          pid.setHigherScoreBetter(false);
        }
        else  // "q-value"
        {
          pid.setScoreType("q-value");
          pid.setHigherScoreBetter(false);
        }
      }
    }

    // Re-rank hits within each touched pid under the new score orientation.
    for (size_t i : pids_touched)
    {
      peptide_ids[i].sort();
    }
  }
}

} // namespace OpenMS
