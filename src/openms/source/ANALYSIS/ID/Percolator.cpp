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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

namespace P = OpenMS::Internal::Percolator;

namespace
{
  // Percolator's sanity check, recording whether it replaced the learned
  // weights by the default direction (the best single feature): validateDirection()
  // does that whenever it returns false, i.e. no target at q < test_fdr, fewer
  // targets than the default direction, or a wrong sign of the initial direction.
  class RecordingSanityCheck : public P::SanityCheck
  {
  public:
    bool validateDirection(std::vector<std::vector<double>>& w) override
    {
      const bool valid = P::SanityCheck::validateDirection(w);
      used_default_direction = used_default_direction || !valid;
      return valid;
    }

    bool used_default_direction = false;
  };
}

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
  // set_handler's destructor cascades: ~SetHandler → ~DataSet → PSM deletion,
  // and ~SetHandler ALSO invokes DataSet::resetFeatureNames() which flips the
  // process-global FeatureNames::numFeatures back to 0. The next train()/
  // score()/rescore() call relies on that reset (see the resetFeatureNames
  // calls at the top of each entry point). If a future change replaces
  // SetHandler with something that doesn't cascade through, FeatureNames
  // must be reset manually here.
  set_handler.reset();

  // svm_weights is intentionally NOT cleared here — getSvmWeights() is
  // documented to retain the last training's weights across subsequent
  // score() calls.
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
                     "Posterior error probability estimator: "
                     "'logistic_regression' (target-decoy competition + isotonic fit; "
                     "legacy name, not actually logistic regression), "
                     "'isotonic' (target-q → pep monotone regression), or "
                     "'nonparametric' (kernel smoothing).");
  defaults_.setValidStrings("pep_method", {"logistic_regression","isotonic","nonparametric"});
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

  defaults_.setValue("num_threads",       3,
                     "Maximum number of OpenMP threads for cross-validation; never more than the "
                     "threads the tool may use (-threads). With the default settings, cross-validation "
                     "trains 9 models (3 folds x 3 Cneg values) in parallel. Default and range match upstream.");
  defaults_.setMinInt("num_threads",      1);
  defaults_.setMaxInt("num_threads",      128);

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
  P::Normalizer::resetNormalizer();
  P::Globals::clean();
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
  impl_->num_threads         = static_cast<int>(param_.getValue("num_threads"));
  impl_->report_as_main_score = param_.getValue("report_as_main_score").toString();
}

const std::vector<std::vector<double>>& Percolator::getSvmWeights() const
{
  return impl_->svm_weights;
}

double Percolator::getPi0() const
{
  return impl_->last_pi0;
}

bool Percolator::usedDefaultDirection() const
{
  return impl_->last_used_default_direction;
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
// The threads of a parallel region whose sorts and pi0 bootstraps run as its tasks: at most 12,
// as the cross-validation's (thirdparty/percolator/CrossValidation.cpp: numTeamThreads), more
// would mostly wait for tasks; never more than OpenMP allows here.
[[maybe_unused]] static int teamThreads(int num_threads)
{
#ifdef _OPENMP
  return std::max(1, std::min({num_threads, 12, omp_get_max_threads()}));
#else
  (void)num_threads;
  return 1;
#endif
}

// Rethrows the first error, in iteration order, caught in the iterations of a parallel loop.
static void rethrowFirstError(const std::vector<std::exception_ptr>& errors)
{
  for (const std::exception_ptr& error : errors)
  {
    if (error) std::rethrow_exception(error);
  }
}

// Fills a PSMDescription from one row of RescoreInput (n_feat values at feats); its feature row
// (psm.features) is allocated already.
static void fillPSMFromRow(
  P::PSMDescription& psm,
  const double* feats,
  size_t n_feat,
  const RowKeys& keys,
  size_t row_index,
  std::string& synth_id)
{
  // Synthesize an opaque row ID as the "peptide" string: "row_%08zu" (written without snprintf, which cost more than
  // the rest of the row); the scores are mapped back to the rows by it.
  char digits[24];
  size_t n_digits = 0;
  for (size_t v = row_index; v != 0 || n_digits < 8; v /= 10) digits[n_digits++] = static_cast<char>('0' + v % 10);
  synth_id.assign("row_");
  while (n_digits > 0) synth_id += digits[--n_digits];
  psm.setId(synth_id);
  psm.setPeptide(synth_id);
  psm.scan        = static_cast<unsigned int>(keys.scan);
  psm.specFileNr  = static_cast<unsigned int>(keys.spec_file);
  psm.expMass     = keys.exp_mass;
  psm.calcMass    = keys.calc_mass;
  psm.retentionTime_ = 0.0;

  double* feat_ptr = psm.features;
  for (size_t j = 0; j < n_feat; ++j)
  {
    feat_ptr[j] = feats[j];
  }
}

static P::PSMDescription* makePSMFromRow(
  const double* feats,
  size_t n_feat,
  const RowKeys& keys,
  size_t row_index,
  std::string& synth_id,
  P::FeatureMemoryPool& pool)
{
  auto* psm = new P::PSMDescription();
  // Allocate feature row from the pool.
  psm->features = pool.allocate();
  fillPSMFromRow(*psm, feats, n_feat, keys, row_index, synth_id);
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

// Build a new RescoreInput by selecting the given row indices. Preserves
// the optional PIN-compat / cv_group_keys fields when they're populated.
static RescoreInput makeSubsetInput_(
  const RescoreInput& input, const std::vector<size_t>& indices)
{
  RescoreInput out;
  out.feature_names = input.feature_names;
  const bool dense = input.features.empty();
  const size_t n_feat = input.featureCount();
  if (dense) out.feature_matrix.reserve(indices.size() * n_feat);
  else       out.features.reserve(indices.size());
  out.is_decoy.reserve(indices.size());
  const bool has_cv       = !input.cv_group_keys.empty();
  const bool has_scan     = !input.scan_numbers.empty();
  const bool has_specfile = !input.spec_file_numbers.empty();
  const bool has_exp      = !input.exp_masses.empty();
  const bool has_calc     = !input.calc_masses.empty();
  if (has_cv)       out.cv_group_keys.reserve(indices.size());
  if (has_scan)     out.scan_numbers.reserve(indices.size());
  if (has_specfile) out.spec_file_numbers.reserve(indices.size());
  if (has_exp)      out.exp_masses.reserve(indices.size());
  if (has_calc)     out.calc_masses.reserve(indices.size());
  for (size_t i : indices)
  {
    if (dense) out.feature_matrix.insert(out.feature_matrix.end(), input.row(i), input.row(i) + n_feat);
    else       out.features.push_back(input.features[i]);
    out.is_decoy.push_back(input.is_decoy[i]);
    if (has_cv)       out.cv_group_keys.push_back(input.cv_group_keys[i]);
    if (has_scan)     out.scan_numbers.push_back(input.scan_numbers[i]);
    if (has_specfile) out.spec_file_numbers.push_back(input.spec_file_numbers[i]);
    if (has_exp)      out.exp_masses.push_back(input.exp_masses[i]);
    if (has_calc)     out.calc_masses.push_back(input.calc_masses[i]);
  }
  return out;
}

// Dimension + optional-field checks shared by rescore/train/score entry points.
static void validateRescoreInput_(const RescoreInput& input)
{
  if (input.features.empty() && input.feature_matrix.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "RescoreInput has empty features matrix", "");
  }
  if (!input.features.empty() && !input.feature_matrix.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "RescoreInput has rows in both features and feature_matrix", "");
  }
  if (input.features.empty() &&
      (input.feature_names.empty() || input.feature_matrix.size() / input.feature_names.size() != input.is_decoy.size() ||
       input.feature_matrix.size() % input.feature_names.size() != 0))
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "feature_matrix must hold is_decoy.size() rows of feature_names.size() values", "");
  }
  if (input.rowCount() != input.is_decoy.size())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "features and is_decoy length mismatch", "");
  }
  const size_t n_rows = input.rowCount();
  const size_t n_feat = input.featureCount();
  for (const auto& row : input.features)
  {
    if (row.size() != n_feat)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "feature rows must all have the same length", "");
    }
  }
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
}

// Feature-count + feature-names compatibility between a scoring input and a
// trained model. BOTH sides must supply feature_names — skipping the check on
// one-sided empty is a silent misalignment risk (the original guard was too
// permissive; this enforces the precondition).
static void validateModelCompat_(
  const RescoreInput& input, const PercolatorModel& model)
{
  const size_t n_feat = input.featureCount();
  if (model.weights.size() != n_feat + 1)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "PercolatorModel weight count does not match input feature count "
      "(expected n_features + 1 for bias)",
      std::to_string(model.weights.size()) + " vs " +
      std::to_string(n_feat + 1));
  }
  if (model.feature_names.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "PercolatorModel has no feature_names — cannot validate input alignment",
      "");
  }
  if (input.feature_names.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "RescoreInput.feature_names required when scoring against a trained model",
      "");
  }
  if (model.feature_names.size() != input.feature_names.size())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "PercolatorModel feature_names size mismatch with input",
      std::to_string(model.feature_names.size()) + " vs " +
      std::to_string(input.feature_names.size()));
  }
  for (size_t i = 0; i < model.feature_names.size(); ++i)
  {
    if (model.feature_names[i] != input.feature_names[i])
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "PercolatorModel feature_names differ from input at column " +
        std::to_string(i),
        static_cast<std::string>(model.feature_names[i]) + " != " +
        static_cast<std::string>(input.feature_names[i]));
    }
  }
}

PercolatorModel Percolator::train(const RescoreInput& input)
{
  // Mark pi0 invalid up front: train() doesn't run scoring / q-value
  // estimation, so it never produces a pi0. Resetting here matches the
  // pattern in rescore()/score() and ensures getPi0() == -1.0 reliably
  // means "no usable pi0 from the current state" (rather than a stale
  // value left over from a prior successful rescore()/score() call).
  // Also covers throws from validateRescoreInput_ below.
  impl_->last_pi0 = -1.0;
  impl_->last_used_default_direction = false;

  validateRescoreInput_(input);

  impl_->reset();
  P::PseudoRandom::setSeed(static_cast<unsigned long>(impl_->seed));

  const size_t n_rows = input.rowCount();
  const size_t n_feat = input.featureCount();

  // SetHandler + feature pool. maxPSMs=0 because we never use SetHandler's own
  // readPSMs path (it's a file-level loader); we register rows directly.
  impl_->set_handler = std::make_unique<P::SetHandler>(0);
  P::FeatureMemoryPool& pool = impl_->set_handler->getFeaturePool();
  pool.createPool(static_cast<size_t>(n_feat));

  // Register feature names. This ALSO sets the static FeatureNames::numFeatures
  // that Scores::scoreAndAddPSM consults — i.e., it is NOT logging-only; wiring
  // it here keeps training and (later) scoring on the same feature count.
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

  impl_->target_ds = std::make_unique<P::DataSet>();
  impl_->decoy_ds  = std::make_unique<P::DataSet>();
  impl_->target_ds->setLabel(P::LabelType::TARGET);
  impl_->decoy_ds->setLabel(P::LabelType::DECOY);

  // Populate ALL rows of the training input (caller handles any subsetting).
  impl_->synthesized_ids.resize(n_rows);
  impl_->owned_psms.reserve(n_rows);
  for (size_t i = 0; i < n_rows; ++i)
  {
    const RowKeys keys = makeRowKeys(input, i);
    auto* psm = makePSMFromRow(input.row(i), n_feat, keys, i,
                               impl_->synthesized_ids[i], pool);
    impl_->owned_psms.push_back(psm);
    if (input.is_decoy[i]) impl_->decoy_ds->registerPsm(psm);
    else                   impl_->target_ds->registerPsm(psm);
  }

  impl_->set_handler->push_back_dataset(impl_->target_ds.release());
  impl_->set_handler->push_back_dataset(impl_->decoy_ds.release());

  int norm_type = P::Normalizer::STDV;
  if      (impl_->normalizer == "uni")  norm_type = P::Normalizer::UNI;
  else if (impl_->normalizer == "none") norm_type = P::Normalizer::NONORM;
  P::Globals::getInstance()->setNoTerminate(true);
  P::Normalizer::resetNormalizer();
  P::Normalizer::setType(norm_type);
  // normalizeFeatures() owns Normalizer construction via its out-parameter.
  P::Normalizer* normalizer = nullptr;

  // SanityCheck: reset static state so successive calls don't leak configuration.
  P::SanityCheck::setInitDefaultDir(0);
  P::SanityCheck::setInitDefaultDirName(impl_->initial_direction);
  // The plain SanityCheck that SanityCheck::initialize() returns, recording
  // whether it replaced the learned weights (see usedDefaultDirection()).
  // unique_ptr releases it on every exit path (including throws).
  auto sanity = std::make_unique<RecordingSanityCheck>();
  try
  {
    sanity->checkAndSetDefaultDir();
  }
  catch (const std::exception& e)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      std::string("initial_direction not resolved: ") + e.what(),
      impl_->initial_direction);
  }

  impl_->set_handler->normalizeFeatures(normalizer);

  P::Scores train_scores(/*usePi0=*/impl_->use_pi0);
  train_scores.populateWithPSMs(*impl_->set_handler);

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
    /*numThreads=*/static_cast<unsigned int>(impl_->num_threads),
    /*skipNormalizeScores=*/false,
    /*decoyFractionTraining=*/1.0,
    /*numFolds=*/3);

  std::vector<double> avg_weights;
  try
  {
    // CV uses sanity transiently per call (no storage); the unique_ptr in
    // this scope outlives cv, so the bare pointer stays valid.
    cv.preIterationSetup(train_scores, sanity.get(), normalizer,
                         impl_->set_handler->getFeaturePool());
    cv.train(normalizer);
    cv.postIterationProcessing(train_scores, sanity.get());
    impl_->last_used_default_direction = sanity->used_default_direction;
    // getAvgWeights returns un-normalized (raw feature space) averaged weights:
    // size = n_features + 1, bias last. This is what score() will use.
    cv.getAvgWeights(avg_weights, normalizer);
  }
  catch (const std::exception& e)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      std::string("Percolator training failed: ") + e.what(), "");
  }

  // Mirror into getSvmWeights() buffer for backward compatibility.
  impl_->svm_weights.clear();
  if (!avg_weights.empty()) impl_->svm_weights.push_back(avg_weights);

  PercolatorModel model;
  model.format_version  = 1;
  model.feature_names   = input.feature_names;
  if (model.feature_names.empty())
  {
    for (size_t j = 0; j < n_feat; ++j)
    {
      model.feature_names.push_back("feat_" + std::to_string(j));
    }
  }
  model.weights         = std::move(avg_weights);
  model.normalizer_type = impl_->normalizer;
  model.seed            = impl_->seed;
  return model;
}

RescoreOutput Percolator::score(const RescoreInput& input,
                                const PercolatorModel& model)
{
  // Mark pi0 invalid up front: any throw — including from input/model
  // validation below or anywhere inside the scoring try block — leaves
  // last_pi0 in a "no usable value from this call" state, instead of stale
  // from a prior successful call.
  impl_->last_pi0 = -1.0;

  validateRescoreInput_(input);
  validateModelCompat_(input, model);

  impl_->reset();

  const size_t n_rows = input.rowCount();
  const size_t n_feat = input.featureCount();

  // Fresh feature pool. Reuses the second_pass_* scratch fields in Impl so
  // the next reset() cleans everything up via RAII.
  impl_->second_pass_pool = std::make_unique<P::FeatureMemoryPool>();
  impl_->second_pass_pool->createPool(static_cast<size_t>(n_feat));
  P::FeatureMemoryPool& pool = *impl_->second_pass_pool;

  // scoreAndAddPSM reads FeatureNames::getNumFeatures() internally. Register
  // the names we know so that static state is consistent with this scoring.
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

  P::Scores all_scores(/*usePi0=*/impl_->use_pi0);
  impl_->synthesized_ids.resize(n_rows);
  impl_->second_pass_psms.reserve(n_rows);

  for (size_t i = 0; i < n_rows; ++i)
  {
    const RowKeys keys = makeRowKeys(input, i);
    auto* psm = makePSMFromRow(input.row(i), n_feat, keys, i,
                               impl_->synthesized_ids[i], pool);
    impl_->second_pass_psms.push_back(psm);

    P::ScoreHolder sh;
    sh.pPSM  = psm;
    sh.label = input.is_decoy[i] ? P::LabelType::DECOY : P::LabelType::TARGET;
    // score = Σ feat[j] * model.weights[j] + model.weights[n_feat]
    // Raw features × raw weights — no normalizer at predict time.
    all_scores.scoreAndAddPSM(sh, model.weights, pool);
  }

  // Scoring pipeline: any of these can throw MyException on pathological
  // distributions (all-decoy, too-good-separation, median-decoy >= fdr-score).
  // Wrap and rebrand as Exception::InvalidValue for API consistency with
  // train()/rescore() — otherwise callers see raw upstream exceptions.
  try
  {
    all_scores.postMergeStep();
    all_scores.calcQvals(impl_->test_fdr);
    // normalizeScores takes a non-const ref and rescales the weight vector in
    // place (endScoreNormalizeWeights). Use a local copy so the caller's model
    // is not mutated.
    std::vector<double> weights_copy = model.weights;
    all_scores.normalizeScores(impl_->test_fdr, weights_copy);

    if (impl_->post_processing_tdc)
    {
      all_scores.weedOutRedundantTDC();
    }

    // Mirror Caller::calculatePSMProb: after any TDC weed-out, recompute
    // q-values so they use the post-TDC pi0 set by weedOutRedundantTDC's
    // internal postMergeStep. Kept unconditional for symmetry with rescore()
    // — when TDC is off this is a no-op (pi0 and scores are unchanged since
    // the earlier calcQvals call) but the cost is negligible on already-
    // sorted scores.
    all_scores.calcQvals(impl_->test_fdr);
    impl_->last_pi0 = all_scores.getPi0();

    if (impl_->pep_method == "isotonic")
    {
      all_scores.calcPep(/*spline=*/false, /*interp=*/false, /*pava=*/true);
    }
    else if (impl_->pep_method == "nonparametric")
    {
      all_scores.calcPep(/*spline=*/true,  /*interp=*/false, /*pava=*/false);
    }
    else
    {
      all_scores.calcPep(/*spline=*/false, /*interp=*/false, /*pava=*/false);
    }
  }
  catch (const std::exception& e)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      std::string("Percolator scoring failed: ") + e.what(), "");
  }

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
    if (id.size() < 12 || id.compare(0, 4, "row_") != 0) continue;
    char* end = nullptr;
    const char* num_begin = id.c_str() + 4;
    const unsigned long row = std::strtoul(num_begin, &end, 10);
    // Reject ids whose digit suffix didn't parse cleanly: no digits consumed
    // or unexpected trailing chars after the number.
    if (end == nullptr || end == num_begin || *end != '\0') continue;
    if (row >= n_rows) continue;
    out.scores[row]   = sh.score;
    out.q_values[row] = sh.q;
    out.peps[row]     = sh.pep;
    ++matched;
    if (sh.score < score_min) score_min = sh.score;
    if (sh.score > score_max) score_max = sh.score;
  }
  OPENMS_LOG_DEBUG << "[Percolator::score] extracted " << matched << " / " << n_rows
                   << " rows from all_scores (size=" << all_scores.size() << "); "
                   << "score range [" << score_min << ", " << score_max << "]"
                   << std::endl;

  return out;
}

RescoreOutput Percolator::rescore(const RescoreInput& input)
{
  // Mark pi0 invalid up front: any throw — including from input validation
  // below or anywhere inside the training/scoring try block — leaves
  // last_pi0 in a "no usable value from this call" state, instead of stale
  // from a prior successful call.
  impl_->last_pi0 = -1.0;
  impl_->last_used_default_direction = false;

  validateRescoreInput_(input);

  const size_t n_rows = input.rowCount();
  const size_t n_feat = input.featureCount();

  impl_->reset();
  P::PseudoRandom::setSeed(static_cast<unsigned long>(impl_->seed));

  // Reservoir-sample the training subset. When subset_max_train is 0 or
  // >= n_rows, training_rows is the identity [0..n).
  const std::vector<size_t> training_rows =
    selectTrainingSubset(input, n_rows, impl_->subset_max_train);
  const bool did_subsample = training_rows.size() < n_rows;

  impl_->set_handler = std::make_unique<P::SetHandler>(0);
  P::FeatureMemoryPool& pool = impl_->set_handler->getFeaturePool();
  pool.createPool(static_cast<size_t>(n_feat));

  // Register feature names. Sets static FeatureNames::numFeatures used by
  // Scores::scoreAndAddPSM in the second-pass branch below.
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

  impl_->target_ds = std::make_unique<P::DataSet>();
  impl_->decoy_ds  = std::make_unique<P::DataSet>();
  impl_->target_ds->setLabel(P::LabelType::TARGET);
  impl_->decoy_ds->setLabel(P::LabelType::DECOY);

  // Populate rows (training subset only — the second pass below scores the
  // full set using the trained model, mirroring upstream Caller.cpp flow).
  // The PSMs and their feature rows are allocated and registered by this
  // thread, in the order of the rows, and filled in parallel. (The data sets
  // delete the PSMs in this thread; PSMs allocated by other threads would take
  // much longer to free.)
  impl_->synthesized_ids.resize(n_rows);
  impl_->owned_psms.reserve(training_rows.size());
  for (size_t i : training_rows)
  {
    auto* psm = new P::PSMDescription();
    psm->features = pool.allocate();
    impl_->owned_psms.push_back(psm);
    if (input.is_decoy[i]) impl_->decoy_ds->registerPsm(psm);
    else                   impl_->target_ds->registerPsm(psm);
  }
  std::vector<std::exception_ptr> row_errors(training_rows.size());
#pragma omp parallel for schedule(dynamic, 1024)
  for (SignedSize k = 0; k < static_cast<SignedSize>(training_rows.size()); ++k)
  {
    try
    {
      const size_t i = training_rows[k];
      fillPSMFromRow(*impl_->owned_psms[k], input.row(i), n_feat, makeRowKeys(input, i), i,
                     impl_->synthesized_ids[i]);
    }
    catch (...)
    {
      row_errors[k] = std::current_exception();
    }
  }
  rethrowFirstError(row_errors);

  impl_->set_handler->push_back_dataset(impl_->target_ds.release());
  impl_->set_handler->push_back_dataset(impl_->decoy_ds.release());

  int norm_type = P::Normalizer::STDV;
  if      (impl_->normalizer == "uni")  norm_type = P::Normalizer::UNI;
  else if (impl_->normalizer == "none") norm_type = P::Normalizer::NONORM;
  P::Globals::getInstance()->setNoTerminate(true);
  P::Normalizer::resetNormalizer();
  P::Normalizer::setType(norm_type);
  // normalizeFeatures() owns Normalizer construction via its out-parameter.
  P::Normalizer* normalizer = nullptr;

  P::SanityCheck::setInitDefaultDir(0);
  P::SanityCheck::setInitDefaultDirName(impl_->initial_direction);
  // The plain SanityCheck that SanityCheck::initialize() returns, recording
  // whether it replaced the learned weights (see usedDefaultDirection()).
  // unique_ptr releases it on every exit path (including throws).
  auto sanity = std::make_unique<RecordingSanityCheck>();
  try
  {
    sanity->checkAndSetDefaultDir();
  }
  catch (const std::exception& e)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      std::string("initial_direction not resolved: ") + e.what(),
      impl_->initial_direction);
  }

  impl_->set_handler->normalizeFeatures(normalizer);

  P::Scores all_scores(/*usePi0=*/impl_->use_pi0);
  all_scores.populateWithPSMs(*impl_->set_handler);

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
    /*numThreads=*/static_cast<unsigned int>(impl_->num_threads),
    /*skipNormalizeScores=*/false,
    /*decoyFractionTraining=*/1.0,
    /*numFolds=*/3);

  std::vector<double> avg_weights;
  try
  {
    // CV uses sanity transiently per call (no storage); the unique_ptr in
    // this scope outlives cv, so the bare pointer stays valid.
    cv.preIterationSetup(all_scores, sanity.get(), normalizer,
                         impl_->set_handler->getFeaturePool());
    cv.train(normalizer);
    cv.postIterationProcessing(all_scores, sanity.get());
    impl_->last_used_default_direction = sanity->used_default_direction;
    cv.getAvgWeights(avg_weights, normalizer);
  }
  catch (const std::exception& e)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      std::string("Percolator training failed: ") + e.what(), "");
  }

  impl_->svm_weights.clear();
  if (!avg_weights.empty()) impl_->svm_weights.push_back(avg_weights);

  try
  {
    // Second pass: when we subsampled for training, score the FULL input set
    // using the trained weights. Mirrors upstream Caller.cpp (rel-3-08-01
    // lines 1196-1236): reset allScores, rebuild via readAndScoreTab, then
    // postMergeStep + calcQvals + normalizeScores. Makes the subset_max_train
    // docstring claim "trained model applied to ALL PSMs" true.
    if (did_subsample)
    {
      impl_->second_pass_pool = std::make_unique<P::FeatureMemoryPool>();
      impl_->second_pass_pool->createPool(static_cast<size_t>(n_feat));
      P::FeatureMemoryPool& full_pool = *impl_->second_pass_pool;

      all_scores = P::Scores(impl_->use_pi0);
      impl_->second_pass_psms.reserve(n_rows);

      for (size_t i = 0; i < n_rows; ++i)
      {
        const RowKeys keys = makeRowKeys(input, i);
        auto* psm = makePSMFromRow(input.row(i), n_feat, keys, i,
                                   impl_->synthesized_ids[i], full_pool);
        impl_->second_pass_psms.push_back(psm);

        P::ScoreHolder sh;
        sh.pPSM  = psm;
        sh.label = input.is_decoy[i] ? P::LabelType::DECOY : P::LabelType::TARGET;
        all_scores.scoreAndAddPSM(sh, avg_weights, full_pool);
      }

      all_scores.postMergeStep();
      all_scores.calcQvals(impl_->test_fdr);
      all_scores.normalizeScores(impl_->test_fdr, avg_weights);
    }

    if (impl_->post_processing_tdc)
    {
      // In a parallel region, as the cross-validation calls Scores::merge: the sorts and the
      // bootstraps of the pi0 estimate run as its tasks, with the same results (see
      // thirdparty/percolator/SortLikeStd.h and PosteriorEstimator::estimatePi0).
      std::exception_ptr tdc_error;
#pragma omp parallel num_threads(teamThreads(impl_->num_threads))
#pragma omp single
      {
        try
        {
          all_scores.weedOutRedundantTDC();
        }
        catch (...)
        {
          tdc_error = std::current_exception();
        }
      }
      if (tdc_error) std::rethrow_exception(tdc_error);
    }

    // Mirror Caller::calculatePSMProb: after the merge, recompute q-values
    // on the merged set so they use the merged-set pi0 (set by postMergeStep
    // inside merge, and re-set by weedOutRedundantTDC's internal
    // postMergeStep when TDC ran), not the stale per-fold pi0 left behind by
    // Scores::merge's per-fold calcQvals loop.
    all_scores.calcQvals(impl_->test_fdr);
    impl_->last_pi0 = all_scores.getPi0();

    if (impl_->pep_method == "isotonic")
    {
      all_scores.calcPep(/*spline=*/false, /*interp=*/false, /*pava=*/true);
    }
    else if (impl_->pep_method == "nonparametric")
    {
      all_scores.calcPep(/*spline=*/true,  /*interp=*/false, /*pava=*/false);
    }
    else
    {
      all_scores.calcPep(/*spline=*/false, /*interp=*/false, /*pava=*/false);
    }
  }
  catch (const std::exception& e)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      std::string("Percolator scoring failed: ") + e.what(), "");
  }

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
    if (id.size() < 12 || id.compare(0, 4, "row_") != 0) continue;
    char* end = nullptr;
    const char* num_begin = id.c_str() + 4;
    const unsigned long row = std::strtoul(num_begin, &end, 10);
    // Reject ids whose digit suffix didn't parse cleanly: no digits consumed
    // or unexpected trailing chars after the number.
    if (end == nullptr || end == num_begin || *end != '\0') continue;
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
  StringList discoverFeatureNames_(const PeptideHit& hit, const std::string& td_meta)
  {
    static const std::vector<std::string> blocklist {
      "percolator_score", "percolator_q_value", "percolator_pep",
      "q-value", "q-value_score", "score", "RT", "MZ",
      "target_decoy", "isDecoy"
    };
    StringList out;
    std::vector<std::string> keys;
    hit.getKeys(keys);
    for (const std::string& k : keys)
    {
      if (k == td_meta) continue;
      bool skipped = false;
      for (const std::string& b : blocklist)
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
  const std::string first_sid = PercolatorInfile::getScanIdentifier(peptide_ids.front(), 0);
  const RegularExpression scan_regex(SpectrumLookup::getRegExFromNativeID(first_sid));

  size_t row = 0;
  size_t pid_index = 0;
  for (const auto& pid : peptide_ids)
  {
    ++pid_index;
    if (pid.getHits().empty()) continue;

    const std::string scan_identifier = PercolatorInfile::getScanIdentifier(pid, pid_index);
    const Int scan_number = SpectrumLookup::extractScanNumber(
        scan_identifier, scan_regex, /*no_error=*/true);

    // the spectrum file, numbered in order of appearance as Percolator numbers the FileName column of a .pin file
    const std::string file_key = PercolatorInfile::getFileIdentifier(pid);

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
        calc_mass = PercolatorInfile::getFeatureValue(hit.getMetaValue("CalcMass"), "CalcMass");
      }
      else
      {
        // Mirror PercolatorInfile::store fallback (PercolatorInfile.cpp:460):
        // m/z at the hit's charge, not the neutral mass. This matters for
        // CV-fold parity with the subprocess path when CalcMass has not
        // already been stamped on the hit.
        try { calc_mass = hit.getSequence().getMZ(hit.getCharge()); }
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

void Percolator::saveModel(const PercolatorModel& model, const std::string& filename)
{
  if (model.weights.size() != model.feature_names.size() + 1)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "PercolatorModel invariant broken: weights.size() != feature_names.size()+1",
      std::to_string(model.weights.size()) + " vs " +
      std::to_string(model.feature_names.size() + 1));
  }
  std::ofstream os(filename.c_str());
  if (!os)
  {
    throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      filename);
  }
  os.precision(std::numeric_limits<double>::max_digits10);
  os << "# OpenMS Percolator model file\n";
  os << "format_version: " << model.format_version << '\n';
  os << "normalizer: "     << model.normalizer_type << '\n';
  os << "seed: "           << model.seed << '\n';
  os << "n_features: "     << model.feature_names.size() << '\n';
  os << "bias: "           << model.weights.back() << '\n';
  os << "# features — one 'name<TAB>weight' per line\n";
  for (size_t j = 0; j < model.feature_names.size(); ++j)
  {
    os << static_cast<std::string>(model.feature_names[j])
       << '\t' << model.weights[j] << '\n';
  }
}

PercolatorModel Percolator::loadModel(const std::string& filename)
{
  std::ifstream is(filename.c_str());
  if (!is)
  {
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
  }
  PercolatorModel model;
  model.format_version = 0;

  // Track which required header keys have been observed. Duplicates are
  // errors (not silent last-wins).
  std::unordered_set<std::string> seen_keys;
  long declared_n_features = -1;
  double bias = 0.0;
  bool bias_seen = false;
  static const std::unordered_set<std::string> known_keys{
    "format_version", "normalizer", "seed", "n_features", "bias"
  };
  static const std::unordered_set<std::string> known_normalizers{
    "stdv", "uni", "none"
  };

  auto starts_with = [](const std::string& s, const std::string& pfx)
  {
    return s.size() >= pfx.size() && s.compare(0, pfx.size(), pfx) == 0;
  };
  auto trim = [](std::string s)
  {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return std::string{};
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
  };

  std::string line;
  while (std::getline(is, line))
  {
    if (line.empty() || starts_with(line, "#")) continue;

    // Header key/value lines have no TAB; data lines have a TAB separator
    // between feature name and weight. A feature name cannot contain a ':'
    // without a TAB (names are written as the left-hand side of TAB-split
    // data rows), so checking for TAB first is unambiguous.
    if (!line.contains('\t'))
    {
      const size_t colon = line.find(':');
      if (colon == std::string::npos)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          filename + ": expected 'key: value' or 'name<TAB>weight', got '" +
          line + "'", line);
      }
      const std::string key = trim(line.substr(0, colon));
      const std::string val = trim(line.substr(colon + 1));
      if (!known_keys.contains(key))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          filename + ": unknown header key '" + key + "'", line);
      }
      if (!seen_keys.insert(key).second)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          filename + ": duplicate header key '" + key + "'", line);
      }
      try
      {
        if      (key == "format_version") model.format_version = std::stoi(val);
        else if (key == "seed")           model.seed = std::stoi(val);
        else if (key == "n_features")     declared_n_features = std::stol(val);
        else if (key == "bias")         { bias = std::stod(val); bias_seen = true; }
        else if (key == "normalizer")
        {
          if (!known_normalizers.contains(val))
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              filename + ": invalid normalizer value (expected stdv|uni|none)",
              val);
          }
          model.normalizer_type = val;
        }
      }
      catch (const std::invalid_argument&)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          filename + ": cannot parse value for '" + key + "'", val);
      }
      continue;
    }

    // Data row: "<feature_name>\t<weight>"
    const size_t tab = line.find('\t');
    const std::string name = line.substr(0, tab);
    double w;
    try { w = std::stod(line.substr(tab + 1)); }
    catch (const std::exception&)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        filename + ": cannot parse weight on line '" + line + "'", line);
    }
    model.feature_names.push_back(name);
    model.weights.push_back(w);
  }

  // Required-key completeness.
  for (const auto& required : known_keys)
  {
    if (!seen_keys.contains(required))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        filename + ": missing required header key '" + required + "'", "");
    }
  }
  if (model.format_version != 1)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      filename + ": unsupported format_version (expected 1)",
      std::to_string(model.format_version));
  }
  if (declared_n_features < 0 ||
      static_cast<size_t>(declared_n_features) != model.feature_names.size())
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      filename + ": declared n_features does not match feature row count",
      std::to_string(declared_n_features) + " vs " +
      std::to_string(model.feature_names.size()));
  }
  // Append bias to weights so the final layout matches what train() produces:
  // weights.size() == feature_names.size() + 1, bias last.
  if (!bias_seen)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      filename + ": 'bias' header key not found", "");
  }
  model.weights.push_back(bias);
  return model;
}

void Percolator::rescore(std::vector<PeptideIdentification>& peptide_ids,
                         const StringList& feature_names)
{
  if (peptide_ids.empty()) return;

  const std::string td_meta = param_.getValue("target_decoy_metavalue").toString();

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
      for (const std::string& f : features)
      {
        if (!hit.metaValueExists(f))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "PeptideHit missing feature meta value", f);
        }
        row.push_back(PercolatorInfile::getFeatureValue(hit.getMetaValue(f), f));
      }
      ri.features.push_back(std::move(row));
      ri.is_decoy.push_back(hit.getMetaValue(td_meta).toString() == "decoy");

      // CV group key: the pid index. All hits from peptide_ids[i] share i,
      // which is the only invariant CV grouping needs (rows in the same
      // group land in the same fold). Equivalent to the previous
      // hash(identifier+RT) but without string allocations or hash
      // collisions in the hot loop.
      ri.cv_group_keys.push_back(static_cast<int>(i));

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

void Percolator::buildPSMRescoreInput(const PeptideIdentificationList& peptide_ids,
                                      const StringList& feature_set,
                                      const std::string& enzyme,
                                      int min_charge,
                                      int max_charge,
                                      RescoreInput& input,
                                      PSMRowLocations& rows)
{
  // SpecId/Label/Peptide/Proteins/ScanNr are strings or bookkeeping, not training features. ExpMass is
  // used by Percolator for the spectrum key (it goes into RescoreInput.exp_masses), not for training.
  const std::set<std::string> pin_metadata_not_feature {
    "SpecId", "Label", "ScanNr", "ExpMass", "Peptide", "Proteins"
  };
  StringList numeric_features;
  for (const std::string& f : feature_set)
  {
    if (pin_metadata_not_feature.contains(f)) continue;
    numeric_features.push_back(f);
  }
  // (computed first: errors of the PSMs come before a missing feature set, as when the PIN features were stamped)
  const Size skipped = PercolatorInfile::buildFeatureMatrix(peptide_ids, numeric_features, enzyme, min_charge, max_charge,
                                                            input, rows);
  OPENMS_LOG_INFO << "Computed the PIN features; " << skipped
                  << " PSMs skipped (no evidence or unknown TD)." << std::endl;
  if (numeric_features.empty())
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "No numeric Percolator features in the feature set.");
  }
  if (rows.empty())
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "No PSM carries all Percolator features.");
  }
}

Size Percolator::rescorePSMs(PeptideIdentificationList& peptide_ids,
                             const StringList& feature_set,
                             const std::string& enzyme,
                             int min_charge,
                             int max_charge,
                             const std::string& score_type,
                             bool best_hit_only)
{
  if (score_type != "q-value" && score_type != "pep" && score_type != "svm")
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "score_type must be 'q-value', 'pep' or 'svm'", score_type);
  }
  if (std::none_of(peptide_ids.begin(), peptide_ids.end(),
                   [](const PeptideIdentification& pid) { return !pid.getHits().empty(); }))
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "No PSMs to rescore.");
  }

  // Training then uses the feature vectors the executable would read from the .pin file, row for row.
  RescoreInput ri;
  PSMRowLocations hit_locs;
  buildPSMRescoreInput(peptide_ids, feature_set, enzyme, min_charge, max_charge, ri, hit_locs);

  // A hit that carries a meta value under the name of a PIN column gets its PIN value, as stamping did (before the
  // training, so also if it fails). After the matrix: its computation succeeded for every hit, the matrix's errors
  // are raised in order first.
  PercolatorInfile::refreshPinMetaValues(peptide_ids, enzyme, min_charge, max_charge);
  OPENMS_LOG_INFO << "Rescoring " << ri.rowCount() << " PSMs in-process with "
                  << ri.featureCount() << " features." << std::endl;

  // Low-level API: the high-level rescore(peptide_ids, ...) ignores the PIN-compatible fields.
  const RescoreOutput ro = rescore(ri);

  // When the PIN features were stamped, the PIN columns that a stamped hit (evidence, target/decoy type) did not have
  // were removed after the write-back below, also the original score that it stores under the old score type if that is
  // the name of a PIN column (e.g. "score"). Kept that way: these hits lose that meta value again.
  // (identification, hit, old score type), removed also if the write-back throws, as by the cleanup of the stamped values
  struct RemoveOldScore
  {
    PeptideIdentificationList& ids;
    std::vector<std::tuple<size_t, size_t, std::string>> hits;
    ~RemoveOldScore()
    {
      for (const auto& [i, j, type] : hits) ids[i].getHits()[j].removeMetaValue(type); // (allocates nothing)
    }
  } remove_old_score{peptide_ids, {}};
  {
    const StringList pin_columns = PercolatorInfile::getPinColumnNames(min_charge, max_charge);
    const std::string* last_type = nullptr;
    bool pin_type = false;
    for (size_t i = 0; i < peptide_ids.size(); ++i)
    {
      const std::string& type = peptide_ids[i].getScoreType();
      if (last_type == nullptr || type != *last_type)
      {
        last_type = &type;
        pin_type = std::find(pin_columns.begin(), pin_columns.end(), type) != pin_columns.end();
      }
      if (!pin_type) continue;
      const std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      for (size_t j = 0; j < hits.size(); ++j)
      {
        if (!hits[j].getPeptideEvidences().empty() && hits[j].getTargetDecoyType() != PeptideHit::TargetDecoyType::UNKNOWN
            && !hits[j].metaValueExists(type))
        {
          remove_old_score.hits.emplace_back(i, j, type);
        }
      }
    }
  }

  // looking names up takes a lock on the registry, which threads would contend for: the parallel loops below use indices
  auto index_of = [](const std::string& name) { return MetaInfoInterface::metaRegistry().getIndex(name); };
  constexpr size_t NO_ROW = std::numeric_limits<size_t>::max();

  // The row of each hit (NO_ROW: not rescored); the Percolator scores of a hit are read from the rows, not back from
  // the meta values they were written to.
  std::vector<std::vector<size_t>> row_of(peptide_ids.size());
  for (size_t i = 0; i < peptide_ids.size(); ++i)
  {
    row_of[i].assign(peptide_ids[i].getHits().size(), NO_ROW);
  }
  for (size_t row = 0; row < ro.scores.size(); ++row)
  {
    row_of[hit_locs[row].first][hit_locs[row].second] = row;
  }

  // The identifications without hits come first (updated serially), then the first one with hits; with several score
  // types (each a meta value name) all spectra are updated serially, registering their names in order.
  size_t first_with_hits = 0;
  while (first_with_hits < peptide_ids.size() && peptide_ids[first_with_hits].getHits().empty()) ++first_with_hits;
  const bool one_score_type = first_with_hits < peptide_ids.size()
    && std::all_of(peptide_ids.begin() + first_with_hits + 1, peptide_ids.end(),
                   [&](const PeptideIdentification& pid) { return pid.getHits().empty() || pid.getScoreType() == peptide_ids[first_with_hits].getScoreType(); });
  // Only the best hits get meta values: with one score type (the names are registered below as the writes to all hits
  // registered them).
  const bool best_only = best_hit_only && one_score_type;

  // The hits are updated in parallel, except for the first row and the first spectrum with hits: these register the
  // new meta value names in the order of a serial loop, and MetaInfo keeps meta values in the order of their name
  // indices. The parallel updates use the name indices: looking names up takes a lock on the registry, which the
  // threads would contend for.
  auto store_row = [&](size_t row, const auto& score_key, const auto& q_value_key, const auto& pep_key)
  {
    const auto [pid_i, hit_i] = hit_locs[row];
    PeptideHit& hit = peptide_ids[pid_i].getHits()[hit_i];
    hit.setMetaValue(score_key,   ro.scores[row]);
    hit.setMetaValue(q_value_key, ro.q_values[row]);
    hit.setMetaValue(pep_key,     ro.peps[row]);
  };
  if (!ro.scores.empty())
  {
    if (best_only)
    {
      // (as the write to the first row registered them)
      MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
      for (const char* name : {"percolator_score", "percolator_q_value", "percolator_pep"}) registry.registerName(name);
    }
    else
    {
      store_row(0, "percolator_score", "percolator_q_value", "percolator_pep");
      const UInt score_index = index_of("percolator_score"), q_value_index = index_of("percolator_q_value"),
                 pep_index = index_of("percolator_pep");
      std::vector<std::exception_ptr> errors(ro.scores.size());
#pragma omp parallel for schedule(dynamic, 256)
      for (SignedSize row = 1; row < static_cast<SignedSize>(ro.scores.size()); ++row)
      {
        try
        {
          store_row(static_cast<size_t>(row), score_index, q_value_index, pep_index);
        }
        catch (...)
        {
          errors[row] = std::current_exception();
        }
      }
      rethrowFirstError(errors);
    }
  }

  // Canonical score fields and PSI-MS CV terms, as PercolatorAdapter writes them for the executable.
  // index: the name indices of the old score type, percolator_score, _q_value, _pep, MS:1001492, 91 and 93, or empty
  auto store_canonical = [&](size_t i, const std::vector<UInt>& index)
  {
    PeptideIdentification& pid = peptide_ids[i];
    const std::string old_score_type = pid.getScoreType();
    if (score_type == "pep")
    {
      pid.setScoreType("Posterior Error Probability");
      pid.setHigherScoreBetter(false);
    }
    else if (score_type == "svm")
    {
      pid.setScoreType("svm");
      pid.setHigherScoreBetter(true);
    }
    else
    {
      pid.setScoreType("q-value");
      pid.setHigherScoreBetter(false);
    }

    auto& hits = pid.getHits();
    // Hits that were not rescored get the values the executable path gives PSMs missing from
    // Percolator's output (q-value = PEP = 1, SVM score = -100: confidently not identified), so
    // their main score matches the new score type of the identification.
    // (an old score type named like a Percolator value: the original score, kept under that name after the Percolator
    // values were written, replaced the value there, and the values were read back from the meta values)
    const bool old_is_score = old_score_type == "percolator_score", old_is_q_value = old_score_type == "percolator_q_value",
               old_is_pep = old_score_type == "percolator_pep";
    auto scores_of = [&](size_t j, double& svm, double& qval, double& pep)
    {
      svm = -100.0; qval = 1.0; pep = 1.0;
      const size_t row = row_of[i][j];
      if (row != NO_ROW)
      {
        svm  = old_is_score   ? hits[j].getScore() : ro.scores[row];
        qval = old_is_q_value ? hits[j].getScore() : ro.q_values[row];
        pep  = old_is_pep     ? hits[j].getScore() : ro.peps[row];
      }
    };
    auto new_score = [&](double svm, double qval, double pep)
    {
      if (score_type == "q-value") return qval;
      if (score_type == "pep") return pep;
      return svm;
    };
    // best_only: the hit that the stable sort by the new score puts first (the first of the lowest q-values or PEPs,
    // highest SVM scores); with a NaN score that order is not defined: all hits are written
    size_t best = NO_ROW;
    if (best_only && !index.empty() && !hits.empty())
    {
      bool nan = false;
      double best_score = 0.0;
      for (size_t j = 0; j < hits.size() && !nan; ++j)
      {
        double svm, qval, pep;
        scores_of(j, svm, qval, pep);
        const double score = new_score(svm, qval, pep);
        if (std::isnan(score)) { nan = true; break; }
        if (j == 0 || (score_type == "svm" ? score > best_score : score < best_score))
        {
          best = j;
          best_score = score;
        }
      }
      if (nan)
      {
        best = NO_ROW;
        for (size_t j = 0; j < hits.size(); ++j)
        {
          if (row_of[i][j] != NO_ROW) store_row(row_of[i][j], index[1], index[2], index[3]);
        }
      }
    }
    for (size_t j = 0; j < hits.size(); ++j)
    {
      PeptideHit& hit = hits[j];
      double svm, qval, pep;
      scores_of(j, svm, qval, pep);
      if (index.empty())
      {
        hit.setMetaValue(old_score_type, hit.getScore());  // preserve original
        hit.setMetaValue("MS:1001492", svm);    // percolator:score
        hit.setMetaValue("MS:1001491", qval);   // percolator:Q value
        hit.setMetaValue("MS:1001493", pep);    // percolator:PEP
      }
      else if (best == NO_ROW || j == best)
      {
        // (in the order of the writes to all hits: the Percolator values first, then the original score)
        if (best != NO_ROW && row_of[i][j] != NO_ROW) store_row(row_of[i][j], index[1], index[2], index[3]);
        hit.setMetaValue(index[0], hit.getScore());
        hit.setMetaValue(index[4], svm);
        hit.setMetaValue(index[5], qval);
        hit.setMetaValue(index[6], pep);
      }
      hit.setScore(new_score(svm, qval, pep));
    }
  };
  size_t first_parallel = 0;
  while (first_parallel < first_with_hits)
  {
    store_canonical(first_parallel++, {});
  }
  if (first_parallel < peptide_ids.size())
  {
    const std::string first_score_type = peptide_ids[first_parallel].getScoreType();
    std::vector<UInt> index;
    if (best_only)
    {
      // the names the serial update of the first spectrum with hits registered, in that order; it is updated in parallel
      MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
      for (const std::string& name : {first_score_type, std::string("MS:1001492"), std::string("MS:1001491"), std::string("MS:1001493")})
      {
        registry.registerName(name);
      }
    }
    else
    {
      store_canonical(first_parallel++, {});
    }
    if (one_score_type)
    {
      for (const std::string& name : {first_score_type, std::string("percolator_score"), std::string("percolator_q_value"),
                                      std::string("percolator_pep"), std::string("MS:1001492"), std::string("MS:1001491"),
                                      std::string("MS:1001493")})
      {
        index.push_back(index_of(name));
      }
    }
    std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(dynamic, 256) if(one_score_type)
    for (SignedSize i = static_cast<SignedSize>(first_parallel); i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      try
      {
        store_canonical(static_cast<size_t>(i), index);
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    rethrowFirstError(errors);
  }
  return ro.scores.size();
}

} // namespace OpenMS
