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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>

namespace OpenMS
{

namespace P = OpenMS::Internal::Percolator;

void Percolator::Impl::reset()
{
  // Delete owned PSMs before resetting containers. Percolator's DataSet holds
  // raw pointers via push_back but doesn't free them.
  for (auto* psm : owned_psms)
  {
    delete psm;
  }
  owned_psms.clear();
  synthesized_ids.clear();
  target_ds.reset();
  decoy_ds.reset();
  set_handler.reset();
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
  defaults_.setValue("initial_direction","",   "Feature name for initial scoring direction. Empty = auto.");
  defaults_.setValue("pep_method",       "logistic_regression",
                     "PEP estimator: 'logistic_regression' or 'isotonic'.");
  defaults_.setValidStrings("pep_method", {"logistic_regression","isotonic"});
  defaults_.setValue("seed",             1,    "Random seed for CV splitting.");
  defaults_.setValue("target_decoy_metavalue", "target_decoy",
                     "Meta-value on PeptideHit indicating target/decoy ('target'/'decoy').");
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
}

// Populate a single PSMDescription from one row of RescoreInput. The row's
// feature pointer comes from SetHandler's FeatureMemoryPool.
static P::PSMDescription* makePSMFromRow(
  const std::vector<double>& feats,
  int cv_key,
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
  psm->scan        = static_cast<unsigned int>(cv_key);   // opaque CV grouping key
  psm->specFileNr  = 0;
  psm->expMass     = 0.0;
  psm->calcMass    = 0.0;
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

  // Set up SetHandler + feature pool.
  impl_->set_handler = std::make_unique<P::SetHandler>(0 /* maxPSMs, 0 = no cap */);
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

  // Populate rows.
  impl_->synthesized_ids.resize(n_rows);
  impl_->owned_psms.reserve(n_rows);
  for (size_t i = 0; i < n_rows; ++i)
  {
    const int cv_key = input.cv_group_keys.empty()
      ? static_cast<int>(i)
      : input.cv_group_keys[i];
    auto* psm = makePSMFromRow(input.features[i], cv_key, i,
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
  P::Normalizer::setType(P::Normalizer::STDV);
  P::Normalizer* normalizer = P::Normalizer::getNormalizer();

  // SanityCheck: "default" fingerprint → generic PSM-level sanity.
  P::SanityCheck* sanity = P::SanityCheck::initialize("");
  if (!sanity)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Could not initialize Percolator SanityCheck", "");
  }

  impl_->set_handler->normalizeFeatures(normalizer);

  // Build the full score set and populate from SetHandler.
  P::Scores all_scores(/*usePi0=*/true);
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
    /*usePi0=*/true,
    /*nestedXvalBins=*/1,
    /*trainBestPositive=*/false,
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

  for (const auto& sh : all_scores)
  {
    const std::string& id = sh.pPSM->getFullPeptide();
    // Parse "row_NNNNNNNN" → row index
    if (id.size() < 12 || id.compare(0, 4, "row_") != 0) continue;
    char* end = nullptr;
    const unsigned long row = std::strtoul(id.c_str() + 4, &end, 10);
    if (row >= n_rows) continue;
    out.scores[row]   = sh.score;
    out.q_values[row] = sh.q;
    out.peps[row]     = sh.pep;
  }

  return out;
}

void Percolator::rescore(std::vector<PeptideIdentification>& /*peptide_ids*/,
                         const StringList& /*feature_names*/)
{
  // High-level PSM API: Phase B in the implementation plan. Not yet implemented.
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
}

} // namespace OpenMS
