// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// Internal pimpl header for OpenMS::Percolator. Not installed publicly.

#pragma once

#include <OpenMS/ANALYSIS/ID/Percolator.h>

// Vendored Percolator headers. The thirdparty/percolator include directory is
// wired PRIVATE to libOpenMS, so these are only visible here.
#include "DataSet.h"
#include "SetHandler.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "CrossValidation.h"
#include "PSMDescription.h"
#include "FeatureNames.h"
#include "PseudoRandom.h"
// OpenMS-local: reservoir-sampling helper derived from SetHandler::readPSMs.
// Needed because we don't vendor Caller.cpp (which orchestrates upstream's
// subset-train + full-set-score flow). See InProcessSampler.h for upstream
// line references and re-sync instructions.
#include "InProcessSampler.h"

#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{
  struct Percolator::Impl
  {
    // Per-call scratch. Reset between rescore() invocations.
    // Percolator's vendored code is designed for one shot per instance.
    std::unique_ptr<OpenMS::Internal::Percolator::SetHandler> set_handler;
    std::unique_ptr<OpenMS::Internal::Percolator::DataSet> target_ds;
    std::unique_ptr<OpenMS::Internal::Percolator::DataSet> decoy_ds;
    std::vector<std::string> synthesized_ids;
    /// PSMs are owned by the DataSets inside set_handler; this vector is for
    /// bookkeeping only (so we can count rows). Do NOT delete these pointers
    /// directly — DataSet::~DataSet will.
    std::vector<OpenMS::Internal::Percolator::PSMDescription*> owned_psms;
    /// PSMs built for the "score the full set" second pass when subsampling is
    /// active. Not owned by any DataSet — reset() must delete them manually.
    std::vector<OpenMS::Internal::Percolator::PSMDescription*> second_pass_psms;
    /// Fresh FeatureMemoryPool used by the second pass (feature slots are
    /// deallocated by scoreAndAddPSM but the pool itself must outlive it).
    std::unique_ptr<OpenMS::Internal::Percolator::FeatureMemoryPool> second_pass_pool;

    // Persisted params (mirrored from DefaultParamHandler param_).
    double c_pos = 0.1;
    double c_neg = 0.0;          ///< 0 = auto
    double test_fdr = 0.01;
    double train_fdr = 0.01;
    int    num_iterations = 10;
    int    seed = 1;
    std::string pep_method = "logistic_regression";
    std::string normalizer = "stdv";          ///< "stdv" | "uni" | "none"
    bool   train_best_positive = true;
    bool   post_processing_tdc = true;
    int    nested_xval_bins = 1;
    int    subset_max_train = 100000;         ///< 0 = no cap; default 100k
    int    num_threads = 3;                   ///< OpenMP threads for CV (one per fold; max useful = 3)
    std::string initial_direction;            ///< feature name; prefix "-" = lower is better; empty = auto
    bool   use_pi0 = true;                    ///< pi0 correction on q-values; false = pure TDA (pi0=1)
    std::string report_as_main_score = "none"; ///< "none"|"q-value"|"pep"|"svm" — which Percolator output to promote to hit.getScore()

    /// Per-fold, per-feature weights populated after rescore(). Read by
    /// Percolator::getSvmWeights(). Empty before the first rescore().
    std::vector<std::vector<double>> svm_weights;

    /// Pi0 (null fraction) estimated on the post-CV, post-TDC merged Scores
    /// set right before calcPep. Populated by rescore() and score(). Read by
    /// Percolator::getPi0(). Negative (-1.0) before the first call.
    double last_pi0 = -1.0;

    void reset();
  };
}
