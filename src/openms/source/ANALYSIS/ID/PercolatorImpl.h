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
    std::vector<OpenMS::Internal::Percolator::PSMDescription*> owned_psms;

    // Persisted params (mirrored from DefaultParamHandler param_).
    double c_pos = 0.1;
    double c_neg = 0.0;          ///< 0 = auto
    double test_fdr = 0.01;
    double train_fdr = 0.01;
    int    num_iterations = 10;
    int    seed = 1;
    std::string pep_method = "logistic_regression";

    void reset();
  };
}
