// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include <ranges>

#include "PeakTypes.h"

namespace OpenMS {
namespace PipEcho {

  class Pep {
  public:
    Pep(const std::vector<Acceptor*>);

    /**
     * Run the PEP analysis.
     *
     * Returns `true` if the analysis ran successfully.
     */
    bool run();

  public:
    /// Everything we need to know about an acceptor feature.
    struct acceptor_t {
      acceptor_t(Acceptor* acceptor, Acceptor::scored_t& scored, bool is_target)
      : acceptor(acceptor),
        donor(scored.second),
        score(scored.first),
        is_target(is_target)
      {};

      Acceptor* acceptor;
      const Donor* donor;
      Score score;
      double pep_score = 1.0;
      bool is_target = false;
    };

    /// Just for cleaner looking code.
    using acceptor_ptr_t = std::shared_ptr<acceptor_t>;

    /// A group of `acceptor_t` objects.
    using group_t = std::vector<acceptor_ptr_t>;

  private:
    /// Group the acceptors into roughly equal sized buckets.
    void group_acceptors(std::vector<group_t>&);

    /// One round of training and predicting for the given groups.
    bool round(std::vector<group_t>&, std::function<double(const acceptor_t&)>);

    /// Train a model using the data in the first view, and generate
    /// predictions for the data in the second parameter..
    bool train_predict(auto&& training,
                       group_t& predict,
                       double cutoff,
                       std::function<double(const acceptor_t&)> get);

    /// Remove duplicate acceptors after they all PEP values.
    //void remove_duplicates();

    /// Compute FDR, reducing the acceptors to those that are above
    /// the cutoff threshold.
    //void perform_fdr(double);

  private:
    group_t acceptors;

  };
}}
