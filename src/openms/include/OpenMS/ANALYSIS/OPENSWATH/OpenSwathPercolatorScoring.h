// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/PercolatorTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/config.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief In-process Percolator scoring for OpenSWATH OSW and OSWPQ results.

    The scorer reads OpenSWATH feature tables, converts them into a
    domain-agnostic Percolator::RescoreInput, performs target/decoy rescoring,
    and writes the resulting score, q-value, PEP, and peak-group rank back to
    the source format. Rows with missing active OpenSWATH features are skipped
    instead of being imputed to zero, matching PyProphet's cleanup step.

    Peak-group-level q-values / PEPs follow the OpenSWATH / PyProphet
    convention: Percolator provides the discriminant, top-ranked target/decoy
    peak groups define the error model, and the resulting statistics are then
    projected back onto all scored candidates.

    Supported levels:
    - MS1: writes SCORE_MS1 (SQLite) or score_ms1_* columns (.oswpq)
    - MS2: writes SCORE_MS2 (SQLite) or score_ms2_* columns (.oswpq)
    - MS1MS2: uses combined MS2 + MS1 features and writes SCORE_MS2 /
      score_ms2_* (matching the historical PyProphet convention)
    - Transition: requires prior MS2/MS1MS2 scoring and writes
      SCORE_TRANSITION / score_transition_*

    For MS1, MS2 and MS1MS2 scoring, Percolator is trained on a subset of the
    peak groups (see trainingRows()) and the model is then applied to all peak
    groups. The subset holds the best peak group of every precursor of every
    run by the parameter training_main_score, which brings in most of the rare
    true peak groups, and all peak groups of precursors drawn at random, up to
    percolator:subset_max_train further peak groups. Within these
    complete precursors, Percolator's choice of the best peak group per
    precursor in each training iteration (percolator:train_best_positive) can
    still promote a peak group that the main score ranked lower.
    percolator:subset_max_train = 0 trains on all peak groups. Transitions
    are always trained on all rows.

    Transition-level scoring reuses the historical OpenMS / PyProphet-style
    candidate filters: only transitions from non-decoy precursors whose
    peak-group rank and PEP pass the configured thresholds and whose
    transition-level isotope-overlap / log-SN filters pass are rescored.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathPercolatorScoring :
    public DefaultParamHandler
  {
  public:
    /// Supported OpenSWATH scoring levels
    enum class Level
    {
      MS1,
      MS2,
      MS1MS2,
      TRANSITION,
      SIZE_OF_LEVEL
    };

    /// String names of the supported scoring levels
    static const std::array<std::string, static_cast<Size>(Level::SIZE_OF_LEVEL)> names_of_level;

    /// Small execution summary returned by score()
    struct OPENMS_DLLAPI ScoreSummary
    {
      Size total_rows = 0;
      Size target_rows = 0;
      Size decoy_rows = 0;
      Size feature_count = 0;
    };

    OpenSwathPercolatorScoring();
    ~OpenSwathPercolatorScoring() override;

    /**
      @brief Score an OpenSWATH OSW / OSWPQ file in place or into a copy.

      Input order is canonicalized by run filename, peptide or compound identity, charge and decoy label,
      then retention time and retained feature values. Generated IDs are used only to preserve the
      one-run/one-precursor grouping and to write results. Transition rows remain separate CV groups.
      Retained rows require run filenames and modified peptide sequences or compound identity.
      Native SQL compounds use name, formula, SMILES and adducts. Native Parquet molecules use the
      source library traml_id and precursor m/z when modified_sequence is empty. Source library
      identities must remain stable between runs. Incomplete active-score rows are discarded first.
      Distinct assay aliases with identical scientific sort keys keep their original partitions,
      but alias-specific CV assignments are not guaranteed invariant.

      @param[in] input_path   Input `.osw` SQLite file or `.oswpq` directory/archive
      @param[in] level        OpenSWATH scoring level to rescore
      @param[in] output_path  Optional output path. If empty, score in place.
                              Non-empty outputs preserve the input container kind
                              (SQLite -> SQLite, directory -> directory, archive -> archive).

      @return Small execution summary for the rescored rows.

      @throws Exception::FileNotFound if the input file does not exist
      @throws Exception::Precondition if required prerequisite scores are missing
      @throws Exception::MissingInformation if a retained row has no run filename, scientific precursor identity, or valid transition decoy label
      @throws Exception::InvalidValue if no usable feature columns are available or a retention time is not finite
    */
    ScoreSummary score(const std::string& input_path,
                       Level level,
                       const std::string& output_path = "");

    /**
      @brief The rows Percolator is trained on for peak-group (ms1, ms2, ms1ms2) rescoring.

      Rows sharing a value of RescoreInput::cv_group_keys form a group (one precursor in one run).
      The training set is the best row of every group by @p main_score (ties keep the first row) plus
      all rows of groups drawn at random, in the order of a shuffle with @p seed, as long as their other
      rows fit into @p max_others (the first group that does not fit ends the drawing).

      @param[in] input       Peak groups; cv_group_keys must hold one key per row
      @param[in] main_score  Feature name; prefix '-' if lower values are better
      @param[in] max_others  Rows besides the best of every group (0: all rows)
      @param[in] seed        Seed of the shuffle that draws the groups

      @return Ascending row indices. All rows if @p max_others is 0, if cv_group_keys does not hold one
              key per row, if @p main_score is not among the features, or if all rows would fit.
    */
    static std::vector<Size> trainingRows(const RescoreInput& input,
                                          const std::string& main_score,
                                          Size max_others,
                                          UInt64 seed);

  protected:
    void updateMembers_() override;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
} // namespace OpenMS
