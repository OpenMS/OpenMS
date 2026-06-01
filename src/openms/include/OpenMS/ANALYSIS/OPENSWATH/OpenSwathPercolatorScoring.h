// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/config.h>

#include <array>
#include <memory>

namespace OpenMS
{
  /**
    @brief In-process Percolator scoring for OpenSWATH OSW and OSWPQ results.

    The scorer reads OpenSWATH feature tables, converts them into a
    domain-agnostic Percolator::RescoreInput, performs target/decoy rescoring,
    and writes the resulting score, q-value, PEP, and peak-group rank back to
    the source format.

    Supported levels:
    - MS1: writes SCORE_MS1 (SQLite) or score_ms1_* columns (.oswpq)
    - MS2: writes SCORE_MS2 (SQLite) or score_ms2_* columns (.oswpq)
    - MS1MS2: uses combined MS2 + MS1 features and writes SCORE_MS2 /
      score_ms2_* (matching the historical PyProphet convention)
    - Transition: requires prior MS2/MS1MS2 scoring and writes
      SCORE_TRANSITION / score_transition_*

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

      @param[in] input_path   Input `.osw` SQLite file or `.oswpq` directory/archive
      @param[in] level        OpenSWATH scoring level to rescore
      @param[in] output_path  Optional output path. If empty, score in place.
                              Non-empty outputs preserve the input container kind
                              (SQLite -> SQLite, directory -> directory, archive -> archive).

      @return Small execution summary for the rescored rows.

      @throws Exception::FileNotFound if the input file does not exist
      @throws Exception::Precondition if required prerequisite scores are missing
      @throws Exception::InvalidValue if no usable feature columns are available
    */
    ScoreSummary score(const String& input_path,
                       Level level,
                       const String& output_path = "");

  protected:
    void updateMembers_() override;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
} // namespace OpenMS
