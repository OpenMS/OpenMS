// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

namespace OpenMS
{
  class ConsensusMap;
  class ExperimentalDesign;

/**
  @brief Whole-collection preflight for the QPX Parquet export (@c -out_qpx)

  A QPX collection is three files written in a fixed order -- @c quantms.feature.parquet, then
  @c quantms.psm.parquet, then @c quantms.pg.parquet. Each view refuses input it cannot
  represent, but it does so before opening <i>its own</i> file, so a refusal in the second or
  third view leaves the earlier files on disk: a partial collection that looks like a
  successful export.

  The exporters cannot avoid that themselves. Both streaming paths open the output and emit a
  valid 0-row file even for empty input, so by the time a view can refuse, its predecessors
  have already written. Only a check that runs before the first write can leave the output
  directory untouched.

  requireExportable() therefore runs every deterministic refusal of all three views up front.
  It does not replace the in-exporter checks -- those remain the contract for callers using an
  exporter directly -- so a condition is validated twice on the TOPP path. That is deliberate:
  the preflight must not become the only thing standing between a library user and a malformed
  table.

  @note Only conditions that are decidable from the input are covered. I/O failures, and
        anything that can only be discovered while building a row, still leave a partial
        collection; that is a narrower hazard and is tracked separately.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI QPXCollectionExport
{
public:
  /**
    @brief Refuse a collection that cannot be written in full, before any file is created

    Performs, on the same input the three views will see:
      - the feature view's divergent-identification and merge-index refusals,
      - the psm view's merge-index refusal over assigned <i>and</i> unassigned identifications,
      - the channel-label check shared by the feature and pg views,
      - the pg view's quantification-unit checks (inconsistent, ragged, ambiguous label,
        duplicated sample) and the per-group sample-abundance length check.

    @param[in] cmap The ConsensusMap about to be exported
    @param[in] design The experimental design that was used to quantify @p cmap
    @return true if all three views can be written; false for the conditions the views report
            by returning false/nullptr, having logged the reason
    @throw Exception::IllegalArgument if a feature carries divergent peptide annotations
    @throw Exception::MissingInformation if an identification of a merged run carries no usable
           @c id_merge_index
    @throw Exception::InvalidValue if two identification runs share an identifier
  */
  static bool requireExportable(const ConsensusMap& cmap, const ExperimentalDesign& design);
};

} // namespace OpenMS
