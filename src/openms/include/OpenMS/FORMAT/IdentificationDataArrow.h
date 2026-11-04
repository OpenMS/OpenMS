// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <map>
#include <memory>

namespace arrow
{
class RecordBatchReader;
class Table;
} // namespace arrow

namespace OpenMS
{
/**
  @brief Bulk access to the matches of an IdentificationData as an Arrow table, and edits as patches

  matchTable() returns one row per match with the columns of matches.parquet in a native bundle
  (IdentificationDataFile): match_id, query_id, the molecule, evidence, annotation and metadata columns, one
  score column per score definition (score_q_value, ...), and run_uuid last (a dictionary-encoded string).
  Runs follow each other, and the matches of a run follow their queries. The schema metadata describes the
  table: "openms:score_definitions" and "openms:metadata_descriptors" (JSON, as in the manifest of a bundle;
  the metadata column refers to the descriptors by index), "openms:primary_score" (index or null) and
  "openms:revisions" (JSON object: run UUID -> Run::getRevision()).

  applyPatch() edits matches from a table keyed by run_uuid and match_id: score columns (named as in
  matchTable(); null removes a supplementary value) and target_decoy. Every key and value is checked first;
  then all edits are applied, or none (the patched runs are edited as copies and swapped in). A patch can add
  score definitions (PatchOptions::add_scores) to fill new columns, e.g. after rescoring, and can require the
  revisions of the runs it was computed from (PatchOptions::expected_revisions; not checked by default, since
  tables rebuilt elsewhere, e.g. by pandas or polars, lose the schema metadata of matchTable()).

  Applying a patch replaces the data of the patched runs, like other structural edits: references to their queries
  and matches are invalidated, and score views bound before are rejected (see IdentificationData::ScoreView). The
  runs themselves stay in place. Queries, the molecules of matches and their other fields are not patched.
  @ingroup FileIO
*/
class OPENMS_DLLAPI IdentificationDataArrow
{
public:
  struct OPENMS_DLLAPI PatchOptions
  {
    /// Score definitions added to every run of the score schema (runs with score definitions or matches,
    /// not sequence catalogs) before the patch is applied; definitions that already exist are reused.
    std::vector<IdentificationData::ScoreDefinition> add_scores;
    /// Revision (Run::getRevision()) that a run must still have, by run UUID; runs not listed are not checked.
    std::map<std::string, UInt64> expected_revisions;
  };
  struct OPENMS_DLLAPI PatchResult
  {
    Size rows = 0;                     ///< Matches edited (one per patch row)
    std::vector<std::string> columns;  ///< Patched columns, in patch order
    std::vector<std::string> added;    ///< Score columns added by PatchOptions::add_scores
  };

  /// One row per match, with the columns of matches.parquet (see the class description).
  static std::shared_ptr<arrow::Table> matchTable(const IdentificationData& data);
  /// Only the columns of @p projection (the key columns match_id, query_id and run_uuid always).
  static std::shared_ptr<arrow::Table> matchTable(const IdentificationData& data, const IdentificationDataFile::Projection& projection);

  /// @throw Exception::InvalidParameter for an invalid patch (unknown column, run, match or a repeated key; a value
  /// that is not finite; a missing primary score; a target_decoy out of range; a changed revision). Nothing changes then.
  static PatchResult applyPatch(IdentificationData& data, const arrow::Table& patch);
  static PatchResult applyPatch(IdentificationData& data, const arrow::Table& patch, const PatchOptions& options);
  /// The same for a stream of record batches (e.g. from the Arrow C stream interface); the batches are read first.
  static PatchResult applyPatch(IdentificationData& data, arrow::RecordBatchReader& patch, const PatchOptions& options);
};
} // namespace OpenMS
