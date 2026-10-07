// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <functional>

namespace OpenMS
{
/**
  @brief Native typed Parquet persistence and bounded sequential access to IdentificationData.

  One manifest contains configuration and row ranges in shared query/match files.
  Compatible runs share row groups; a run_uuid (or inference_identifier) column records which run (or inference
  result) owns each row. Scores stay in match rows. Growing inference values use typed tables.
  Persistent UUIDs and record IDs survive filtering. No revision or freshness policy
  is implicit. Publication rejects an existing destination unless Options::replace_existing
  allows replacing a native bundle; owning load is transactional.
  The experimental schema has no compatibility reader for discarded prototypes.
  @ingroup FileIO
*/
class OPENMS_DLLAPI IdentificationDataFile
{
public:
  struct OPENMS_DLLAPI Options
  {
    Size batch_rows = 16384;
    Size row_group_rows = 65536;
    Size batch_bytes = 8 * 1024 * 1024;
    Size row_group_bytes = 64 * 1024 * 1024;
    Size max_record_bytes = 64 * 1024 * 1024;
    /// Positive CPU worker limit per operation; 1 is serial. Does not change Arrow's global pool.
    Size threads = 1;
    /// Writing replaces an existing native identification bundle at the destination (moved aside,
    /// then deleted after publication). Any other existing file or directory is still rejected.
    bool replace_existing = false;
  };
  struct OPENMS_DLLAPI Projection
  {
    bool molecule = true;
    bool evidence = true;
    bool annotations = true;
    bool metadata = true;
    bool all_scores = true;
    std::vector<UInt32> score_ids;
  };
  struct OPENMS_DLLAPI ScanOptions
  {
    Options buffering;
    Projection projection;
    /// Empty selects every run; entries are display identifiers or UUIDs.
    std::vector<std::string> runs;
    /// Optional strict uniqueness validation uses sets proportional to the scanned run.
    bool validate_unique_ids = false;
  };
  struct OPENMS_DLLAPI RunDescriptor
  {
    std::string identifier;
    std::string uuid;
    IdentificationData::MoleculeKind molecule_kind;
    std::vector<IdentificationData::ScoreDefinition> scores;
    /// Names of the score columns in matches.parquet, parallel to @p scores (e.g. score_pep, score_q_value).
    std::vector<std::string> score_columns;
    std::vector<IdentificationData::SourceFile> sources;
    std::vector<IdentificationData::Database> databases;
    std::optional<UInt32> primary_score;
    UInt64 query_count = 0;
    UInt64 match_count = 0;
    UInt64 next_query_id = 1;
    UInt64 next_match_id = 1;
  };
  struct OPENMS_DLLAPI QueryRecord
  {
    UInt64 query_id = 0;
    UInt32 source_id = 0;
    IdentificationData::Observation data;
    std::optional<UInt64> selected_match_id;
  };
  struct OPENMS_DLLAPI MatchRecord
  {
    UInt64 match_id = 0;
    UInt64 query_id = 0;
    IdentificationData::MatchData data;
    /// Dense run-score positions. Unprojected positions are null; consult Projection.
    std::vector<std::optional<double>> scores;
  };
  struct OPENMS_DLLAPI ScanStatistics
  {
    UInt64 queries = 0;
    UInt64 matches = 0;
    UInt64 descriptor_bytes = 0; ///< Resident manifest text size, apart from decoded descriptor overhead.
  };
  using QueryCallback = std::function<void(const std::string&, const std::vector<QueryRecord>&)>;
  using MatchCallback = std::function<void(const std::string&, const std::vector<MatchRecord>&)>;
  using MatchPredicate = std::function<bool(const std::string&, const MatchRecord&)>;

  /// Write through a temporary sibling directory, closing every table before publication.
  static void store(const std::string& path, const IdentificationData& data);
  static void store(const std::string& path, const IdentificationData& data, const Options& options);
  /// Validate into a temporary collection; destination is unchanged on any failure.
  static void load(const std::string& path, IdentificationData& data);
  static void load(const std::string& path, IdentificationData& data, const Options& options);
  /// Load only the selected run. Dataset-level inference is not implicitly loaded.
  static IdentificationData::Run loadRun(const std::string& path, const std::string& run);
  static IdentificationData::Run loadRun(const std::string& path, const std::string& run, const Options& options);
  /// Recognize the manifest format; malformed manifests report an error.
  static bool isNativeFile(const std::string& path);
  /// Inspect configuration without opening any Parquet table.
  static std::vector<RunDescriptor> inspect(const std::string& path);
  /**
    @brief Scan projected typed values; callbacks own no Arrow objects.

    Query and match callbacks are independent batch streams. Candidates for a query
    may span batches. IDs associate the streams. Callback exceptions propagate; a late
    input error can occur after prior callbacks. Buffers are valid during the callback.
  */
  static ScanStatistics scan(const std::string& path, const ScanOptions& options, const QueryCallback& queries, const MatchCallback& matches);
  /// Stream a self-contained subset, retaining IDs/counters. No candidate-sized lookup map.
  static void filter(const std::string& input,
                     const std::string& output,
                     const MatchPredicate& keep,
                     IdentificationData::InferencePolicy inference_policy,
                     bool keep_empty_queries = false);
  static void filter(const std::string& input,
                     const std::string& output,
                     const MatchPredicate& keep,
                     IdentificationData::InferencePolicy inference_policy,
                     bool keep_empty_queries,
                     const Options& options);
};
} // namespace OpenMS
