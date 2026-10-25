// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace OpenDIA
{
  using RecordId = std::uint64_t;

  /// Unstable v0 context. IDs are caller-assigned, never inferred or renumbered.
  struct RecordContext
  {
    unsigned schema_version = 0;
    RecordId run_id = 0;
    RecordId lane_id = 0;
    RecordId calibration_snapshot_id = 0;
    bool operator==(const RecordContext&) const = default;
  };

  /// Caller-specified maximum record and flattened value counts (zero is valid).
  struct RecordLimits
  {
    OpenMS::Size records = 0;
    OpenMS::Size values = 0;
  };

  enum class MobilityUnit { MILLISECOND, INVERSE_REDUCED };
  enum class EdgePolicy { OPEN, CLOSED, LOWER_CLOSED, UPPER_CLOSED };
  /// Absolute physical bounds. Equal endpoints are allowed, including empty open intervals.
  struct QueryInterval
  {
    double lower = 0;
    double upper = 0;
    EdgePolicy edges = EdgePolicy::OPEN;
    bool operator==(const QueryInterval&) const = default;
  };
  struct MobilityInterval
  {
    QueryInterval interval;
    MobilityUnit unit = MobilityUnit::MILLISECOND;
    bool operator==(const MobilityInterval&) const = default;
  };
  enum class ToleranceUnit { TH, PPM };
  struct QueryTolerance
  {
    double full_width = 0;
    ToleranceUnit unit = ToleranceUnit::TH;
    bool operator==(const QueryTolerance&) const = default;
  };
  enum class QueryAggregation { RAW_POINTS, SUM };
  enum class IsotopeMode { NONE, SELECTED };
  /// SELECTED records a nonnegative isotope index and positive charge magnitude.
  /// NONE requires index and charge to be zero. Bounds are never derived here.
  struct IsotopePolicy
  {
    IsotopeMode mode = IsotopeMode::NONE;
    unsigned index = 0;
    unsigned charge = 0;
    bool operator==(const IsotopePolicy&) const = default;
  };
  /// Unstable v0 tophat-only request. Producers supply canonical physical bounds.
  struct ExtractionQuerySpec
  {
    RecordContext context;
    /// Absent means all acquisitions, present empty means none. IDs are opaque pairs.
    std::optional<std::vector<std::pair<RecordId, RecordId>>> acquisitions;
    unsigned ms_level = 1;
    std::optional<QueryInterval> isolation_mz;
    std::optional<QueryInterval> rt_seconds;
    QueryInterval mz;
    std::optional<MobilityInterval> mobility;
    std::optional<double> faims_voltage;
    QueryTolerance tolerance;
    QueryAggregation aggregation = QueryAggregation::RAW_POINTS;
    IsotopePolicy isotope;
    bool operator==(const ExtractionQuerySpec&) const = default;
  };
  /** @brief Checked, tool-private v0 value key. No extraction or cache integration. */
  class ExtractionQuery
  {
  public:
    /**
      @brief Validate a specification and sort/deduplicate only its acquisition set.
      @param[in] specification Physical request and immutable caller snapshot identity.
      @throws OpenMS::Exception::InvalidValue for malformed input.
    */
    explicit ExtractionQuery(ExtractionQuerySpec specification);
    ExtractionQuerySpec specification() const;
    bool operator==(const ExtractionQuery&) const = default;
  private:
    ExtractionQuerySpec specification_;
  };

  struct Mobility
  {
    double value = 0;
    MobilityUnit unit = MobilityUnit::MILLISECOND;
  };

  /// Exact measured coordinates, not lookup bins. Isolation bounds are absolute m/z.
  struct Coordinates
  {
    double rt_seconds = 0;
    double mz = 0;
    RecordId frame_id = 0;
    RecordId scan_id = 0;
    unsigned ms_level = 1;
    double isolation_lower_mz = 0;
    double isolation_upper_mz = 0;
    std::optional<Mobility> mobility;
    std::optional<double> faims_voltage;
    /// Exact source acquisition ID, independent of the caller's numeric namespace.
    std::string native_id;
  };

  struct RawPoint
  {
    Coordinates coordinates;
    double intensity = 0;
  };

  enum class TargetDecoy { TARGET, DECOY };
  struct Provenance
  {
    std::string library;
    std::optional<std::string> prediction;
  };
  struct Hypothesis
  {
    RecordId family_id = 0;
    RecordId precursor_id = 0;
    double precursor_mz = 0;
    TargetDecoy target_decoy = TargetDecoy::TARGET;
    Provenance provenance;
    std::optional<double> predicted_rt_seconds;
    std::optional<Mobility> predicted_mobility;
  };
  struct Fragment
  {
    RecordId fragment_id = 0;
    double mz = 0;
    std::optional<double> predicted_intensity;
  };
  struct EvidencePoint
  {
    RecordId family_id = 0;
    RecordId precursor_id = 0;
    RecordId fragment_id = 0;
    Coordinates coordinates;
    std::optional<double> value;
  };

  /** @brief Unstable v0 immutable owning raw lease. All access returns values.
   * Copies retain storage. Moved-from access throws OpenMS::Exception::InvalidValue.
   */
  class RawTile
  {
  public:
    /**
      @brief Validate and own a bounded raw tile.
      @param[in] context Schema and caller identity namespace.
      @param[in] limits Maximum rows and values (both bound points).
      @param[in] points Exact raw measurements, moved into immutable storage.
      @throws OpenMS::Exception::InvalidValue for malformed input.
    */
    RawTile(RecordContext context, RecordLimits limits, std::vector<RawPoint> points);
    RecordContext context() const;
    OpenMS::Size size() const;
    RawPoint point(OpenMS::Size index) const;
  private:
    struct Storage;
    std::shared_ptr<const Storage> storage_;
    void check_() const;
  };

  /** @brief Unstable v0 immutable owning hypothesis lease with CSR fragment offsets (size + 1). */
  class HypothesisPage
  {
  public:
    /**
      @brief Validate and own a bounded sparse hypothesis page.
      @param[in] context Schema and caller identity namespace.
      @param[in] limits Maximum hypotheses and flattened fragments.
      @param[in] hypotheses Assay identities and predictions.
      @param[in] fragment_offsets Monotone CSR offsets, starting at zero and ending at fragment count.
      @param[in] fragments Flattened sparse fragment records.
      @throws OpenMS::Exception::InvalidValue for malformed input.
    */
    HypothesisPage(RecordContext context, RecordLimits limits, std::vector<Hypothesis> hypotheses,
                   std::vector<OpenMS::Size> fragment_offsets, std::vector<Fragment> fragments);
    RecordContext context() const;
    OpenMS::Size size() const;
    Hypothesis hypothesis(OpenMS::Size index) const;
    OpenMS::Size fragmentCount(OpenMS::Size hypothesis_index) const;
    Fragment fragment(OpenMS::Size hypothesis_index, OpenMS::Size fragment_index) const;
  private:
    struct Storage;
    std::shared_ptr<const Storage> storage_;
    void check_() const;
  };

  /**
    @brief Unstable v0 sparse evidence owning lease.

    An absent row is not a measured zero or an explicit missing row.
  */
  class EvidenceBlock
  {
  public:
    /**
      @brief Validate and own bounded sparse support.
      @param[in] context Schema and caller identity namespace.
      @param[in] limits Maximum rows and values (both bound points).
      @param[in] points Explicit observations, including missing values.
      @throws OpenMS::Exception::InvalidValue for malformed input.
    */
    EvidenceBlock(RecordContext context, RecordLimits limits, std::vector<EvidencePoint> points);
    RecordContext context() const;
    OpenMS::Size size() const;
    EvidencePoint point(OpenMS::Size index) const;
  private:
    struct Storage;
    std::shared_ptr<const Storage> storage_;
    void check_() const;
  };
}
