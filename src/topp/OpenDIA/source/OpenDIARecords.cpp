// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARecords.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <cmath>
#include <limits>
#include <utility>

namespace OpenDIA
{
  namespace
  {
    void require(bool condition, const char* message)
    {
      if (!condition)
      {
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, "OpenDIA v0");
      }
    }
    void finite(double value)
    {
      require(std::isfinite(value), "Nonfinite record value");
    }
    void optionalFinite(const std::optional<double>& value)
    {
      if (value)
      {
        finite(*value);
      }
    }
    void mobility(const std::optional<Mobility>& value)
    {
      if (!value)
      {
        return;
      }
      require(value->unit == MobilityUnit::MILLISECOND || value->unit == MobilityUnit::INVERSE_REDUCED,
              "Unknown mobility unit");
      finite(value->value);
    }
    void coordinates(const Coordinates& value)
    {
      finite(value.rt_seconds);
      finite(value.mz);
      finite(value.isolation_lower_mz);
      finite(value.isolation_upper_mz);
      require(value.ms_level > 0, "MS level must be nonzero");
      require(value.isolation_lower_mz <= value.isolation_upper_mz, "Reversed isolation geometry");
      mobility(value.mobility);
      optionalFinite(value.faims_voltage);
    }
    void context(const RecordContext& value)
    {
      require(value.schema_version == 0, "Unsupported record schema");
    }
    void bounds(OpenMS::Size index, OpenMS::Size size)
    {
      require(index < size, "Record index out of bounds");
    }
  }

  struct RawTile::Storage
  {
    RecordContext context;
    std::vector<RawPoint> points;
  };
  RawTile::RawTile(RecordContext ctx, RecordLimits limits, std::vector<RawPoint> points)
  {
    OpenDIA::context(ctx);
    require(points.size() <= limits.records && points.size() <= limits.values, "Raw tile cap exceeded");
    for (const auto& point : points)
    {
      coordinates(point.coordinates);
      finite(point.intensity);
    }
    storage_ = std::make_shared<const Storage>(Storage{ctx, std::move(points)});
  }
  void RawTile::check_() const { require(bool(storage_), "Invalid raw tile lease"); }
  RecordContext RawTile::context() const { check_(); return storage_->context; }
  OpenMS::Size RawTile::size() const { check_(); return storage_->points.size(); }
  RawPoint RawTile::point(OpenMS::Size index) const
  {
    bounds(index, size());
    return storage_->points[index];
  }

  struct HypothesisPage::Storage
  {
    RecordContext context;
    std::vector<Hypothesis> hypotheses;
    std::vector<OpenMS::Size> offsets;
    std::vector<Fragment> fragments;
  };
  HypothesisPage::HypothesisPage(RecordContext ctx, RecordLimits limits, std::vector<Hypothesis> hypotheses,
                                 std::vector<OpenMS::Size> offsets, std::vector<Fragment> fragments)
  {
    OpenDIA::context(ctx);
    require(hypotheses.size() <= limits.records && fragments.size() <= limits.values, "Hypothesis page cap exceeded");
    // Check before adding one, and reject malformed CSR without indexing its endpoints.
    require(hypotheses.size() < std::numeric_limits<OpenMS::Size>::max(), "Offset shape overflow");
    require(offsets.size() == hypotheses.size() + 1, "Invalid fragment offset shape");
    require(offsets.front() == 0 && offsets.back() == fragments.size(), "Invalid fragment offset endpoints");
    for (OpenMS::Size i = 1; i < offsets.size(); ++i)
    {
      require(offsets[i - 1] <= offsets[i] && offsets[i] <= fragments.size(), "Invalid fragment offsets");
    }
    for (const auto& hypothesis : hypotheses)
    {
      finite(hypothesis.precursor_mz);
      require(hypothesis.target_decoy == TargetDecoy::TARGET || hypothesis.target_decoy == TargetDecoy::DECOY,
              "Unknown target/decoy identity");
      optionalFinite(hypothesis.predicted_rt_seconds);
      mobility(hypothesis.predicted_mobility);
    }
    for (const auto& fragment : fragments)
    {
      finite(fragment.mz);
      optionalFinite(fragment.predicted_intensity);
    }
    storage_ = std::make_shared<const Storage>(Storage{ctx, std::move(hypotheses), std::move(offsets), std::move(fragments)});
  }
  void HypothesisPage::check_() const { require(bool(storage_), "Invalid hypothesis page lease"); }
  RecordContext HypothesisPage::context() const { check_(); return storage_->context; }
  OpenMS::Size HypothesisPage::size() const { check_(); return storage_->hypotheses.size(); }
  Hypothesis HypothesisPage::hypothesis(OpenMS::Size index) const
  {
    bounds(index, size());
    return storage_->hypotheses[index];
  }
  OpenMS::Size HypothesisPage::fragmentCount(OpenMS::Size index) const
  {
    bounds(index, size());
    return storage_->offsets[index + 1] - storage_->offsets[index];
  }
  Fragment HypothesisPage::fragment(OpenMS::Size index, OpenMS::Size fragment_index) const
  {
    bounds(fragment_index, fragmentCount(index));
    // The validated count bounds this addition by fragments.size(), preventing overflow.
    return storage_->fragments[storage_->offsets[index] + fragment_index];
  }

  struct EvidenceBlock::Storage
  {
    RecordContext context;
    std::vector<EvidencePoint> points;
  };
  EvidenceBlock::EvidenceBlock(RecordContext ctx, RecordLimits limits, std::vector<EvidencePoint> points)
  {
    OpenDIA::context(ctx);
    require(points.size() <= limits.records && points.size() <= limits.values, "Evidence block cap exceeded");
    for (const auto& point : points)
    {
      coordinates(point.coordinates);
      optionalFinite(point.value);
    }
    storage_ = std::make_shared<const Storage>(Storage{ctx, std::move(points)});
  }
  void EvidenceBlock::check_() const { require(bool(storage_), "Invalid evidence block lease"); }
  RecordContext EvidenceBlock::context() const { check_(); return storage_->context; }
  OpenMS::Size EvidenceBlock::size() const { check_(); return storage_->points.size(); }
  EvidencePoint EvidenceBlock::point(OpenMS::Size index) const
  {
    bounds(index, size());
    return storage_->points[index];
  }
}
