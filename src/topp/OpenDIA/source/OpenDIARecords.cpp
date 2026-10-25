// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARecords.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
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
    void interval(const QueryInterval& value)
    {
      finite(value.lower);
      finite(value.upper);
      require(value.lower <= value.upper, "Reversed query interval");
      require(value.edges == EdgePolicy::OPEN || value.edges == EdgePolicy::CLOSED ||
              value.edges == EdgePolicy::LOWER_CLOSED || value.edges == EdgePolicy::UPPER_CLOSED,
              "Unknown interval edge policy");
    }
  }

  ExtractionQuery::ExtractionQuery(ExtractionQuerySpec specification) : specification_(std::move(specification))
  {
    auto& value = specification_;
    OpenDIA::context(value.context);
    require(value.ms_level > 0, "MS level must be nonzero");
    interval(value.mz);
    if (value.isolation_mz) { interval(*value.isolation_mz); }
    if (value.rt_seconds) { interval(*value.rt_seconds); }
    if (value.mobility)
    {
      interval(value.mobility->interval);
      mobility(Mobility{value.mobility->interval.lower, value.mobility->unit});
    }
    optionalFinite(value.faims_voltage);
    finite(value.tolerance.full_width);
    require(value.tolerance.full_width >= 0, "Negative full-width tolerance");
    require(value.tolerance.unit == ToleranceUnit::TH || value.tolerance.unit == ToleranceUnit::PPM,
            "Unknown tolerance unit");
    require(value.aggregation == QueryAggregation::RAW_POINTS || value.aggregation == QueryAggregation::SUM,
            "Unknown aggregation");
    require((value.isotope.mode == IsotopeMode::NONE && value.isotope.index == 0 && value.isotope.charge == 0) ||
            (value.isotope.mode == IsotopeMode::SELECTED && value.isotope.charge > 0), "Invalid isotope policy");
    // Acquisition order has no scientific meaning. No physical value is rounded or converted.
    if (value.acquisitions)
    {
      auto& selection = *value.acquisitions;
      std::sort(selection.begin(), selection.end());
      selection.erase(std::unique(selection.begin(), selection.end()), selection.end());
    }
  }
  ExtractionQuerySpec ExtractionQuery::specification() const { return specification_; }

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

  namespace
  {
    // Subtract validated counts instead of accumulating potentially overflowing sums.
    void consume(OpenMS::Size& remaining, OpenMS::Size count)
    {
      require(count <= remaining, "Flattened value cap exceeded");
      remaining -= count;
    }
  }
  struct FamilyBatch::Storage
  {
    RecordContext context;
    std::vector<CandidateFamily> families;
  };
  FamilyBatch::FamilyBatch(RecordContext ctx, RecordLimits limits, std::vector<CandidateFamily> families)
  {
    OpenDIA::context(ctx);
    require(families.size() <= limits.records, "Family record cap exceeded");
    std::set<RecordId> family_ids;
    auto remaining = limits.values;
    for (const auto& family : families)
    {
      require(family_ids.insert(family.family_id).second, "Duplicate family ID");
      require(family.region.specification().context == ctx, "Family region context mismatch");
      require(family.completeness == CandidateCompleteness::UNKNOWN ||
              family.completeness == CandidateCompleteness::PARTIAL ||
              family.completeness == CandidateCompleteness::COMPLETE, "Unknown candidate completeness");
      consume(remaining, family.members.size());
      consume(remaining, family.alternatives.size());
      consume(remaining, family.fragments.size());
      std::set<RecordId> members;
      for (const auto& member : family.members)
      {
        require(members.insert(member.precursor_id).second, "Duplicate family member ID");
        require(member.target_decoy == TargetDecoy::TARGET || member.target_decoy == TargetDecoy::DECOY,
                "Unknown target/decoy identity");
        require(!member.charge || *member.charge != 0, "Present charge must be nonzero");
        require(!member.adduct || !member.adduct->empty(), "Present adduct must be nonempty");
      }
      std::set<RecordId> alternatives;
      for (const auto& alternative : family.alternatives)
      {
        require(alternatives.insert(alternative.alternative_id).second, "Duplicate alternative ID");
        consume(remaining, alternative.precursor_ids.size());
        std::set<RecordId> references;
        for (const auto id : alternative.precursor_ids)
        {
          require(members.contains(id), "Dangling alternative member");
          require(references.insert(id).second, "Duplicate alternative member reference");
        }
      }
      std::set<std::pair<RecordId, RecordId>> incidences;
      for (const auto& fragment : family.fragments)
      {
        require(members.contains(fragment.precursor_id), "Dangling fragment member");
        require(fragment.role == FragmentRole::SHARED || fragment.role == FragmentRole::DISCRIMINATING,
                "Unknown fragment role");
        require(incidences.emplace(fragment.precursor_id, fragment.fragment_id).second, "Duplicate fragment incidence");
      }
    }
    storage_ = std::make_shared<const Storage>(Storage{ctx, std::move(families)});
  }
  void FamilyBatch::check_() const { require(bool(storage_), "Invalid family batch lease"); }
  RecordContext FamilyBatch::context() const { check_(); return storage_->context; }
  OpenMS::Size FamilyBatch::size() const { check_(); return storage_->families.size(); }
  CandidateFamily FamilyBatch::family(OpenMS::Size index) const
  {
    bounds(index, size());
    return storage_->families[index];
  }

  struct ScoreBatch::Storage
  {
    FamilyBatch families;
    EvidenceBlock evidence;
    std::vector<ScoreRecord> scores;
  };
  ScoreBatch::ScoreBatch(RecordLimits limits, FamilyBatch families, EvidenceBlock evidence, std::vector<ScoreRecord> scores)
  {
    require(families.context() == evidence.context(), "Score input context mismatch");
    require(scores.size() <= limits.records, "Score record cap exceeded");
    // Index once. Validation stays sparse and does not repeatedly scan all families.
    struct References
    {
      std::set<RecordId> members;
      std::set<RecordId> alternatives;
    };
    std::map<RecordId, References> references;
    for (OpenMS::Size i = 0; i < families.size(); ++i)
    {
      const auto family = families.family(i);
      auto& refs = references[family.family_id];
      for (const auto& member : family.members) { refs.members.insert(member.precursor_id); }
      for (const auto& alternative : family.alternatives) { refs.alternatives.insert(alternative.alternative_id); }
    }
    std::set<std::pair<RecordId, RecordId>> rows;
    auto remaining = limits.values;
    for (const auto& score : scores)
    {
      const auto family = references.find(score.family_id);
      require(family != references.end(), "Dangling score family");
      require(family->second.alternatives.contains(score.alternative_id), "Dangling score alternative");
      require(rows.emplace(score.family_id, score.alternative_id).second, "Duplicate score row");
      require(score.ambiguity == ScoreAmbiguity::UNRESOLVED || score.ambiguity == ScoreAmbiguity::RESOLVED,
              "Unknown score ambiguity");
      require(!score.scoring_method.empty(), "Empty scoring method");
      optionalFinite(score.identification);
      optionalFinite(score.localization);
      if (score.quantity)
      {
        finite(score.quantity->value);
        require(!score.quantity->unit.empty(), "Empty quantity unit");
      }
      consume(remaining, score.evidence_rows.size());
      for (const auto row : score.evidence_rows)
      {
        const auto point = evidence.point(row);
        require(point.family_id == score.family_id && family->second.members.contains(point.precursor_id),
                "Score evidence family/member mismatch");
      }
    }
    storage_ = std::make_shared<const Storage>(Storage{std::move(families), std::move(evidence), std::move(scores)});
  }
  void ScoreBatch::check_() const { require(bool(storage_), "Invalid score batch lease"); }
  RecordContext ScoreBatch::context() const { check_(); return storage_->families.context(); }
  OpenMS::Size ScoreBatch::size() const { check_(); return storage_->scores.size(); }
  ScoreRecord ScoreBatch::score(OpenMS::Size index) const
  {
    bounds(index, size());
    return storage_->scores[index];
  }
  FamilyBatch ScoreBatch::families() const { check_(); return storage_->families; }
  EvidenceBlock ScoreBatch::evidence() const { check_(); return storage_->evidence; }
}
