/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/data/array_index.h"

#include <algorithm>
#include <ranges>

#include "mzpeak/exception.h"
#include "mzpeak/schema/buffer_format.h"
#include "mzpeak/schema/entity_type.h"

namespace MzPeak::Data {

/******************************************************************************/
// How to compare entries.  Must match EntryChunkFn below.
//
// TODO: Unify these two types.
struct EntryCmpFn {
  bool operator()(const ArrayIndex::Entry& a, const ArrayIndex::Entry& b) const
  {
    // Lexicographic, NOT the &&-chained form.  `a.x < b.x && a.y < b.y && ...`
    // is not a strict weak ordering -- two entries that differ in opposite
    // directions compare false BOTH ways, so the sort treats them as
    // equivalent and equivalence stops being transitive.  std::ranges::sort
    // requires a strict weak order, so that form is undefined behaviour, and
    // this sort decides which entry a dimension sees first (i.e. which one
    // counts as primary when coalescing).  Upstream landed the same fix in
    // 4637f85.
    if (a.array_name != b.array_name) return a.array_name < b.array_name;
    if (a.array_type != b.array_type) return a.array_type < b.array_type;
    if (a.data_type != b.data_type) return a.data_type < b.data_type;
    // Primary entries first.
    return a.buffer_priority > b.buffer_priority;
  }
};

// How to chunk entries.  Must match EntryCmpFn above.
//
// NOTE: buffer_priority is deliberately NOT part of this key, although the sort
// comparator orders on it.  Upstream groups on it so that a primary and a
// non-primary entry for the same array become separate dimensions and the
// primary one can be preferred -- correct when the two are REDUNDANT copies of
// the same data.  But they can also be COMPLEMENTARY: has_uv.mzpeak stores one
// chromatogram's intensities in detector counts and another's in absorbance, in
// two columns, each null where the other has values.  Splitting those into two
// dimensions and keeping only the primary silently drops every row carried by
// the other.  Grouping them together lets the coalesced-point path merge them,
// and that path still throws if two columns claim the same row (i.e. if they
// really are redundant and there is no basis for preferring one).
struct EntryChunkFn {
  bool operator()(const ArrayIndex::Entry& a, const ArrayIndex::Entry& b) const
  {
    return a.array_name == b.array_name && a.array_type == b.array_type &&
           a.data_type == b.data_type;
  }
};

/******************************************************************************/
ArrayIndex::Layout group_name_to_layout(const std::string& name)
{
  if (name == "point") {
    return ArrayIndex::Layout::Point;
  } else if (name == "chunk") {
    return ArrayIndex::Layout::Chunked;
  } else {
    return ArrayIndex::Layout::Unknown;
  }
}

/******************************************************************************/
bool ArrayIndex::Entry::needed_for_decoding() const
{
  switch (buffer_format) {
  case MzPeak::Schema::BufferFormat::Point:
    return true;
  case MzPeak::Schema::BufferFormat::ChunkStart:
    return true;
  case MzPeak::Schema::BufferFormat::ChunkEnd:
    // Projected, even though no value is read OUT of it: the chunked decoder
    // validates the decode against chunk_end (start <= end, chunks ascending
    // and non-overlapping, and the start==end==0 empty-chunk sentinel), and
    // every one of those checks is guarded by `ends != nullptr`.  Returning
    // false here left the column unprojected, so `ends` was ALWAYS null and the
    // entire validation silently never ran -- while the docs claimed it did.
    // The specification requires chunks to be ascending by chunk_start and
    // non-overlapping, so this is a conformance check, not an optimisation.
    return true;
  case MzPeak::Schema::BufferFormat::ChunkValues:
    return true;
  case MzPeak::Schema::BufferFormat::ChunkEncoding:
    return true;
  case MzPeak::Schema::BufferFormat::ChunkSecondary:
    return true;
  case MzPeak::Schema::BufferFormat::ChunkTransform:
    return true;
  }

  std::unreachable();
}

/******************************************************************************/
bool ArrayIndex::Entry::is_value_entry() const
{
  switch (buffer_format) {
  case MzPeak::Schema::BufferFormat::Point:
    return true;
  case MzPeak::Schema::BufferFormat::ChunkStart:
    return false;
  case MzPeak::Schema::BufferFormat::ChunkEnd:
    // Projected, even though no value is read OUT of it: the chunked decoder
    // validates the decode against chunk_end (start <= end, chunks ascending
    // and non-overlapping, and the start==end==0 empty-chunk sentinel), and
    // every one of those checks is guarded by `ends != nullptr`.  Returning
    // false here left the column unprojected, so `ends` was ALWAYS null and the
    // entire validation silently never ran -- while the docs claimed it did.
    // The specification requires chunks to be ascending by chunk_start and
    // non-overlapping, so this is a conformance check, not an optimisation.
    return true;
  case MzPeak::Schema::BufferFormat::ChunkValues:
    return true;
  case MzPeak::Schema::BufferFormat::ChunkEncoding:
    return false;
  case MzPeak::Schema::BufferFormat::ChunkSecondary:
    return true;
  case MzPeak::Schema::BufferFormat::ChunkTransform:
    return true;
  }

  std::unreachable();
}

/******************************************************************************/
bool ArrayIndex::Dimension::is_main_axis() const
{
  for (const auto& entry : entries) {
    switch (entry.buffer_format) {
    case MzPeak::Schema::BufferFormat::Point:
      return true;
    case MzPeak::Schema::BufferFormat::ChunkStart:
      return true;
    case MzPeak::Schema::BufferFormat::ChunkEnd:
      return true;
    case MzPeak::Schema::BufferFormat::ChunkValues:
      return true;
    case MzPeak::Schema::BufferFormat::ChunkEncoding:
      return true;
    case MzPeak::Schema::BufferFormat::ChunkSecondary:
      return false;
    case MzPeak::Schema::BufferFormat::ChunkTransform:
      continue; // Could be main or secondary.
    }
  }

  return false;
}

/******************************************************************************/
bool ArrayIndex::Dimension::needs_delta_model() const
{
  // Null-marking reconstruction applies to the sorting-rank-0 array (m/z) ONLY;
  // parallel nulls in intensity arrays are read as 0 (spec signal-data.md).
  //
  // Deliberately NOT keyed off the transform accession.  The PSI-MS CV defines
  // MS:1003902 as the m/z-interpolating variant of MS:1003901, but the mzPeak
  // reference implementation writes them with the opposite sense (its
  // NULL_INTERPOLATE is MS:1003901, NULL_ZERO is MS:1003902) and every bundled
  // fixture follows suit — tagging the INTENSITY array MS:1003902.  Trusting
  // the accession therefore ran intensity through the delta-model interpolator
  // and produced negative intensities at null positions.  sorting_rank is the
  // one signal both sides agree on.
  return std::ranges::any_of(entries, [](const auto& e) {
    return e.sorting_rank.has_value() && e.sorting_rank.value() == 0;
  });
}

/******************************************************************************/
Util::Type ArrayIndex::Dimension::type_or_throw() const
{
  std::optional<Util::Type> type = data_type.as_type();
  if (type.has_value()) return type.value();
  throw TypeError("dimension " + name + " does not have a type set!");
}

/******************************************************************************/
const ArrayIndex::Entry& ArrayIndex::Dimension::values_entry() const
{
  // There are a few buffer formats that indicate that an entry is
  // definitely the column that stores dimension values.  However,
  // some of them (i.e. `ChunkTransform`) are ambitious so we need to
  // consider them after all other entries have been considered.
  //
  // We don't assume the `entries` vector is in an particular order
  // here.
  Entry const* chunk_transform = nullptr;
  Entry const* chunk_values = nullptr;

  for (const auto& entry : entries) {
    switch (entry.buffer_format) {
    case MzPeak::Schema::BufferFormat::Point:
      return entry;
    case MzPeak::Schema::BufferFormat::ChunkStart:
      continue;
    case MzPeak::Schema::BufferFormat::ChunkEnd:
      continue;
    case MzPeak::Schema::BufferFormat::ChunkValues:
      chunk_values = &entry;
      continue;
    case MzPeak::Schema::BufferFormat::ChunkEncoding:
      continue;
    case MzPeak::Schema::BufferFormat::ChunkSecondary:
      return entry;
    case MzPeak::Schema::BufferFormat::ChunkTransform:
      chunk_transform = &entry;
      continue;
    }
  }

  if (chunk_transform != nullptr) {
    return *chunk_transform;
  } else if (chunk_values != nullptr) {
    return *chunk_values;
  } else {
    std::string msg("dimension " + name + " lacks a data values column");
    throw InvalidFormatError(msg);
  }
}

/******************************************************************************/
std::optional<ArrayIndex::Entry>
ArrayIndex::Dimension::entry_with(BufferFormat format) const
{
  auto it = std::ranges::find(entries, format, &ArrayIndex::Entry::buffer_format);

  if (it == entries.end()) {
    return {};
  } else {
    return *it;
  }
}

/******************************************************************************/
ArrayIndex::ArrayIndex(EntityType entity_type, const json::object& obj)
    : entity_type_(entity_type)
    , prefix_(obj.at("prefix").as_string())
    , layout_(group_name_to_layout(prefix_))
    , entries_()
    , num_entities_()
{
  auto entries = obj.find("entries");

  if (entries != obj.end() && entries->value().is_array()) {
    auto entries_ary(entries->value().as_array());
    entries_.reserve(entries_ary.size());

    for (const auto& entry_obj : entries_ary) {
      if (entry_obj.is_object()) {
        const auto& eo(entry_obj.as_object());
        Entry entry;
        entry.array_name = eo.at("array_name").as_string();
        entry.buffer_format =
            buffer_format_from_string(eo.at("buffer_format").as_string());
        entry.context = EntityType(eo.at("context").as_string());
        entry.path = eo.at("path").as_string();

        // Every entry path is "<prefix>.<column>".  Deriving the column name by
        // cutting the prefix length off the front assumes that holds; when it
        // does not, substr throws a bare std::out_of_range whose message is
        // "basic_string" and which is not even an MzPeak exception, so the
        // caller learns nothing about what is wrong with the file.
        if (entry.path.size() <= prefix_.size() + 1 ||
            entry.path.compare(0, prefix_.size(), prefix_) != 0 ||
            entry.path[prefix_.size()] != '.') {
          throw JsonError("array index entry path '" + entry.path +
                          "' does not begin with the declared prefix '" + prefix_ +
                          ".'");
        }
        entry.name = entry.path.substr(prefix_.size() + 1);

        std::optional<Schema::CV> data_type_cv =
            Schema::CV::from_string(eo.at("data_type").as_string());

        if (data_type_cv.has_value()) {
          entry.data_type = Schema::PSI::DataType(data_type_cv.value());
        } else {
          std::string msg("invalid CV: ");
          msg += eo.at("data_type").as_string();
          throw JsonError(msg);
        }

        entry.array_type =
            PSI::array_type_from_string(eo.at("array_type").as_string());
        entry.unit = eo.at("unit").as_string();

        if (auto bp = eo.find("buffer_priority");
            bp != eo.end() && bp->value().is_string()) {
          entry.buffer_priority = bp->value().as_string() == "primary";
        }

        if (auto sr = eo.find("sorting_rank");
            sr != eo.end() && sr->value().is_number()) {
          if (sr->value().is_int64()) {
            entry.sorting_rank = sr->value().as_int64();
          } else {
            entry.sorting_rank = sr->value().as_uint64();
          }
        }

        if (auto dpi = eo.find("data_processing_id");
            dpi != eo.end() && dpi->value().is_string()) {
          entry.data_processing_id = dpi->value().as_string();
        }

        if (auto tr = eo.find("transform");
            tr != eo.end() && tr->value().is_string()) {
          std::optional<Schema::CV> maybe_cv =
              Schema::CV::from_string(tr->value().as_string());

          entry.transform = maybe_cv.and_then(
              [](const auto& cv) -> std::optional<Schema::PSI::Transform> {
                return Schema::PSI::Transform(cv);
              });
        }

        entries_.push_back(std::move(entry));
      }
    }
  }

  std::ranges::sort(entries_, EntryCmpFn());
}

/******************************************************************************/
EntityType ArrayIndex::entity_type() const { return entity_type_; }

/******************************************************************************/
const std::string& ArrayIndex::prefix() const { return prefix_; }

/******************************************************************************/
ArrayIndex::Layout ArrayIndex::layout() const { return layout_; }

/******************************************************************************/
const std::vector<ArrayIndex::Entry>& ArrayIndex::entries() const
{
  return entries_;
}

/******************************************************************************/
void ArrayIndex::num_entities(const std::optional<std::size_t>& ne)
{
  num_entities_ = ne;
}

/******************************************************************************/
std::optional<std::size_t> ArrayIndex::num_entities() const { return num_entities_; }

/******************************************************************************/
std::vector<ArrayIndex::Dimension> ArrayIndex::dimensions() const
{
  // Adjacent entries with the same key form one dimension. Written out rather
  // than as views::chunk_by | ranges::to: those are C++23 LIBRARY features
  // that Apple's clang 15/16 libc++ (the GitHub macOS runners) and MSVC do not
  // all have, and this is the only thing that needed them.
  std::vector<std::vector<Entry>> groups;
  const EntryChunkFn same_group;
  for (const Entry& entry : entries_) {
    if (groups.empty() || !same_group(groups.back().back(), entry))
      groups.emplace_back();
    groups.back().push_back(entry);
  }

  std::vector<Dimension> result;
  result.reserve(groups.size());

  for (auto& group : groups) {
    if (group.empty()) continue;
    auto& head = group[0];

    result.push_back({head.name, head.data_type, head.array_type,
                      head.buffer_priority, head.transform, std::move(group)});
  }

  return result;
}

/******************************************************************************/
std::optional<Schema::Column> ArrayIndex::entry_column(const Schema::GroupMap& map,
                                                       const Entry& col)
{
  auto it = map.find(prefix_);
  if (it == map.end()) return {};

  auto field = it->second->field(col.name);
  if (!field.has_value()) return {};

  return std::make_pair(it->second, field.value());
}

} // namespace MzPeak::Data
