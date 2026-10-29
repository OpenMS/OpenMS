/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include "mzpeak/data/signals.h"

#include <chrono>

#include <arrow/record_batch.h>
#include <bit>
#include <memory>
#include <parquet/arrow/reader.h>
#include <parquet/statistics.h>

#include "mzpeak/data/array_index.h"
#include "mzpeak/schema/entity_type.h"
#include "mzpeak/util/executor.h"
#include "mzpeak/util/planner.h"
#include "mzpeak/util/projection.h"

namespace MzPeak::Data {

/******************************************************************************/
struct Signals::Impl {
  Impl(std::unique_ptr<Util::Parquet> parquet)
      : parquet_(std::move(parquet))
      , array_index_(parse_array_index())
  {
  }

  std::shared_ptr<ArrayIndex> parse_array_index() const;

  std::unique_ptr<Util::Parquet> parquet_;
  std::shared_ptr<ArrayIndex> array_index_;
};

/******************************************************************************/
std::shared_ptr<ArrayIndex> Signals::Impl::parse_array_index() const
{
  Util::Parquet::file_metadata_t fmd(parquet_->file_metadata());
  EntityType entity_type = parquet_->index_file().entity_type();

  std::string num_key(entity_type.metadata_count_key());
  std::optional<std::size_t> num_entities(parquet_->kv_size_t(fmd, num_key));

  std::string index_key(entity_type.array_index_name());
  auto index_str(parquet_->kv_string(fmd, index_key));
  if (!index_str.has_value()) throw InvalidFormatError("missing array_index");

  namespace json = boost::json;
  boost::system::error_code ec;
  json::value v(json::parse(index_str.value()));
  if (ec) throw MzPeak::JsonError(ec.message());

  if (v.is_object()) {
    auto ai = std::make_shared<ArrayIndex>(entity_type, v.as_object());
    ai->num_entities(num_entities);
    return ai;
  } else {
    throw JsonError("array_index should be a JSON object");
  }
}

/******************************************************************************/
Signals::Signals(std::unique_ptr<Util::Parquet> parquet)
    : impl_(std::make_unique<Impl>(std::move(parquet)))
{
}

/******************************************************************************/
Signals::~Signals() = default;

/******************************************************************************/
const std::shared_ptr<ArrayIndex>& Signals::array_index() const
{
  return impl_->array_index_;
}

/******************************************************************************/
std::size_t Signals::record_count() const
{
  auto ne(impl_->array_index_->num_entities());
  if (ne.has_value()) return *ne;

  // Fall back to scanning row-group max statistics on the entity-index column.
  // The index is 0-based and contiguous so COUNT = max_index + 1.
  const std::string prefix = impl_->array_index_->prefix();
  const std::string index_col = [&]() -> std::string {
    switch (impl_->array_index_->entity_type().type().value_or(
        Schema::EntityType::Spectrum)) {
    case Schema::EntityType::Chromatogram:
      return "chromatogram_index";
    case Schema::EntityType::WavelengthSpectrum:
      return "wavelength_spectrum_index";
    default:
      return "spectrum_index";
    }
  }();

  auto dest = impl_->parquet_->field(prefix, index_col);
  if (!dest.has_value()) {
    throw ParquetError("record_count: missing " + prefix + "." + index_col +
                       " column");
  }
  const int col = dest->second->absolute_index();
  auto fmd = impl_->parquet_->file_metadata();
  if (fmd->num_row_groups() == 0) return 0;

  int64_t max_signed = -1;
  for (int g = 0; g < fmd->num_row_groups(); ++g) {
    auto chunk = fmd->RowGroup(g)->ColumnChunk(col);
    if (!chunk->is_stats_set()) continue;
    auto stats = chunk->statistics();
    if (!stats || !stats->HasMinMax()) continue;
    auto typed = std::dynamic_pointer_cast<parquet::Int64Statistics>(stats);
    if (!typed) continue;
    max_signed = std::max(max_signed, typed->max());
  }

  if (max_signed < 0) {
    if (fmd->num_rows() == 0) return 0;
    throw ParquetError("record_count: no statistics for " + prefix + "." +
                       index_col);
  }
  // The index column is unsigned 64-bit stored with a signed physical type.
  return static_cast<std::size_t>(std::bit_cast<uint64_t>(max_signed)) + 1;
}

/******************************************************************************/
std::optional<Schema::Column> Signals::column(std::string_view name) const
{
  return impl_->parquet_->field(impl_->array_index_->prefix(), name);
}

/******************************************************************************/
std::optional<Schema::Column> Signals::column(const ArrayIndex::Entry& entry) const
{
  return impl_->array_index_->entry_column(*impl_->parquet_->groups(), entry);
}

/******************************************************************************/
std::optional<Schema::Column> Signals::column(const ArrayIndex::Dimension& dim,
                                              Schema::BufferFormat format) const
{
  return dim.entry_with(format).and_then(
      [this](const auto& entry) { return column(entry); });
}

/******************************************************************************/
const std::shared_ptr<Schema::GroupMap>& Signals::groups() const
{
  return impl_->parquet_->groups();
}

/******************************************************************************/
Util::Query::Builder Signals::index() const
{
  auto entity_type = impl_->array_index_->entity_type();
  auto field_name = entity_type.index_column_name();
  auto index_field = column(field_name);

  if (!index_field.has_value()) {
    throw InvalidFormatError("parquet file is missing the index column: " +
                             field_name);
  }

  return Util::Query::Builder(index_field.value());
}

/******************************************************************************/
std::unique_ptr<Util::Slice>
Signals::select(const std::vector<ArrayIndex::Dimension>& projection,
                const Util::Query& query)
{
  Util::Projection columns;

  for (const auto& dim : projection) {
    for (const auto& entry : dim.entries) {
      if (!entry.needed_for_decoding()) continue;

      auto field =
          impl_->array_index_->entry_column(*impl_->parquet_->groups(), entry);
      if (!field.has_value()) {
        throw InvalidFormatError("array entry not present in schema: " + entry.name);
      } else {
        columns.project(field.value());
      }
    }
  }

#ifdef MZPEAK_READ_COUNTERS
  const auto t_plan0 = std::chrono::steady_clock::now();
#endif
  Util::Planner planner = impl_->parquet_->planner(query);
#ifdef MZPEAK_READ_COUNTERS
  Util::count_ns_plan_ctor(std::chrono::duration_cast<std::chrono::nanoseconds>(
                               std::chrono::steady_clock::now() - t_plan0)
                               .count());
#endif
  auto plan = planner.plan();
#ifdef MZPEAK_READ_COUNTERS
  Util::count_ns_plan(std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - t_plan0)
                          .count());
#endif

  Util::Executor executor = impl_->parquet_->executor(columns);
  return executor.execute(plan);
}

} // namespace MzPeak::Data
