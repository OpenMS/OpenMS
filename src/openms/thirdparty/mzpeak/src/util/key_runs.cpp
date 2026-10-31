/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/key_runs.h"

#include <algorithm>
#include <arrow/array.h>
#include <arrow/buffer.h>
#include <arrow/record_batch.h>
#include <arrow/type.h>
#include <bit>
#include <parquet/api/reader.h>
#include <parquet/arrow/reader.h>
#include <parquet/statistics.h>

#include "mzpeak/exception.h"
#include "mzpeak/util/row_group_cache.h"

namespace MzPeak::Util {

namespace
{
  /// Is @p a strictly before @p b in the key's own order?
  bool before(bool is_unsigned, uint64_t a, uint64_t b)
  {
    return is_unsigned ? a < b : std::bit_cast<int64_t>(a) < std::bit_cast<int64_t>(b);
  }
} // namespace

/******************************************************************************/
std::optional<std::pair<int64_t, int64_t>>
KeyRuns::rows(const any_value_type& value) const
{
  // The value must be of the column's own C++ type, exactly as the executor's
  // binary search requires; anything else is not answered here at all.
  uint64_t wanted = 0;
  if (is_unsigned) {
    const auto* v = std::get_if<uint64_t>(&value);
    if (v == nullptr) return std::nullopt;
    wanted = *v;
  } else {
    const auto* v = std::get_if<int64_t>(&value);
    if (v == nullptr) return std::nullopt;
    wanted = std::bit_cast<uint64_t>(*v);
  }

  const auto it = std::lower_bound(
      key.begin(), key.end(), wanted,
      [this](uint64_t a, uint64_t b) { return before(is_unsigned, a, b); });
  if (it == key.end() || *it != wanted) return std::pair<int64_t, int64_t>{0, 0};
  const auto i = static_cast<std::size_t>(it - key.begin());
  return std::pair<int64_t, int64_t>{start[i], start[i + 1]};
}

/******************************************************************************/
namespace
{
/// The declared-sorted key leaf, when it has a shape runs can stand in for.
/// Structure only: whether it is sorted and null-free is for scan_key_runs().
std::optional<int32_t> key_leaf_shape(const parquet::RowGroupMetaData& group,
                                      const parquet::SchemaDescriptor& schema)
{
  // Only the FIRST sorting column, as the rest of the read path reads it: a
  // second one is sorted only within runs of the first.
  const std::vector<parquet::SortingColumn> sorting = group.sorting_columns();
  if (sorting.empty() || sorting.front().descending || sorting.front().nulls_first) {
    return std::nullopt;
  }
  const int32_t leaf = sorting.front().column_idx;
  if (leaf < 0 || leaf >= schema.num_columns() || schema.num_columns() < 2) {
    return std::nullopt;
  }

  const parquet::ColumnDescriptor* column = schema.Column(leaf);
  if (column->physical_type() != parquet::Type::INT64 ||
      column->max_repetition_level() != 0) {
    return std::nullopt;
  }
  const auto& logical = column->logical_type();
  if (logical && !logical->is_none() && !logical->is_int()) return std::nullopt;

  // Top-level, or a direct child of a struct that has other children: those
  // are the two shapes a placeholder can stand in for without moving anything.
  const std::size_t depth = column->path()->ToDotVector().size();
  if (depth == 2) {
    const auto* parent = column->schema_node()->parent();
    if (parent == nullptr || !parent->is_group() ||
        static_cast<const parquet::schema::GroupNode*>(parent)->field_count() < 2) {
      return std::nullopt;
    }
  } else if (depth != 1) {
    return std::nullopt;
  }
  return leaf;
}
} // namespace

/******************************************************************************/
std::optional<int32_t> run_key_leaf(const parquet::RowGroupMetaData& group,
                                    const parquet::SchemaDescriptor& schema)
{
  const std::optional<int32_t> leaf = key_leaf_shape(group, schema);
  if (!leaf.has_value()) return std::nullopt;

  // For SIZING, a nullable key must also be SAID to have no nulls.  The decode
  // takes the same word when the footer gives it, and reads the definition
  // levels when it does not, so a footer without statistics only makes the
  // charge an over-count, never an under-count.
  const parquet::ColumnDescriptor* column = schema.Column(*leaf);
  if (column->max_definition_level() > 0) {
    auto chunk = group.ColumnChunk(*leaf);
    if (!chunk || !chunk->is_stats_set()) return std::nullopt;
    std::shared_ptr<parquet::Statistics> stats = chunk->statistics();
    if (!stats || !stats->HasNullCount() || stats->null_count() != 0) {
      return std::nullopt;
    }
  }
  return leaf;
}

/******************************************************************************/
std::shared_ptr<KeyRuns> scan_key_runs(parquet::arrow::FileReader& reader, int32_t group)
{
  parquet::ParquetFileReader& file = *reader.parquet_reader();
  const std::shared_ptr<parquet::FileMetaData> metadata = file.metadata();
  const auto row_group = metadata->RowGroup(group);
  const std::optional<int32_t> leaf = key_leaf_shape(*row_group, *metadata->schema());
  if (!leaf.has_value()) return nullptr;
  const int64_t rows = row_group->num_rows();
  if (rows <= 0) return nullptr;

  // Where a batch will carry the key, and what Arrow decodes it to.
  std::shared_ptr<arrow::Schema> schema;
  if (!reader.GetSchema(&schema).ok() || !schema) return nullptr;
  const parquet::ColumnDescriptor* descriptor = metadata->schema()->Column(*leaf);
  const std::vector<std::string> path = descriptor->path()->ToDotVector();

  auto runs = std::make_shared<KeyRuns>();
  runs->leaf = *leaf;
  runs->top = schema->GetFieldIndex(path.front());
  if (runs->top < 0) return nullptr;
  runs->top_field = schema->field(runs->top);
  if (path.size() == 1) {
    runs->type = runs->top_field->type();
  } else {
    const auto& type = runs->top_field->type();
    if (type->id() != arrow::Type::STRUCT || type->num_fields() < 2) return nullptr;
    runs->child = static_cast<const arrow::StructType&>(*type).GetFieldIndex(path[1]);
    if (runs->child < 0) return nullptr;
    runs->type = type->field(runs->child)->type();
  }
  if (runs->type->id() == arrow::Type::UINT64) {
    runs->is_unsigned = true;
  } else if (runs->type->id() != arrow::Type::INT64) {
    return nullptr;
  }

  // A nullable key must have no null, and reading its VALUES alone does not
  // show that.  A chunk stores one value per non-null row, but a dictionary
  // page may end its indices in a bit-packed group zero-padded to eight, and
  // the padding reads back as the dictionary's first value: `rows` values
  // from a chunk with a null.  So the footer's null count decides when the
  // footer gives one Arrow accepts, as it does for sizing (run_key_leaf()):
  // above zero keeps the column, and zero is believed -- a footer that says
  // zero and holds a null is trusted here as it is there.  Without one (no
  // statistics, statistics is_stats_set() rejects or that do not decode, or
  // no null count in them), the definition levels are read with the values,
  // and a row without a value keeps the column.  Decoding the levels was a
  // third of the scan, so only a footer without a count pays for it.
  constexpr int64_t kChunk = 64 * 1024;
  std::vector<int16_t> levels;
  if (descriptor->max_definition_level() > 0) {
    std::shared_ptr<parquet::Statistics> stats;
    try {
      const auto chunk = row_group->ColumnChunk(*leaf);
      if (chunk && chunk->is_stats_set()) stats = chunk->statistics();
    } catch (const std::exception&) {
      // Statistics that do not decode (a min or max cut short) give no count,
      // and the levels decide.  Keeping the column instead would not spare
      // the group: Arrow's full decode of the key decodes them too, and fails.
    }
    if (stats && stats->HasNullCount()) {
      if (stats->null_count() != 0) return nullptr;
    } else {
      levels.resize(kChunk);
    }
  }

  // The column itself, 64 K values at a time into one reused buffer: it is
  // never materialised, and a sorted key's DELTA_BINARY_PACKED pages decode to
  // the buffer at memory speed.
  try {
    auto column = std::static_pointer_cast<parquet::Int64Reader>(
        file.RowGroup(group)->Column(*leaf));
    std::vector<int64_t> values(kChunk);
    const int64_t* data = values.data();

    int64_t row = 0;
    while (row < rows) {
      int64_t read = 0;
      const int64_t n = column->ReadBatch(std::min(kChunk, rows - row),
                                          levels.empty() ? nullptr : levels.data(), nullptr,
                                          values.data(), &read);
      // With levels, `n` counts rows and `read` only those with a value, so a
      // null makes the two differ.  A chunk that runs out short of `rows`
      // fails the row count below.
      if (n <= 0 || read != n) break;

      int64_t i = 0;
      if (runs->key.empty()) {
        runs->key.push_back(std::bit_cast<uint64_t>(data[0]));
        runs->start.push_back(row);
        i = 1;
      }
      uint64_t last = runs->key.back();
      while (i < n) {
        // Through the current run eight values at a time (a spectrum is
        // hundreds of rows), then one at a time to its end.
        while (i + 8 <= n) {
          uint64_t diff = 0;
          for (int j = 0; j < 8; ++j) diff |= std::bit_cast<uint64_t>(data[i + j]) ^ last;
          if (diff != 0) break;
          i += 8;
        }
        while (i < n && std::bit_cast<uint64_t>(data[i]) == last) ++i;
        if (i == n) break;

        const auto v = std::bit_cast<uint64_t>(data[i]);
        // Declared sorted, and is not: keep the column, and with it exactly
        // the behaviour such a file had before.
        if (before(runs->is_unsigned, v, last)) return nullptr;
        runs->key.push_back(v);
        runs->start.push_back(row + i);
        last = v;
        ++i;
      }
      row += n;
    }
    if (row != rows) return nullptr;
  } catch (const std::exception&) {
    // The values ran out under a null the footer denied, or the chunk is
    // damaged: the full decode meets the same fault, if it is one, and
    // reports it as it always has.
    return nullptr;
  }

  runs->start.push_back(rows);
  return runs;
}

/******************************************************************************/
std::shared_ptr<arrow::RecordBatch>
insert_key_placeholder(const std::shared_ptr<arrow::RecordBatch>& batch,
                       const KeyRuns& runs)
{
  const int64_t n = batch->num_rows();
  auto placeholder = std::make_shared<arrow::NullArray>(n);
  std::vector<std::shared_ptr<arrow::Array>> columns = batch->columns();

  if (runs.child < 0) {
    if (runs.top > static_cast<int>(columns.size())) return nullptr;
    columns.insert(columns.begin() + runs.top, placeholder);
    auto schema = batch->schema()->AddField(
        runs.top, arrow::field(runs.top_field->name(), arrow::null()));
    if (!schema.ok()) return nullptr;
    return arrow::RecordBatch::Make(*schema, n, std::move(columns));
  }

  if (runs.top >= static_cast<int>(columns.size())) return nullptr;
  const std::shared_ptr<arrow::Array>& column = columns[static_cast<std::size_t>(runs.top)];
  if (column->type_id() != arrow::Type::STRUCT || column->offset() != 0 ||
      column->length() != n) {
    return nullptr;
  }
  const auto& without = static_cast<const arrow::StructArray&>(*column);
  const auto& original = static_cast<const arrow::StructType&>(*runs.top_field->type());
  if (without.num_fields() + 1 != original.num_fields()) return nullptr;

  arrow::ArrayVector children = without.fields();
  children.insert(children.begin() + runs.child, placeholder);
  arrow::FieldVector fields = without.struct_type()->fields();
  fields.insert(fields.begin() + runs.child,
                arrow::field(original.field(runs.child)->name(), arrow::null()));
  auto type = arrow::struct_(fields);

  columns[static_cast<std::size_t>(runs.top)] = std::make_shared<arrow::StructArray>(
      type, n, children, without.data()->buffers[0], without.data()->null_count.load(), 0);
  auto schema = batch->schema()->SetField(
      runs.top, arrow::field(runs.top_field->name(), type, runs.top_field->nullable(),
                             runs.top_field->metadata()));
  if (!schema.ok()) return nullptr;
  return arrow::RecordBatch::Make(*schema, n, std::move(columns));
}

/******************************************************************************/
std::shared_ptr<const RowGroupBatches> restore_key_column(const RowGroupBatches& stripped)
{
  const KeyRuns& runs = *stripped.key_runs;
  auto full = std::make_shared<RowGroupBatches>();
  full->row_offset = stripped.row_offset;

  bool keyed = true;
  std::size_t run = 0;
  for (std::size_t bi = 0; bi < stripped.batches.size(); ++bi) {
    const std::shared_ptr<arrow::RecordBatch>& batch = stripped.batches[bi];
    const int64_t n = batch->num_rows();
    const int64_t base = stripped.row_offset[bi];

    auto allocated = arrow::AllocateBuffer(n * static_cast<int64_t>(sizeof(uint64_t)));
    if (!allocated.ok()) {
      throw ParquetError("restoring a key column: " + allocated.status().ToString());
    }
    std::shared_ptr<arrow::Buffer> values(std::move(*allocated));
    auto* out = reinterpret_cast<uint64_t*>(values->mutable_data());

    // Runs are strictly increasing in row, and batches are in row order, so
    // one cursor serves the whole group.
    while (run + 1 < runs.start.size() && runs.start[run + 1] <= base) ++run;
    for (int64_t r = base; r < base + n;) {
      const int64_t end = std::min(runs.start[run + 1], base + n);
      std::fill(out + (r - base), out + (end - base), runs.key[run]);
      if (end == runs.start[run + 1]) ++run;
      r = end;
    }

    auto key = arrow::MakeArray(arrow::ArrayData::Make(runs.type, n, {nullptr, values}, 0));
    std::vector<std::shared_ptr<arrow::Array>> columns = batch->columns();
    auto& slot = columns[static_cast<std::size_t>(runs.top)];
    if (runs.child < 0) {
      slot = key;
    } else {
      const auto& with = static_cast<const arrow::StructArray&>(*slot);
      arrow::ArrayVector children = with.fields();
      children[static_cast<std::size_t>(runs.child)] = key;
      slot = std::make_shared<arrow::StructArray>(runs.top_field->type(), n, children,
                                                  with.data()->buffers[0],
                                                  with.data()->null_count.load(), 0);
    }
    auto schema = batch->schema()->SetField(runs.top, runs.top_field);
    if (!schema.ok()) {
      throw ParquetError("restoring a key column: " + schema.status().ToString());
    }
    full->batches.push_back(arrow::RecordBatch::Make(*schema, n, std::move(columns)));

    // The per-batch key spans a full decode records (see decode_group_).
    if (n > 0) {
      full->key_first.push_back(static_cast<int64_t>(out[0]));
      full->key_last.push_back(static_cast<int64_t>(out[n - 1]));
    } else {
      keyed = false;
    }
  }

  if (keyed && !full->batches.empty()) {
    full->key_leaf = runs.leaf;
  } else {
    full->key_first.clear();
    full->key_last.clear();
  }
  return full;
}

} // namespace MzPeak::Util
