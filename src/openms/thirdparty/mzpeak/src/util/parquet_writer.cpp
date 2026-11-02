/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/parquet_writer.h"

#include <algorithm>
#include <arrow/array/array_nested.h>
#include <arrow/array/builder_binary.h>
#include <arrow/array/builder_primitive.h>
#include <arrow/buffer.h>
#include <arrow/io/file.h>
#include <arrow/io/memory.h>
#include <arrow/result.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <arrow/util/key_value_metadata.h>
#include <cmath>
#include <limits>
#include <memory>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>
#include <parquet/types.h>

#include "mzpeak/exception.h"

namespace MzPeak::Util {

/******************************************************************************/
double minutes_from_seconds(double seconds)
{
  const double m = seconds / 60.0;
  if (m * 60.0 == seconds || !std::isfinite(m)) return m;
  const double up = std::nextafter(m, std::numeric_limits<double>::infinity());
  if (up * 60.0 == seconds) return up;
  const double down = std::nextafter(m, -std::numeric_limits<double>::infinity());
  if (down * 60.0 == seconds) return down;
  return m;
}

namespace {

/******************************************************************************/
// Throw a ParquetError if the given status is not OK.
void check(const arrow::Status& status, const std::string& what)
{
  if (!status.ok()) throw ParquetError(what + ": " + status.ToString());
}

/******************************************************************************/
// Build a typed primitive Arrow array from a vector of values.
template <typename Builder, typename T>
std::shared_ptr<arrow::Array> build_array(const std::vector<T>& values)
{
  Builder builder;
  check(builder.Reserve(static_cast<int64_t>(values.size())), "reserve array");
  check(builder.AppendValues(values), "append values");

  std::shared_ptr<arrow::Array> array;
  check(builder.Finish(&array), "finish array");
  return array;
}

/******************************************************************************/
// Build a nullable Arrow array from a vector of optionals.  Missing values
// become Arrow nulls; present values are appended as-is.
template <typename Builder, typename T>
std::shared_ptr<arrow::Array>
build_optional_array(const std::vector<std::optional<T>>& values)
{
  Builder builder;
  check(builder.Reserve(static_cast<int64_t>(values.size())), "reserve array");
  for (const auto& v : values) {
    if (v.has_value()) {
      check(builder.Append(static_cast<typename Builder::value_type>(*v)),
            "append value");
    } else {
      check(builder.AppendNull(), "append null");
    }
  }
  std::shared_ptr<arrow::Array> array;
  check(builder.Finish(&array), "finish array");
  return array;
}

/******************************************************************************/
// Number of DFS leaves under a type.  A primitive is one leaf; a group is the
// sum of its children.
//
// A LIST counts as its element's leaves, which is what Parquet does too: a
// `large_list<struct<a:int32, b:int32>>` written by Arrow occupies leaves 0 and
// 1, putting a following `index` column at leaf 2.  Verified against a written
// file rather than reasoned about, because the metadata tables carry exactly
// that shape and a miscount here declares the wrong column sorted.
int leaf_count(const arrow::DataType& type)
{
  if (type.num_fields() == 0) return 1;
  int leaves = 0;
  for (const auto& field : type.fields())
    leaves += leaf_count(*field->type());
  return leaves;
}

/******************************************************************************/
// True for the name this writer gives an entity index: `index` on a flat
// metadata table, `<entity>_index` on a signal table's point/chunk struct.
bool is_entity_index(const std::string& name)
{
  return name == "index" || name.ends_with("_index");
}

/******************************************************************************/
// The DFS leaf position of the entity index, or nothing when the table is not
// shaped the way this writer emits -- a leading primitive `index`, or a leading
// struct whose first child is one.
//
// Returning nothing is the point.  A caller that cannot identify the sorted
// leaf must decline to declare one rather than guess at zero.
std::optional<int> entity_index_leaf_of(const arrow::Schema& schema)
{
  int leaves_before = 0;

  for (int i = 0; i < schema.num_fields(); ++i) {
    const auto& field = schema.field(i);
    const auto& type = *field->type();

    // A flat table: the index is a top-level primitive column.
    if (type.num_fields() == 0) {
      if (is_entity_index(field->name())) return leaves_before;
      leaves_before += 1;
      continue;
    }

    // A signal table: one struct whose first child is the index.
    if (type.id() == arrow::Type::STRUCT) {
      const auto& first = type.field(0);
      if (first->type()->num_fields() == 0 && is_entity_index(first->name())) {
        return leaves_before;
      }
    }

    leaves_before += leaf_count(type);
  }

  return std::nullopt;
}

/******************************************************************************/
// Is the entity-index column actually ascending?
//
// Deriving WHICH leaf is the index is not the same as knowing it is sorted.
// Every table this writer builds happens to be emitted in ascending index
// order, but that is an ordering property of a dozen separate builders rather
// than something the schema guarantees -- and a reader may binary search on the
// declaration.  One pass over a column the writer already holds in memory
// settles it.
bool entity_index_is_ascending(const arrow::Table& table)
{
  if (table.num_rows() == 0) return true;
  if (table.num_columns() == 0) return false;

  auto combined = table.CombineChunks();
  if (!combined.ok()) return false;
  auto column = (*combined)->column(0);
  if (column->num_chunks() != 1) return false;

  std::shared_ptr<arrow::Array> values = column->chunk(0);
  if (values->type_id() == arrow::Type::STRUCT) {
    values = std::static_pointer_cast<arrow::StructArray>(values)->field(0);
  }

  // (by the type id: UInt64Array is a template instance, whose type_info can differ between the images of a process;
  // OpenMS change, see README.OpenMS.md)
  if (values->type_id() != arrow::Type::UINT64) return false;
  auto index = std::static_pointer_cast<arrow::UInt64Array>(values);
  if (index->null_count() != 0) return false; // nulls_first=false would be a lie

  for (int64_t i = 1; i < index->length(); ++i) {
    if (index->Value(i) < index->Value(i - 1)) return false;
  }
  return true;
}

/******************************************************************************/
// The `point` struct table of a spectra data/peaks file: one struct column
// whose three leaves mirror the reference file, each a nullable field of the
// matching type.  Shared by the one-shot writer and the streaming one so the
// two can never drift.
std::shared_ptr<arrow::Table>
point_spectra_table(const std::vector<uint64_t>& spectrum_index,
                    const std::vector<double>& mz,
                    const std::vector<float>& intensity)
{
  auto index_array(build_array<arrow::UInt64Builder>(spectrum_index));
  auto mz_array(build_array<arrow::DoubleBuilder>(mz));
  auto intensity_array(build_array<arrow::FloatBuilder>(intensity));
  arrow::FieldVector point_fields{
      arrow::field("spectrum_index", arrow::uint64(), /*nullable=*/true),
      arrow::field("mz", arrow::float64(), /*nullable=*/true),
      arrow::field("intensity", arrow::float32(), /*nullable=*/true),
  };
  auto point_result(arrow::StructArray::Make(
      {index_array, mz_array, intensity_array}, point_fields));
  if (!point_result.ok()) {
    throw ParquetError("build point struct: " + point_result.status().ToString());
  }
  std::shared_ptr<arrow::Array> point_array(point_result.ValueOrDie());
  auto schema(arrow::schema(
      {arrow::field("point", point_array->type(), /*nullable=*/true)}));
  return arrow::Table::Make(schema, {point_array});
}

// The project's standard writer properties: ZSTD, statistics, page index,
// and -- when @p entity_index_leaf names one -- a sorting column on it.
std::shared_ptr<parquet::WriterProperties>
standard_writer_properties(std::optional<int> entity_index_leaf,
                           std::optional<int64_t> max_row_group_length = std::nullopt)
{
  parquet::WriterProperties::Builder props_builder;
  props_builder.compression(arrow::Compression::ZSTD);
  props_builder.enable_statistics();
  props_builder.enable_write_page_index();
  // Parquet caps a row group at 1M rows by default and splits a larger
  // WriteTable across groups.  A streaming writer sizes its own groups and
  // hands over whole spectra, so it raises the cap and gets exactly the
  // groups it asked for instead of a 1M group plus a ragged remainder.
  if (max_row_group_length.has_value()) {
    props_builder.max_row_group_length(*max_row_group_length);
  }
  if (entity_index_leaf.has_value()) {
    props_builder.set_sorting_columns({parquet::SortingColumn{
        /*column_idx=*/*entity_index_leaf, /*descending=*/false,
        /*nulls_first=*/false}});
  }
  return props_builder.build();
}

void add_file_kv(parquet::arrow::FileWriter& writer,
                 const std::map<std::string, std::string>& file_kv)
{
  if (file_kv.empty()) return;
  std::vector<std::string> keys;
  std::vector<std::string> values;
  keys.reserve(file_kv.size());
  values.reserve(file_kv.size());
  for (const auto& [key, value] : file_kv) {
    keys.push_back(key);
    values.push_back(value);
  }
  auto kv(std::make_shared<arrow::KeyValueMetadata>(keys, values));
  check(writer.AddKeyValueMetadata(kv), "add key/value metadata");
}

// Write an Arrow table to a Parquet sink with the project's standard
// properties: ZSTD, statistics, page index, store_schema, a sorting column
// on the first leaf (the entity index), a bounded row-group size, and
// optional file-level key/value metadata.  Shared by the data and metadata
// table writers so their on-disk shape stays identical (the reference reader
// selects rows via the index column's page index in both).
void write_table_to_sink(const std::shared_ptr<arrow::io::OutputStream>& sink,
                         const std::shared_ptr<arrow::Table>& table,
                         const std::map<std::string, std::string>& file_kv,
                         std::optional<int64_t> max_row_group = std::nullopt)
{
  // Which DFS leaf is the entity index?
  //
  // A `SortingColumn` names a leaf by its position in the file's depth-first
  // leaf order, so this has to be counted, not assumed.  An earlier version
  // looked like it counted -- it inspected the schema and then assigned 0 on
  // every branch -- which read as a safeguard while being none.
  //
  // When the shape is not the one this writer emits, NO sorting column is
  // declared.  Declaring the wrong one is worse than declaring none: a reader
  // may believe it and binary search a column that is not sorted, which finds
  // one run and silently misses the rest.
  std::optional<int> entity_index_leaf = entity_index_leaf_of(*table->schema());

  // Deriving WHICH leaf is the index is not the same as knowing it is sorted.
  // If the data is not actually ascending, declare no sorting column rather than
  // a false one -- a reader may believe it and binary search unsorted data.
  if (entity_index_leaf.has_value() && !entity_index_is_ascending(*table)) {
    entity_index_leaf.reset();
  }

  auto writer_props(standard_writer_properties(entity_index_leaf));

  // Store the Arrow schema so the struct/types round-trip exactly.
  auto arrow_props(
      parquet::ArrowWriterProperties::Builder().store_schema()->build());

  auto writer_result(parquet::arrow::FileWriter::Open(
      *table->schema(), arrow::default_memory_pool(), sink, writer_props,
      arrow_props));
  if (!writer_result.ok()) {
    throw ParquetError("open parquet writer: " + writer_result.status().ToString());
  }
  std::unique_ptr<parquet::arrow::FileWriter> writer(
      std::move(writer_result).ValueOrDie());

  // Bound the row-group size: one giant row group defeats page/row-group
  // pruning and is memory-hungry for large files.  Keep it positive even
  // for an empty table (WriteTable rejects a zero chunk size).
  // Callers may request a smaller cap; a test needs several row groups without
  // writing millions of rows to get them.
  constexpr int64_t kMaxRowGroup = 1 << 20; // ~1M rows
  const int64_t cap = max_row_group.value_or(kMaxRowGroup);
  int64_t row_group_size =
      table->num_rows() > 0 ? std::min<int64_t>(table->num_rows(), cap) : cap;
  check(writer->WriteTable(*table, row_group_size), "write table");

  // Embed file-level key/value metadata after the data, before Close().
  add_file_kv(*writer, file_kv);

  check(writer->Close(), "close parquet writer");
}

/******************************************************************************/
// Read the bytes of an in-memory Parquet buffer sink into a string.
std::string
finish_to_string(const std::shared_ptr<arrow::io::BufferOutputStream>& sink)
{
  auto buffer_result(sink->Finish());
  if (!buffer_result.ok()) {
    throw ParquetError("finish in-memory buffer sink: " +
                       buffer_result.status().ToString());
  }
  std::shared_ptr<arrow::Buffer> buffer(buffer_result.ValueOrDie());
  return std::string(reinterpret_cast<const char*>(buffer->data()),
                     static_cast<std::size_t>(buffer->size()));
}

/******************************************************************************/
// Build a nullable string array from a vector of optionals.  The numeric
// build_optional_array cannot serve here: arrow::StringBuilder has no nested
// value_type, so the cast in that template does not compile for strings.
template <typename Builder>
std::shared_ptr<arrow::Array>
build_optional_string_array(const std::vector<std::optional<std::string>>& values)
{
  Builder builder;
  check(builder.Reserve(static_cast<int64_t>(values.size())), "reserve array");
  for (const auto& v : values) {
    if (v.has_value()) {
      check(builder.Append(*v), "append value");
    } else {
      check(builder.AppendNull(), "append null");
    }
  }
  std::shared_ptr<arrow::Array> array;
  check(builder.Finish(&array), "finish array");
  return array;
}

/******************************************************************************/
// Build the point-layout spectra table and write it to the given Arrow
// output sink.  Shared by the file-path and in-memory-buffer entry points
// so the schema, writer properties and metadata handling stay identical.
void write_point_spectra_data_to_sink(
    const std::shared_ptr<arrow::io::OutputStream>& sink,
    const std::vector<uint64_t>& spectrum_index,
    const std::vector<double>& mz,
    const std::vector<float>& intensity,
    const std::map<std::string, std::string>& file_kv,
    std::optional<int64_t> max_row_group)
{
  if (spectrum_index.size() != mz.size() || mz.size() != intensity.size()) {
    throw ParquetError("write_point_spectra_data: input vectors must have "
                       "the same length");
  }

  write_table_to_sink(sink, point_spectra_table(spectrum_index, mz, intensity),
                      file_kv, max_row_group);
}

/******************************************************************************/
// Build the spectra metadata table (single `spectrum` struct column whose
// first child is the uint64 index) and write it to the given sink.
void write_spectra_metadata_to_sink(
    const std::shared_ptr<arrow::io::OutputStream>& sink,
    const std::vector<SpectrumMetaRow>& rows)
{
  std::vector<uint64_t> index;
  std::vector<std::string> id;
  std::vector<uint8_t> ms_level;
  std::vector<std::optional<double>> time;
  std::vector<std::optional<int>> polarity;
  std::vector<uint64_t> n_points;
  std::vector<uint64_t> n_peaks;
  std::vector<std::string> representation;
  index.reserve(rows.size());
  id.reserve(rows.size());
  ms_level.reserve(rows.size());
  time.reserve(rows.size());
  polarity.reserve(rows.size());
  n_points.reserve(rows.size());
  n_peaks.reserve(rows.size());
  representation.reserve(rows.size());
  for (const auto& r : rows) {
    index.push_back(r.index);
    id.push_back(r.id);
    ms_level.push_back(r.ms_level);
    // SpectrumMetaRow::retention_time is SECONDS (see writer.h); the spec's
    // `time` column is MINUTES (UO_0000031), and the reader multiplies by 60 on
    // the way back.  Storing seconds here made every round trip 60x too large.
    time.push_back(r.retention_time ? std::optional<double>(minutes_from_seconds(*r.retention_time))
                                    : std::nullopt);
    polarity.push_back(r.polarity);
    n_points.push_back(r.number_of_data_points);
    n_peaks.push_back(r.number_of_peaks);
    representation.push_back(r.representation);
  }

  // Emit the FLAT layout: plain column names at the top level, with the CV term
  // declared in the index's `column_mapping` instead of baked into the name.
  // `index` MUST come first — the reference reader accesses it positionally and
  // selects rows via its page index.
  auto schema(arrow::schema({
      arrow::field("index", arrow::uint64(), /*nullable=*/true),
      arrow::field("id", arrow::large_utf8(), /*nullable=*/true),
      arrow::field("ms_level", arrow::uint8(), /*nullable=*/true),
      arrow::field("time", arrow::float64(), /*nullable=*/true),
      arrow::field("scan_polarity", arrow::int8(), /*nullable=*/true),
      arrow::field("number_of_data_points", arrow::uint64(), /*nullable=*/true),
      arrow::field("number_of_peaks", arrow::uint64(), /*nullable=*/true),
      arrow::field("spectrum_representation", arrow::utf8(),
                   /*nullable=*/true),
  }));

  std::vector<std::shared_ptr<arrow::Array>> columns{
      build_array<arrow::UInt64Builder>(index),
      build_array<arrow::LargeStringBuilder>(id),
      build_array<arrow::UInt8Builder>(ms_level),
      build_optional_array<arrow::DoubleBuilder>(time),
      build_optional_array<arrow::Int8Builder>(polarity),
      build_array<arrow::UInt64Builder>(n_points),
      build_array<arrow::UInt64Builder>(n_peaks),
      build_array<arrow::StringBuilder>(representation),
  };

  auto table(arrow::Table::Make(schema, columns));

  // Upstream carries spectrum_count on the metadata members too, not just on
  // the signal file.  Without it the reference reader has no count to size the
  // spectrum set from.
  std::map<std::string, std::string> file_kv{
      {"spectrum_count", std::to_string(rows.size())}};
  write_table_to_sink(sink, table, file_kv);
}

} // namespace

/******************************************************************************/
void write_point_spectra_data(const std::string& path,
                              const std::vector<uint64_t>& spectrum_index,
                              const std::vector<double>& mz,
                              const std::vector<float>& intensity,
                              const std::map<std::string, std::string>& file_kv,
                              std::optional<int64_t> max_row_group)
{
  // Output sink.
  auto sink_result(arrow::io::FileOutputStream::Open(path));
  if (!sink_result.ok()) {
    throw ParquetError("open output file " + path + ": " +
                       sink_result.status().ToString());
  }
  std::shared_ptr<arrow::io::FileOutputStream> sink(sink_result.ValueOrDie());

  write_point_spectra_data_to_sink(sink, spectrum_index, mz, intensity, file_kv,
                                   max_row_group);
}

/******************************************************************************/
std::string
point_spectra_data_bytes(const std::vector<uint64_t>& spectrum_index,
                         const std::vector<double>& mz,
                         const std::vector<float>& intensity,
                         const std::map<std::string, std::string>& file_kv)
{
  // In-memory sink: the same table/schema/properties path as the file
  // writer, but the bytes are returned instead of landing on disk.
  auto sink_result(arrow::io::BufferOutputStream::Create());
  if (!sink_result.ok()) {
    throw ParquetError("create in-memory buffer sink: " +
                       sink_result.status().ToString());
  }
  std::shared_ptr<arrow::io::BufferOutputStream> sink(sink_result.ValueOrDie());

  write_point_spectra_data_to_sink(sink, spectrum_index, mz, intensity, file_kv,
                                   std::nullopt);

  return finish_to_string(sink);
}

/******************************************************************************/
void write_spectra_metadata(const std::string& path,
                            const std::vector<SpectrumMetaRow>& rows)
{
  auto sink_result(arrow::io::FileOutputStream::Open(path));
  if (!sink_result.ok()) {
    throw ParquetError("open output file " + path + ": " +
                       sink_result.status().ToString());
  }
  std::shared_ptr<arrow::io::FileOutputStream> sink(sink_result.ValueOrDie());

  write_spectra_metadata_to_sink(sink, rows);
}

/******************************************************************************/
std::string spectra_metadata_bytes(const std::vector<SpectrumMetaRow>& rows)
{
  auto sink_result(arrow::io::BufferOutputStream::Create());
  if (!sink_result.ok()) {
    throw ParquetError("create in-memory buffer sink: " +
                       sink_result.status().ToString());
  }
  std::shared_ptr<arrow::io::BufferOutputStream> sink(sink_result.ValueOrDie());

  write_spectra_metadata_to_sink(sink, rows);

  return finish_to_string(sink);
}

/******************************************************************************/
namespace {

/// Build the scans facet: one row per spectrum, carrying retention time.
std::shared_ptr<arrow::Table> scans_table(const std::vector<SpectrumMetaRow>& rows)
{
  std::vector<uint64_t> source_index, scan_index;
  std::vector<std::optional<double>> scan_start_time;
  std::vector<std::optional<double>> ion_mobility_value;
  std::vector<std::optional<std::string>> ion_mobility_type;
  for (const auto& r : rows) {
    source_index.push_back(r.index);
    scan_index.push_back(0);
    scan_start_time.push_back(
        r.retention_time
            // Seconds -> minutes, as for `time` above.
            ? std::optional<double>(minutes_from_seconds(*r.retention_time))
            : std::nullopt);
    ion_mobility_value.push_back(r.ion_mobility);
    ion_mobility_type.push_back(r.ion_mobility_type);
  }
  // ion_mobility_value / ion_mobility_type: the scan-level ion mobility columns the reader joins by
  // source_index (OpenMS addition; the reference reader ignores columns it does not know).
  auto schema(arrow::schema({
      arrow::field("source_index", arrow::uint64(), true),
      arrow::field("scan_index", arrow::uint64(), true),
      arrow::field("scan_start_time", arrow::float64(), true),
      arrow::field("ion_mobility_value", arrow::float64(), true),
      arrow::field("ion_mobility_type", arrow::utf8(), true),
  }));
  return arrow::Table::Make(
      schema, {build_array<arrow::UInt64Builder>(source_index),
               build_array<arrow::UInt64Builder>(scan_index),
               build_optional_array<arrow::DoubleBuilder>(scan_start_time),
               build_optional_array<arrow::DoubleBuilder>(ion_mobility_value),
               build_optional_string_array<arrow::StringBuilder>(ion_mobility_type)});
}

/// Precursor facet: one row per (spectrum, precursor).  Column names are the
/// ones the reference writer uses -- the reader resolves
/// MS_1000827_isolation_window_target_mz to `isolation_window_target` through
/// its alias table and, failing that, through the index's column_mapping, so
/// both spellings the ecosystem produces land on these columns.
///
/// Emitted with zero rows for an MS1-only run: the reference reader opens the
/// member unconditionally and fails on an absent one.
std::shared_ptr<arrow::Table>
precursors_table(const std::vector<SpectrumMetaRow>& rows)
{
  std::vector<uint64_t> source_index, precursor_index;
  std::vector<std::optional<float>> target, lower, upper;
  for (const auto& r : rows) {
    for (std::size_t k = 0; k < r.precursors.size(); ++k) {
      const PrecursorData& p = r.precursors[k];
      source_index.push_back(r.index);
      precursor_index.push_back(static_cast<uint64_t>(k));
      target.push_back(p.isolation_target_mz);
      lower.push_back(p.isolation_lower_offset);
      upper.push_back(p.isolation_upper_offset);
    }
  }

  auto window(arrow::StructArray::Make(
      {build_optional_array<arrow::FloatBuilder>(target),
       build_optional_array<arrow::FloatBuilder>(lower),
       build_optional_array<arrow::FloatBuilder>(upper)},
      {"isolation_window_target", "isolation_window_lower_offset",
       "isolation_window_upper_offset"}));
  check(window.status(), "build isolation_window struct");

  auto schema(arrow::schema({
      arrow::field("source_index", arrow::uint64(), true),
      arrow::field("precursor_index", arrow::uint64(), true),
      arrow::field("isolation_window", (*window)->type(), true),
  }));
  return arrow::Table::Make(schema, {build_array<arrow::UInt64Builder>(source_index),
                                     build_array<arrow::UInt64Builder>(precursor_index),
                                     *window});
}

/// Selected-ion facet: one row per (spectrum, precursor, ion), attached to its
/// precursor by (source_index, precursor_index) -- the join the reader makes.
std::shared_ptr<arrow::Table>
selected_ions_table(const std::vector<SpectrumMetaRow>& rows)
{
  std::vector<uint64_t> source_index, precursor_index;
  std::vector<std::optional<double>> mz;
  std::vector<std::optional<int32_t>> charge;
  std::vector<std::optional<float>> intensity;
  for (const auto& r : rows) {
    for (std::size_t k = 0; k < r.precursors.size(); ++k) {
      for (const SelectedIonData& ion : r.precursors[k].selected_ions) {
        source_index.push_back(r.index);
        precursor_index.push_back(static_cast<uint64_t>(k));
        mz.push_back(ion.mz);
        charge.push_back(ion.charge ? std::optional<int32_t>(*ion.charge)
                                    : std::nullopt);
        intensity.push_back(ion.intensity);
      }
    }
  }
  auto schema(arrow::schema({
      arrow::field("source_index", arrow::uint64(), true),
      arrow::field("precursor_index", arrow::uint64(), true),
      arrow::field("selected_ion_mz", arrow::float64(), true),
      arrow::field("charge_state", arrow::int32(), true),
      arrow::field("peak_intensity", arrow::float32(), true),
  }));
  return arrow::Table::Make(
      schema, {build_array<arrow::UInt64Builder>(source_index),
               build_array<arrow::UInt64Builder>(precursor_index),
               build_optional_array<arrow::DoubleBuilder>(mz),
               build_optional_array<arrow::Int32Builder>(charge),
               build_optional_array<arrow::FloatBuilder>(intensity)});
}

std::string table_bytes(const std::shared_ptr<arrow::Table>& table,
                        const std::map<std::string, std::string>& file_kv)
{
  auto sink_result(arrow::io::BufferOutputStream::Create());
  if (!sink_result.ok()) {
    throw ParquetError("create buffer: " + sink_result.status().ToString());
  }
  auto sink(sink_result.ValueOrDie());
  write_table_to_sink(sink, table, file_kv);
  auto buf(sink->Finish());
  if (!buf.ok()) throw ParquetError("finish buffer: " + buf.status().ToString());
  return (*buf)->ToString();
}

void write_table_to_path(const std::string& path,
                         const std::shared_ptr<arrow::Table>& table,
                         const std::map<std::string, std::string>& file_kv)
{
  auto sink_result(arrow::io::FileOutputStream::Open(path));
  if (!sink_result.ok()) {
    throw ParquetError("open output file " + path + ": " +
                       sink_result.status().ToString());
  }
  write_table_to_sink(sink_result.ValueOrDie(), table, file_kv);
}

} // namespace

/******************************************************************************/
void write_spectra_metadata_facets(const std::string& dir,
                                   const std::vector<SpectrumMetaRow>& rows)
{
  std::map<std::string, std::string> kv{
      {"spectrum_count", std::to_string(rows.size())}};
  write_table_to_path(dir + "/spectra_metadata_scans.parquet", scans_table(rows),
                      kv);
  write_table_to_path(dir + "/spectra_metadata_precursors.parquet",
                      precursors_table(rows), kv);
  write_table_to_path(dir + "/spectra_metadata_selected_ions.parquet",
                      selected_ions_table(rows), kv);
}

/******************************************************************************/
std::array<std::string, 3>
spectra_metadata_facet_bytes(const std::vector<SpectrumMetaRow>& rows)
{
  std::map<std::string, std::string> kv{
      {"spectrum_count", std::to_string(rows.size())}};
  return {table_bytes(scans_table(rows), kv),
          table_bytes(precursors_table(rows), kv),
          table_bytes(selected_ions_table(rows), kv)};
}

/******************************************************************************/
namespace {

/// Build a `point` struct table from three parallel leaf arrays and write it.
/// Shared by both new entity types so their files differ only in column names
/// and types, never in writer properties.
void write_point_table_to_sink(const std::shared_ptr<arrow::io::OutputStream>& sink,
                               const std::string& index_name,
                               const std::shared_ptr<arrow::Array>& index_array,
                               const std::string& axis_name,
                               const std::shared_ptr<arrow::Array>& axis_array,
                               const std::shared_ptr<arrow::Array>& intensity_array,
                               const std::map<std::string, std::string>& file_kv)
{
  arrow::FieldVector point_fields{
      arrow::field(index_name, index_array->type(), /*nullable=*/true),
      arrow::field(axis_name, axis_array->type(), /*nullable=*/true),
      arrow::field("intensity", intensity_array->type(), /*nullable=*/true),
  };

  auto point_result(arrow::StructArray::Make(
      {index_array, axis_array, intensity_array}, point_fields));
  if (!point_result.ok()) {
    throw ParquetError("build point struct: " + point_result.status().ToString());
  }
  std::shared_ptr<arrow::Array> point_array(point_result.ValueOrDie());

  auto schema(arrow::schema(
      {arrow::field("point", point_array->type(), /*nullable=*/true)}));
  write_table_to_sink(sink, arrow::Table::Make(schema, {point_array}), file_kv);
}

/// Open a file sink or throw.
std::shared_ptr<arrow::io::OutputStream> file_sink(const std::string& path)
{
  auto result(arrow::io::FileOutputStream::Open(path));
  if (!result.ok()) {
    throw ParquetError("open output file " + path + ": " +
                       result.status().ToString());
  }
  return result.ValueOrDie();
}

/// Create an in-memory sink or throw.
std::shared_ptr<arrow::io::BufferOutputStream> buffer_sink()
{
  auto result(arrow::io::BufferOutputStream::Create());
  if (!result.ok()) {
    throw ParquetError("create in-memory buffer sink: " +
                       result.status().ToString());
  }
  return result.ValueOrDie();
}

void chromatograms_data_to_sink(const std::shared_ptr<arrow::io::OutputStream>& sink,
                                const std::vector<uint64_t>& chromatogram_index,
                                const std::vector<double>& time,
                                const std::vector<float>& intensity,
                                const std::map<std::string, std::string>& file_kv)
{
  if (chromatogram_index.size() != time.size() || time.size() != intensity.size()) {
    throw ParquetError("write_point_chromatograms_data: input vectors must have "
                       "the same length");
  }
  write_point_table_to_sink(sink, "chromatogram_index",
                            build_array<arrow::UInt64Builder>(chromatogram_index),
                            "time", build_array<arrow::DoubleBuilder>(time),
                            build_array<arrow::FloatBuilder>(intensity), file_kv);
}

void wavelength_data_to_sink(const std::shared_ptr<arrow::io::OutputStream>& sink,
                             const std::vector<uint64_t>& spectrum_index,
                             const std::vector<float>& wavelength,
                             const std::vector<float>& intensity,
                             const std::map<std::string, std::string>& file_kv)
{
  if (spectrum_index.size() != wavelength.size() ||
      wavelength.size() != intensity.size()) {
    throw ParquetError("write_point_wavelength_data: input vectors must have "
                       "the same length");
  }
  write_point_table_to_sink(sink, "wavelength_spectrum_index",
                            build_array<arrow::UInt64Builder>(spectrum_index),
                            "wavelength",
                            build_array<arrow::FloatBuilder>(wavelength),
                            build_array<arrow::FloatBuilder>(intensity), file_kv);
}

/// The chromatogram metadata table: flat columns, split layout.  Column names
/// are plain and their CV terms are declared in the index's column_mapping,
/// which is what the current reference reader resolves through.
std::shared_ptr<arrow::Table>
chromatograms_metadata_table(const std::vector<ChromatogramMetaRow>& rows)
{
  std::vector<uint64_t> index;
  std::vector<std::string> id;
  std::vector<std::optional<std::string>> type;
  std::vector<std::optional<int>> polarity;
  std::vector<uint64_t> npoints;

  for (const auto& r : rows) {
    index.push_back(r.index);
    id.push_back(r.id);
    type.push_back(r.chromatogram_type.empty()
                       ? std::nullopt
                       : std::optional<std::string>(r.chromatogram_type));
    polarity.push_back(r.polarity);
    npoints.push_back(r.number_of_data_points);
  }

  auto schema(arrow::schema({
      arrow::field("index", arrow::uint64(), true),
      arrow::field("id", arrow::large_utf8(), true),
      arrow::field("chromatogram_type", arrow::utf8(), true),
      arrow::field("scan_polarity", arrow::int8(), true),
      arrow::field("number_of_data_points", arrow::uint64(), true),
  }));

  return arrow::Table::Make(schema,
                            {build_array<arrow::UInt64Builder>(index),
                             build_array<arrow::LargeStringBuilder>(id),
                             build_optional_string_array<arrow::StringBuilder>(type),
                             build_optional_array<arrow::Int8Builder>(polarity),
                             build_array<arrow::UInt64Builder>(npoints)});
}

/// The wavelength metadata table: flat columns, split layout.
std::shared_ptr<arrow::Table>
wavelength_metadata_table(const std::vector<WavelengthMetaRow>& rows)
{
  std::vector<uint64_t> index;
  std::vector<std::string> id;
  std::vector<std::optional<double>> time;
  std::vector<std::optional<std::string>> type;
  std::vector<std::optional<std::string>> representation;
  std::vector<uint64_t> npoints;
  std::vector<std::optional<double>> low;
  std::vector<std::optional<double>> high;
  std::vector<std::optional<double>> lambda_max;
  std::vector<std::optional<float>> bpi;
  std::vector<std::optional<float>> tic;

  for (const auto& r : rows) {
    index.push_back(r.index);
    id.push_back(r.id);
    time.push_back(r.time);
    type.push_back(r.spectrum_type.empty()
                       ? std::nullopt
                       : std::optional<std::string>(r.spectrum_type));
    representation.push_back(r.representation.empty()
                                 ? std::nullopt
                                 : std::optional<std::string>(r.representation));
    npoints.push_back(r.number_of_data_points);
    low.push_back(r.lowest_observed_wavelength);
    high.push_back(r.highest_observed_wavelength);
    lambda_max.push_back(r.lambda_max);
    bpi.push_back(r.base_peak_intensity);
    tic.push_back(r.total_ion_current);
  }

  auto schema(arrow::schema({
      arrow::field("index", arrow::uint64(), true),
      arrow::field("id", arrow::large_utf8(), true),
      arrow::field("time", arrow::float64(), true),
      arrow::field("spectrum_type", arrow::utf8(), true),
      arrow::field("spectrum_representation", arrow::utf8(), true),
      arrow::field("number_of_data_points", arrow::uint64(), true),
      arrow::field("lowest_observed_wavelength", arrow::float64(), true),
      arrow::field("highest_observed_wavelength", arrow::float64(), true),
      arrow::field("lambda_max", arrow::float64(), true),
      arrow::field("base_peak_intensity", arrow::float32(), true),
      arrow::field("total_ion_current", arrow::float32(), true),
  }));

  return arrow::Table::Make(
      schema, {build_array<arrow::UInt64Builder>(index),
               build_array<arrow::LargeStringBuilder>(id),
               build_optional_array<arrow::DoubleBuilder>(time),
               build_optional_string_array<arrow::StringBuilder>(type),
               build_optional_string_array<arrow::StringBuilder>(representation),
               build_array<arrow::UInt64Builder>(npoints),
               build_optional_array<arrow::DoubleBuilder>(low),
               build_optional_array<arrow::DoubleBuilder>(high),
               build_optional_array<arrow::DoubleBuilder>(lambda_max),
               build_optional_array<arrow::FloatBuilder>(bpi),
               build_optional_array<arrow::FloatBuilder>(tic)});
}

} // namespace

/******************************************************************************/
void write_point_chromatograms_data(
    const std::string& path,
    const std::vector<uint64_t>& chromatogram_index,
    const std::vector<double>& time,
    const std::vector<float>& intensity,
    const std::map<std::string, std::string>& file_kv)
{
  chromatograms_data_to_sink(file_sink(path), chromatogram_index, time, intensity,
                             file_kv);
}

/******************************************************************************/
std::string
point_chromatograms_data_bytes(const std::vector<uint64_t>& chromatogram_index,
                               const std::vector<double>& time,
                               const std::vector<float>& intensity,
                               const std::map<std::string, std::string>& file_kv)
{
  auto sink(buffer_sink());
  chromatograms_data_to_sink(sink, chromatogram_index, time, intensity, file_kv);
  return finish_to_string(sink);
}

/******************************************************************************/
void write_point_wavelength_data(const std::string& path,
                                 const std::vector<uint64_t>& spectrum_index,
                                 const std::vector<float>& wavelength,
                                 const std::vector<float>& intensity,
                                 const std::map<std::string, std::string>& file_kv)
{
  wavelength_data_to_sink(file_sink(path), spectrum_index, wavelength, intensity,
                          file_kv);
}

/******************************************************************************/
std::string
point_wavelength_data_bytes(const std::vector<uint64_t>& spectrum_index,
                            const std::vector<float>& wavelength,
                            const std::vector<float>& intensity,
                            const std::map<std::string, std::string>& file_kv)
{
  auto sink(buffer_sink());
  wavelength_data_to_sink(sink, spectrum_index, wavelength, intensity, file_kv);
  return finish_to_string(sink);
}

/******************************************************************************/
void write_chromatograms_metadata(const std::string& path,
                                  const std::vector<ChromatogramMetaRow>& rows,
                                  const std::map<std::string, std::string>& file_kv)
{
  write_table_to_sink(file_sink(path), chromatograms_metadata_table(rows), file_kv);
}

/******************************************************************************/
std::array<std::string, 2>
chromatogram_facet_bytes(const std::map<std::string, std::string>& file_kv)
{
  // ChromatogramData carries no precursor, so both facets have zero rows.
  return {table_bytes(precursors_table({}), file_kv),
          table_bytes(selected_ions_table({}), file_kv)};
}

/******************************************************************************/
std::string
chromatograms_metadata_bytes(const std::vector<ChromatogramMetaRow>& rows,
                             const std::map<std::string, std::string>& file_kv)
{
  auto sink(buffer_sink());
  write_table_to_sink(sink, chromatograms_metadata_table(rows), file_kv);
  return finish_to_string(sink);
}

/******************************************************************************/
void write_wavelength_metadata(const std::string& path,
                               const std::vector<WavelengthMetaRow>& rows,
                               const std::map<std::string, std::string>& file_kv)
{
  write_table_to_sink(file_sink(path), wavelength_metadata_table(rows), file_kv);
}

/******************************************************************************/
std::string wavelength_scans_bytes(const std::vector<WavelengthMetaRow>& rows,
                                   const std::map<std::string, std::string>& file_kv)
{
  std::vector<uint64_t> source_index;
  std::vector<uint64_t> scan_index;
  std::vector<std::optional<double>> start_time;

  for (const auto& r : rows) {
    source_index.push_back(r.index);
    scan_index.push_back(0);
    // Stored in minutes, like the primary column it mirrors.
    start_time.push_back(r.time.has_value() ? std::optional<double>(*r.time)
                                            : std::nullopt);
  }

  auto schema(arrow::schema({
      arrow::field("source_index", arrow::uint64(), true),
      arrow::field("scan_index", arrow::uint64(), true),
      arrow::field("scan_start_time", arrow::float64(), true),
  }));

  auto table(arrow::Table::Make(
      schema, {build_array<arrow::UInt64Builder>(source_index),
               build_array<arrow::UInt64Builder>(scan_index),
               build_optional_array<arrow::DoubleBuilder>(start_time)}));
  return table_bytes(table, file_kv);
}

/******************************************************************************/
std::string
wavelength_metadata_bytes(const std::vector<WavelengthMetaRow>& rows,
                          const std::map<std::string, std::string>& file_kv)
{
  auto sink(buffer_sink());
  write_table_to_sink(sink, wavelength_metadata_table(rows), file_kv);
  return finish_to_string(sink);
}

/******************************************************************************/
struct PointTableStream::Impl {
  std::shared_ptr<arrow::io::FileOutputStream> sink;
  std::unique_ptr<parquet::arrow::FileWriter> writer;
  int64_t rows = 0;
  uint64_t last_index = 0;
  bool closed = false;
};

PointTableStream::PointTableStream(const std::string& path)
    : impl_(std::make_unique<Impl>())
{
  auto sink_result(arrow::io::FileOutputStream::Open(path));
  if (!sink_result.ok()) {
    throw ParquetError("open point table for writing: " + path + ": " +
                       sink_result.status().ToString());
  }
  impl_->sink = sink_result.ValueOrDie();

  // The schema is fixed, so the sorting-column leaf can be derived from an
  // empty table up front; whether the data honours it is enforced per group.
  auto schema(point_spectra_table({}, {}, {})->schema());
  auto writer_props(standard_writer_properties(entity_index_leaf_of(*schema),
                                               std::numeric_limits<int64_t>::max()));
  auto arrow_props(
      parquet::ArrowWriterProperties::Builder().store_schema()->build());
  auto writer_result(parquet::arrow::FileWriter::Open(
      *schema, arrow::default_memory_pool(), impl_->sink, writer_props,
      arrow_props));
  if (!writer_result.ok()) {
    throw ParquetError("open parquet writer: " + writer_result.status().ToString());
  }
  impl_->writer = std::move(writer_result).ValueOrDie();
}

PointTableStream::~PointTableStream()
{
  // An unclosed stream is an abandoned file: release the handle, write no
  // footer, throw nothing.
  if (impl_ && !impl_->closed && impl_->sink) {
    (void)impl_->sink->Close();
  }
}

void PointTableStream::write_row_group(const std::vector<uint64_t>& spectrum_index,
                                       const std::vector<double>& mz,
                                       const std::vector<float>& intensity)
{
  if (spectrum_index.size() != mz.size() || mz.size() != intensity.size()) {
    throw ParquetError("PointTableStream: input vectors must have the same length");
  }
  if (spectrum_index.empty()) return;
  if (impl_->closed) throw ParquetError("PointTableStream: write after close");
  // The file declares the index sorted; keep that true across groups.
  if (impl_->rows > 0 && spectrum_index.front() < impl_->last_index) {
    throw ParquetError("PointTableStream: spectrum index " +
                       std::to_string(spectrum_index.front()) +
                       " arrives after " + std::to_string(impl_->last_index) +
                       "; groups must be appended in ascending index order");
  }
  for (std::size_t i = 1; i < spectrum_index.size(); ++i) {
    if (spectrum_index[i] < spectrum_index[i - 1]) {
      throw ParquetError("PointTableStream: spectrum index not ascending within a group");
    }
  }
  auto table(point_spectra_table(spectrum_index, mz, intensity));
  check(impl_->writer->WriteTable(*table, table->num_rows()), "write row group");
  impl_->rows += table->num_rows();
  impl_->last_index = spectrum_index.back();
}

int64_t PointTableStream::rows() const { return impl_->rows; }

void PointTableStream::close(const std::map<std::string, std::string>& file_kv)
{
  if (impl_->closed) return;
  add_file_kv(*impl_->writer, file_kv);
  check(impl_->writer->Close(), "close parquet writer");
  check(impl_->sink->Close(), "close point table file");
  impl_->closed = true;
}

} // namespace MzPeak::Util
