// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include "IdentificationDataFileSupport.h"

#include <algorithm>
#include <bit>
#include <limits>
#include <numeric>
#include <set>

namespace OpenMS::Internal::IdentificationDataIO
{
[[noreturn]] void invalid(const std::string& message)
{ throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid native identification data", message); }
void check(const arrow::Status& status)
{
  if (! status.ok()) invalid(status.ToString());
}
void validateText(const std::string& input)
{
  // Strict UTF-8: reject overlong encodings, surrogate code points and values above U+10FFFF.
  const auto* bytes = reinterpret_cast<const unsigned char*>(input.data());
  for (Size i = 0; i < input.size();)
  {
    unsigned char lead = bytes[i++];
    if (lead < 0x80) continue;
    unsigned count;
    UInt32 code;
    if (lead >= 0xc2 && lead <= 0xdf)
    {
      count = 1;
      code = lead & 0x1f;
    }
    else if (lead >= 0xe0 && lead <= 0xef)
    {
      count = 2;
      code = lead & 0x0f;
    }
    else if (lead >= 0xf0 && lead <= 0xf4)
    {
      count = 3;
      code = lead & 0x07;
    }
    else
      invalid("Text is not valid UTF-8");
    if (count > input.size() - i) invalid("Truncated UTF-8 text");
    for (unsigned n = 0; n < count; ++n)
    {
      unsigned char part = bytes[i++];
      if ((part & 0xc0) != 0x80) invalid("Invalid UTF-8 continuation");
      code = (code << 6) | (part & 0x3f);
    }
    if ((count == 2 && code < 0x800) || (count == 3 && code < 0x10000) || (code >= 0xd800 && code <= 0xdfff) || code > 0x10ffff)
      invalid("Invalid UTF-8 code point");
  }
}
void validateOptions(const Options& options)
{
  if (! options.batch_rows || ! options.row_group_rows || ! options.batch_bytes || ! options.row_group_bytes || ! options.max_record_bytes)
    invalid("Buffer limits must be positive");
  if (options.batch_rows > static_cast<Size>(std::numeric_limits<int64_t>::max())
      || options.row_group_rows > static_cast<Size>(std::numeric_limits<int64_t>::max()))
    invalid("Row limit exceeds int64 range");
}
std::filesystem::path tablePath(const std::filesystem::path& root, const std::string& relative)
{
  validateText(relative);
  if (relative.find('\0') != std::string::npos) invalid("NUL is not permitted in table paths");
  const std::filesystem::path name(relative);
  if (name.empty() || name.is_absolute()) invalid("Expected a relative table path");
  for (const auto& component : name)
    if (component == ".." || component == ".") invalid("Invalid table path component");
  return root / name;
}
void appendText(arrow::ArrayBuilder& builder, const std::string& input)
{
  validateText(input);
  append<arrow::StringBuilder>(builder, input);
}
std::string text(const arrow::Array& array, int64_t row)
{
  if (array.IsNull(row)) invalid("Unexpected null text");
  return static_cast<const arrow::StringArray&>(array).GetString(row);
}
UInt32 Dictionary::add(const std::string& name, const DataValue& item)
{
  validateText(name);
  Key key {name, item.valueType(), item.getUnitType(), item.getUnit()};
  auto pos = ids_.find(key);
  if (pos != ids_.end()) return pos->second;
  if (values_.size() >= std::numeric_limits<UInt32>::max()) invalid("Metadata descriptor ID space exhausted");
  UInt32 id = static_cast<UInt32>(values_.size());
  ids_.emplace(std::move(key), id);
  values_.push_back({name, item.valueType(), item.getUnitType(), item.getUnit()});
  return id;
}
UInt32 Dictionary::find(const std::string& name, const DataValue& item) const
{
  auto pos = ids_.find({name, item.valueType(), item.getUnitType(), item.getUnit()});
  if (pos == ids_.end()) invalid("Metadata descriptor not registered before table writing");
  return pos->second;
}
void Dictionary::collect(const MetaInfoInterface& metadata)
{
  std::vector<std::string> names;
  metadata.getKeys(names);
  for (const auto& name : names)
    add(name, metadata.getMetaValue(name));
}
Json Dictionary::toJson() const
{
  Json out = Json::array();
  for (const auto& d : values_)
    out.push_back({{"name", d.name}, {"type", d.type}, {"unit_type", d.unit_type}, {"unit", d.unit}});
  return out;
}
void Dictionary::load(const Json& descriptors)
{
  if (! descriptors.is_array()) invalid("Metadata dictionary must be an array");
  ids_.clear();
  values_.clear();
  for (const auto& d : descriptors)
  {
    const auto type = integer<unsigned>(d.at("type"));
    const auto unit_type = integer<unsigned>(d.at("unit_type"));
    if (type >= DataValue::SIZE_OF_DATATYPE || unit_type > DataValue::OTHER) invalid("Unknown metadata type or unit ontology");
    Descriptor desc {d.at("name").get<std::string>(), static_cast<DataValue::DataType>(type), static_cast<DataValue::UnitType>(unit_type),
                     integer<Int32>(d.at("unit"))};
    validateText(desc.name);
    Key key {desc.name, type, unit_type, desc.unit};
    if (ids_.contains(key)) invalid("Duplicate metadata descriptor");
    if (values_.size() >= std::numeric_limits<UInt32>::max()) invalid("Too many metadata descriptors");
    ids_.emplace(std::move(key), static_cast<UInt32>(values_.size()));
    values_.push_back(std::move(desc));
  }
}
const Dictionary::Descriptor& Dictionary::at(UInt32 index) const
{
  if (index >= values_.size()) invalid("Unknown metadata descriptor ID");
  return values_[index];
}
std::shared_ptr<arrow::DataType> metadataType()
{
  return arrow::list(arrow::field("item",
                                  arrow::struct_({arrow::field("descriptor", arrow::uint32(), false), arrow::field("string_value", arrow::utf8()),
                                                  arrow::field("int_value", arrow::int64()), arrow::field("double_bits", arrow::uint64()),
                                                  arrow::field("string_list", arrow::list(arrow::field("item", arrow::utf8(), false))),
                                                  arrow::field("int_list", arrow::list(arrow::field("item", arrow::int64(), false))),
                                                  arrow::field("double_list_bits", arrow::list(arrow::field("item", arrow::uint64(), false)))}),
                                  false));
}
void appendMetadata(arrow::ArrayBuilder& builder, const MetaInfoInterface& metadata, const Dictionary& dictionary)
{
  auto& list = static_cast<arrow::ListBuilder&>(builder);
  auto& entry = *static_cast<arrow::StructBuilder*>(list.value_builder());
  check(list.Append());
  std::vector<std::string> names;
  metadata.getKeys(names);
  for (const auto& name : names)
  {
    const auto& item = metadata.getMetaValue(name);
    check(entry.Append());
    append<arrow::UInt32Builder>(*entry.field_builder(0), dictionary.find(name, item));
    const int active = item.valueType() == DataValue::EMPTY_VALUE ? 0 : static_cast<int>(item.valueType()) + 1;
    for (int field = 1; field < 7; ++field)
    {
      auto& child = *entry.field_builder(field);
      if (field != active)
      {
        check(child.AppendNull());
        continue;
      }
      switch (item.valueType())
      {
        case DataValue::STRING_VALUE:
          appendText(child, static_cast<std::string>(item));
          break;
        case DataValue::INT_VALUE:
          append<arrow::Int64Builder>(child, static_cast<int64_t>(item));
          break;
        case DataValue::DOUBLE_VALUE:
          append<arrow::UInt64Builder>(child, std::bit_cast<UInt64>(static_cast<double>(item)));
          break;
        case DataValue::STRING_LIST: {
          auto& values = static_cast<arrow::ListBuilder&>(child);
          check(values.Append());
          for (const auto& v : item.toStringList())
            appendText(*values.value_builder(), v);
          break;
        }
        case DataValue::INT_LIST: {
          auto& values = static_cast<arrow::ListBuilder&>(child);
          check(values.Append());
          for (const auto v : item.toIntList())
            append<arrow::Int64Builder>(*values.value_builder(), v);
          break;
        }
        case DataValue::DOUBLE_LIST: {
          auto& values = static_cast<arrow::ListBuilder&>(child);
          check(values.Append());
          for (const auto v : item.toDoubleList())
            append<arrow::UInt64Builder>(*values.value_builder(), std::bit_cast<UInt64>(v));
          break;
        }
        default:
          invalid("Unsupported metadata value type");
      }
    }
  }
}
void readMetadata(const arrow::Array& array, int64_t row, MetaInfoInterface& metadata, const Dictionary& dictionary)
{
  if (array.IsNull(row)) invalid("Null metadata list");
  const auto& list = static_cast<const arrow::ListArray&>(array);
  const auto& entries = static_cast<const arrow::StructArray&>(*list.values());
  std::set<std::string> names;
  for (int64_t i = list.value_offset(row), end = i + list.value_length(row); i < end; ++i)
  {
    if (entries.IsNull(i)) invalid("Null metadata entry");
    const auto& descriptor = dictionary.at(number<arrow::UInt32Array>(*entries.field(0), i));
    if (! names.insert(descriptor.name).second) invalid("Duplicate metadata name in record");
    int active = descriptor.type == DataValue::EMPTY_VALUE ? 0 : static_cast<int>(descriptor.type) + 1;
    for (int f = 1; f < 7; ++f)
      if (entries.field(f)->IsNull(i) == (f == active)) invalid("Metadata value does not match its descriptor type");
    DataValue item;
    switch (descriptor.type)
    {
      case DataValue::STRING_VALUE:
        item = text(*entries.field(1), i);
        break;
      case DataValue::INT_VALUE: {
        int64_t v = number<arrow::Int64Array>(*entries.field(2), i);
        if (v < std::numeric_limits<SignedSize>::min() || v > std::numeric_limits<SignedSize>::max())
          invalid("Metadata integer exceeds platform range");
        item = static_cast<long long>(v);
        break;
      }
      case DataValue::DOUBLE_VALUE:
        item = std::bit_cast<double>(number<arrow::UInt64Array>(*entries.field(3), i));
        break;
      case DataValue::STRING_LIST: {
        const auto& values = static_cast<const arrow::ListArray&>(*entries.field(4));
        StringList target;
        for (int64_t n = values.value_offset(i), end = n + values.value_length(i); n < end; ++n)
          target.push_back(text(*values.values(), n));
        item = target;
        break;
      }
      case DataValue::INT_LIST: {
        const auto& values = static_cast<const arrow::ListArray&>(*entries.field(5));
        IntList target;
        for (int64_t n = values.value_offset(i), end = n + values.value_length(i); n < end; ++n)
        {
          const int64_t v = number<arrow::Int64Array>(*values.values(), n);
          if (v < std::numeric_limits<Int>::min() || v > std::numeric_limits<Int>::max()) invalid("Integer-list entry exceeds model range");
          target.push_back(static_cast<Int>(v));
        }
        item = target;
        break;
      }
      case DataValue::DOUBLE_LIST: {
        const auto& values = static_cast<const arrow::ListArray&>(*entries.field(6));
        DoubleList target;
        for (int64_t n = values.value_offset(i), end = n + values.value_length(i); n < end; ++n)
          target.push_back(std::bit_cast<double>(number<arrow::UInt64Array>(*values.values(), n)));
        item = target;
        break;
      }
      case DataValue::EMPTY_VALUE:
        break;
      default:
        invalid("Unsupported metadata value type");
    }
    item.setUnitType(descriptor.unit_type);
    item.setUnit(descriptor.unit);
    metadata.setMetaValue(descriptor.name, item);
  }
}
Size metadataBytes(const MetaInfoInterface& metadata)
{
  Size size = 0;
  std::vector<std::string> names;
  metadata.getKeys(names);
  for (const auto& name : names)
  {
    const auto& item = metadata.getMetaValue(name);
    size += name.size() + 32;
    switch (item.valueType())
    {
      case DataValue::STRING_VALUE:
        size += static_cast<std::string>(item).size();
        break;
      case DataValue::STRING_LIST:
        for (const auto& s : item.toStringList())
          size += s.size() + 8;
        break;
      case DataValue::INT_LIST:
        size += item.toIntList().size() * 8;
        break;
      case DataValue::DOUBLE_LIST:
        size += item.toDoubleList().size() * 8;
        break;
      default:
        break;
    }
  }
  return size;
}
Json metadataJson(const MetaInfoInterface& metadata)
{
  Json out = Json::array();
  std::vector<std::string> names;
  metadata.getKeys(names);
  for (const auto& name : names)
  {
    validateText(name);
    const auto& v = metadata.getMetaValue(name);
    Json data;
    switch (v.valueType())
    {
      case DataValue::STRING_VALUE: {
        auto s = static_cast<std::string>(v);
        validateText(s);
        data = s;
        break;
      }
      case DataValue::INT_VALUE:
        data = static_cast<int64_t>(v);
        break;
      case DataValue::DOUBLE_VALUE:
        data = std::bit_cast<UInt64>(static_cast<double>(v));
        break;
      case DataValue::STRING_LIST: {
        auto strings = v.toStringList();
        for (const auto& s : strings)
          validateText(s);
        data = strings;
        break;
      }
      case DataValue::INT_LIST:
        data = v.toIntList();
        break;
      case DataValue::DOUBLE_LIST: {
        data = Json::array();
        for (double d : v.toDoubleList())
          data.push_back(std::bit_cast<UInt64>(d));
        break;
      }
      case DataValue::EMPTY_VALUE:
        break;
      default:
        invalid("Unsupported metadata value type");
    }
    out.push_back({{"name", name}, {"type", v.valueType()}, {"unit_type", v.getUnitType()}, {"unit", v.getUnit()}, {"value", std::move(data)}});
  }
  return out;
}
void readMetadataJson(const Json& json, MetaInfoInterface& metadata)
{
  if (! json.is_array()) invalid("Configuration metadata must be an array");
  std::set<std::string> names;
  for (const auto& entry : json)
  {
    const auto name = entry.at("name").get<std::string>();
    validateText(name);
    if (! names.insert(name).second) invalid("Duplicate configuration metadata key");
    const auto type = integer<unsigned>(entry.at("type"));
    const auto unit_type = integer<unsigned>(entry.at("unit_type"));
    if (unit_type > DataValue::OTHER) invalid("Unknown unit ontology");
    const auto& payload = entry.at("value");
    DataValue v;
    switch (type)
    {
      case DataValue::STRING_VALUE: {
        auto s = payload.get<std::string>();
        validateText(s);
        v = s;
        break;
      }
      case DataValue::INT_VALUE:
        v = static_cast<long long>(integer<int64_t>(payload));
        break;
      case DataValue::DOUBLE_VALUE:
        v = std::bit_cast<double>(integer<UInt64>(payload));
        break;
      case DataValue::STRING_LIST: {
        auto strings = payload.get<StringList>();
        for (const auto& s : strings)
          validateText(s);
        v = strings;
        break;
      }
      case DataValue::INT_LIST: {
        if (! payload.is_array()) invalid("Expected integer list");
        IntList values;
        for (const auto& i : payload)
          values.push_back(integer<Int>(i));
        v = values;
        break;
      }
      case DataValue::DOUBLE_LIST: {
        if (! payload.is_array()) invalid("Expected floating-point list");
        DoubleList values;
        for (const auto& d : payload)
          values.push_back(std::bit_cast<double>(integer<UInt64>(d)));
        v = values;
        break;
      }
      case DataValue::EMPTY_VALUE:
        if (! payload.is_null()) invalid("EMPTY_VALUE has a payload");
        break;
      default:
        invalid("Unknown configuration metadata type");
    }
    v.setUnitType(static_cast<DataValue::UnitType>(unit_type));
    v.setUnit(integer<Int32>(entry.at("unit")));
    metadata.setMetaValue(name, v);
  }
}
Json scoreJson(const ID::ScoreDefinition& s)
{
  return {{"name", s.name},
          {"accession", s.accession},
          {"higher_better", s.higher_better},
          {"scope", s.scope},
          {"software", s.software},
          {"software_version", s.software_version},
          {"parameters", metadataJson(s.parameters)},
          {"calibration", s.calibration},
          {"aggregation", s.aggregation}};
}
ID::ScoreDefinition readScoreJson(const Json& j)
{
  ID::ScoreDefinition s;
  s.name = j.at("name");
  s.accession = j.at("accession");
  s.higher_better = j.at("higher_better");
  unsigned scope = integer<unsigned>(j.at("scope"));
  if (scope > static_cast<unsigned>(ID::ScoreScope::OTHER)) invalid("Unknown score scope");
  s.scope = static_cast<ID::ScoreScope>(scope);
  s.software = j.at("software");
  s.software_version = j.at("software_version");
  readMetadataJson(j.at("parameters"), s.parameters);
  s.calibration = j.at("calibration");
  s.aggregation = j.at("aggregation");
  return s;
}

TableWriter::TableWriter(const std::filesystem::path& path, std::shared_ptr<arrow::Schema> schema, const Options& options):
    schema_(std::move(schema)),
    options_(options)
{
  validateOptions(options_);
  if (options.output)
  {
    delegate_ = &options.output->acquire(path, schema_, options);
    physical_path_ = delegate_->physical_path_;
    start_ = delegate_->rows();
    partition_ = options.partition;
    return;
  }
  sink_ = value(arrow::io::FileOutputStream::Open(path.string()));
  parquet::WriterProperties::Builder properties_builder;
  properties_builder.compression(parquet::Compression::ZSTD);
  // Floating-point dictionary equality can collapse -0/+0. Plain encoding keeps the
  // original bits while finite analytical columns remain directly queryable doubles.
  const auto preserve_float_bits = [&](const auto& self, const std::shared_ptr<arrow::DataType>& type, const std::string& path) -> void {
    if (type->id() == arrow::Type::DOUBLE || type->id() == arrow::Type::FLOAT) properties_builder.disable_dictionary(path);
    else if (type->id() == arrow::Type::STRUCT)
      for (const auto& field : type->fields())
        self(self, field->type(), path + "." + field->name());
    else if (type->id() == arrow::Type::LIST)
    {
      const auto& list = static_cast<const arrow::ListType&>(*type);
      self(self, list.value_type(), path + ".list.element");
      self(self, list.value_type(), path + ".list.item");
    }
  };
  for (const auto& field : schema_->fields())
    preserve_float_bits(preserve_float_bits, field->type(), field->name());
  auto properties = properties_builder.build();
  auto arrow_properties = parquet::ArrowWriterProperties::Builder().store_schema()->build();
  writer_ = value(parquet::arrow::FileWriter::Open(*schema_, arrow::default_memory_pool(), sink_, properties, arrow_properties));
  for (const auto& field : schema_->fields())
    builders_.push_back(value(arrow::MakeBuilder(field->type(), arrow::default_memory_pool())));
}
TableWriter::~TableWriter() = default;
void TableWriter::finishRow(Size bytes)
{
  if (logical_closed_) invalid("Appending to a closed table slice");
  if (delegate_)
  {
    append<arrow::UInt64Builder>(delegate_->column(schema_->num_fields()), partition_);
    delegate_->finishRow(bytes + sizeof(UInt64));
    ++total_;
    return;
  }
  if (bytes > options_.max_record_bytes) invalid("A record exceeds max_record_bytes");
  ++buffered_;
  ++total_;
  bytes_ += bytes;
  if (buffered_ >= options_.row_group_rows || bytes_ >= options_.row_group_bytes) flush_();
}
void TableWriter::flush_()
{
  if (! buffered_) return;
  std::vector<std::shared_ptr<arrow::Array>> arrays;
  for (auto& builder : builders_)
    arrays.push_back(value(builder->Finish()));
  auto batch = arrow::RecordBatch::Make(schema_, static_cast<int64_t>(buffered_), arrays);
  check(batch->ValidateFull());
  auto table = value(arrow::Table::FromRecordBatches({batch}));
  check(writer_->WriteTable(*table, static_cast<int64_t>(buffered_)));
  buffered_ = 0;
  bytes_ = 0;
}
void TableWriter::close()
{
  if (delegate_) { logical_closed_ = true; return; }
  if (! writer_) return;
  flush_();
  check(writer_->Close());
  writer_.reset();
  check(sink_->Close());
  sink_.reset();
}
Json TableWriter::reference() const
{
  return {{"path", physical_path_}, {"start", start_}, {"count", total_}, {"partition", partition_}};
}
TableWriter& WritePool::acquire(const std::filesystem::path& logical, const std::shared_ptr<arrow::Schema>& schema, const Options& options)
{
  const auto stem = logical.stem().string();
  const auto key = stem + "\n" + schema->ToString() + (stem == "matches" ? options.score_schema : "");
  auto& writer = writers_[key];
  if (!writer)
  {
    const auto filename = stem + "-" + std::to_string(variants_[stem]++) + ".parquet";
    const auto id = stem == "queries" || stem == "matches" || stem == "parents" ? "run_id" : "inference_id";
    auto physical = value(schema->AddField(schema->num_fields(), arrow::field(id, arrow::uint64(), false)));
    Options standalone(options);
    standalone.output.reset();
    standalone.input.reset();
    writer = std::make_unique<TableWriter>(root_ / filename, physical, standalone);
    writer->physical_path_ = filename;
  }
  return *writer;
}
void WritePool::close()
{
  for (auto& [key, writer] : writers_) writer->close();
}
std::shared_ptr<ReadPool::Entry> ReadPool::open(const std::filesystem::path& path)
{
  // Bound cached decoded groups and open files when many supplementary schemas occur.
  if (!files_.contains(path.string()) && files_.size() >= 16)
    for (auto it = files_.begin(); it != files_.end() && files_.size() >= 16;)
      if (it->second.use_count() == 1) it = files_.erase(it);
      else ++it;
  auto& entry = files_[path.string()];
  if (!entry)
  {
    entry = std::make_shared<Entry>();
    entry->source = value(arrow::io::ReadableFile::Open(path.string()));
    entry->reader = value(parquet::arrow::OpenFile(entry->source, arrow::default_memory_pool()));
    check(entry->reader->GetSchema(&entry->schema));
    const auto metadata = entry->reader->parquet_reader()->metadata();
    entry->starts.push_back(0);
    for (int i = 0; i < metadata->num_row_groups(); ++i)
      entry->starts.push_back(entry->starts.back() + metadata->RowGroup(i)->num_rows());
  }
  return entry;
}
TableReader::TableReader(const std::filesystem::path& root, const Json& reference,
                         const std::shared_ptr<arrow::Schema>& expected, const Options& options,
                         const std::vector<std::string>& columns): options_(options)
{
  validateOptions(options);
  auto pool = options.input ? options.input : std::make_shared<ReadPool>();
  entry_ = pool->open(tablePath(root, reference.at("path").get<std::string>()));
  start_ = integer<UInt64>(reference.at("start"));
  total_rows_ = integer<UInt64>(reference.at("count"));
  partition_ = integer<UInt64>(reference.at("partition"));
  if (start_ > entry_->starts.back() || total_rows_ > entry_->starts.back() - start_) invalid("Table slice outside physical row range");
  end_ = start_ + total_rows_;
  const auto schema = entry_->schema;
  if (schema->num_fields() != expected->num_fields() + 1) invalid("Unexpected shared table schema");
  id_column_ = schema->field(expected->num_fields())->name();
  if (id_column_ != "run_id" && id_column_ != "inference_id") invalid("Missing table partition column");
  if (!schema->Equals(*value(expected->AddField(expected->num_fields(), arrow::field(id_column_, arrow::uint64(), false))), false))
    invalid("Unexpected shared table schema");
  std::set<std::string> selected(columns.begin(), columns.end());
  for (const auto& name : selected)
    if (expected->GetFieldIndex(name) < 0) invalid("Unknown projected field " + name);
  selected.insert(id_column_);
  const auto* parquet_schema = entry_->reader->parquet_reader()->metadata()->schema();
  for (int i = 0; i < parquet_schema->num_columns(); ++i)
    if (columns.empty() || selected.contains(parquet_schema->Column(i)->path()->ToDotVector().front())) leaves_.push_back(i);
  if (total_rows_)
  {
    auto first = std::upper_bound(entry_->starts.begin(), entry_->starts.end(), start_);
    for (Size i = static_cast<Size>(first - entry_->starts.begin() - 1); i + 1 < entry_->starts.size() && entry_->starts[i] < end_; ++i)
      groups_.push_back(static_cast<int>(i));
  }
}
bool TableReader::next()
{
  ++row_;
  while (!batch_ || row_ >= batch_->num_rows())
  {
    if (batches_)
    {
      check(batches_->ReadNext(&batch_));
      row_ = 0;
      if (batch_ && batch_->num_rows()) return true;
      batches_.reset();
    }
    if (next_group_ == groups_.size()) return false;
    const int group = groups_[next_group_++];
    if (entry_->cached_group != group || entry_->cached_leaves != leaves_)
    {
      check(entry_->reader->ReadRowGroup(group, leaves_, &entry_->cached_table));
      check(entry_->cached_table->ValidateFull());
      entry_->cached_group = group;
      entry_->cached_leaves = leaves_;
    }
    const auto start = std::max(start_, entry_->starts[group]);
    const auto end = std::min(end_, entry_->starts[group + 1]);
    table_ = entry_->cached_table->Slice(start - entry_->starts[group], end - start);
    // Validate the ownership column even for projections that omit payloads.
    for (const auto& chunk : table_->GetColumnByName(id_column_)->chunks())
    {
      const auto& ids = static_cast<const arrow::UInt64Array&>(*chunk);
      for (int64_t i = 0; i < ids.length(); ++i)
        if (ids.IsNull(i) || ids.Value(i) != partition_) invalid("Table slice contains a different run/result partition");
    }
    batches_ = std::make_unique<arrow::TableBatchReader>(table_);
    Size batch_rows = options_.batch_rows;
    const auto metadata = entry_->reader->parquet_reader()->metadata()->RowGroup(group);
    long double bytes = 0;
    for (int leaf : leaves_) bytes += metadata->ColumnChunk(leaf)->total_uncompressed_size();
    if (bytes > 0 && metadata->num_rows() > 0)
      batch_rows = std::min(batch_rows, std::max<Size>(1, static_cast<Size>(static_cast<long double>(options_.batch_bytes) * metadata->num_rows() / bytes)));
    batches_->set_chunksize(static_cast<int64_t>(batch_rows));
  }
  return true;
}
const arrow::Array& TableReader::column(const std::string& name) const
{
  auto array = batch_->GetColumnByName(name);
  if (! array) invalid("Requested unprojected column " + name);
  return *array;
}
bool TableReader::hasColumn(const std::string& name) const
{ return batch_ && batch_->schema()->GetFieldIndex(name) >= 0; }
} // namespace OpenMS::Internal::IdentificationDataIO
