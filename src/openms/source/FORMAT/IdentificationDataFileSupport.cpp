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
  if (! writer_) return;
  flush_();
  check(writer_->Close());
  writer_.reset();
  check(sink_->Close());
  sink_.reset();
}
TableReader::TableReader(const std::filesystem::path& path,
                         const std::shared_ptr<arrow::Schema>& expected,
                         const Options& options,
                         const std::vector<std::string>& columns)
{
  validateOptions(options);
  source_ = value(arrow::io::ReadableFile::Open(path.string()));
  reader_ = value(parquet::arrow::OpenFile(source_, arrow::default_memory_pool()));
  std::shared_ptr<arrow::Schema> schema;
  check(reader_->GetSchema(&schema));
  if (! schema->Equals(*expected, false)) invalid("Unexpected schema in " + path.string());
  total_rows_ = static_cast<UInt64>(reader_->parquet_reader()->metadata()->num_rows());
  const auto metadata = reader_->parquet_reader()->metadata();
  const auto* parquet_schema = metadata->schema();
  std::vector<int> row_groups(reader_->num_row_groups());
  std::iota(row_groups.begin(), row_groups.end(), 0);
  std::set<std::string> selected(columns.begin(), columns.end());
  for (const auto& name : selected)
    if (schema->GetFieldIndex(name) < 0) invalid("Unknown projected field " + name);
  // Arrow's Parquet projection uses leaf indices, not top-level column indices.
  std::vector<int> leaves;
  for (int i = 0; i < parquet_schema->num_columns(); ++i)
    if (columns.empty() || selected.contains(parquet_schema->Column(i)->path()->ToDotVector().front())) leaves.push_back(i);
  // Footer byte densities provide a conservative batch-row target across the selected
  // row groups. Pages and a single large value may still exceed this target; this is
  // not a hard process-memory limit and does not include resident descriptors.
  Size batch_rows = options.batch_rows;
  for (int group_index : row_groups)
  {
    const auto group = metadata->RowGroup(group_index);
    long double bytes = 0;
    for (int leaf : leaves)
      bytes += group->ColumnChunk(leaf)->total_uncompressed_size();
    if (bytes > 0 && group->num_rows() > 0)
    {
      const long double target = static_cast<long double>(options.batch_bytes) * group->num_rows() / bytes;
      if (target < static_cast<long double>(batch_rows)) batch_rows = std::max<Size>(1, static_cast<Size>(target));
    }
  }
  reader_->set_batch_size(static_cast<int64_t>(batch_rows));
  batches_ = value(reader_->GetRecordBatchReader(row_groups, leaves));
}
bool TableReader::next()
{
  ++row_;
  while (! batch_ || row_ >= batch_->num_rows())
  {
    check(batches_->ReadNext(&batch_));
    row_ = 0;
    if (! batch_) return false;
    check(batch_->ValidateFull());
    if (batch_->num_rows()) return true;
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
