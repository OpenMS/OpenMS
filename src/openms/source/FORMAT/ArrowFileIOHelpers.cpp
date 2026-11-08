// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "ArrowFileIOHelpers.h"

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <cstdio>
#include <cstring>
#include <random>

namespace OpenMS::Internal::ArrowFileIO
{
// ==================== Export helpers ====================

/// Append all MetaValues (excluding specified keys) from a MetaInfoInterface to the struct builders.
void appendMetaValues(
  const MetaInfoInterface& mii,
  std::shared_ptr<arrow::StringBuilder>& name_b,
  std::shared_ptr<arrow::StringBuilder>& value_b,
  std::shared_ptr<arrow::StringBuilder>& type_b,
  std::shared_ptr<arrow::StructBuilder>& struct_b,
  const std::unordered_set<std::string>& excluded_keys)
{
  std::vector<std::string> keys;
  mii.getKeys(keys);
  for (const auto& key : keys)
  {
    if (excluded_keys.contains(key)) continue;
    const DataValue& val = mii.getMetaValue(key);
    (void)struct_b->Append();
    (void)name_b->Append(key);
    (void)value_b->Append(ArrowIOHelpers::metaValueString(val));
    switch (val.valueType())
    {
      case DataValue::INT_VALUE: (void)type_b->Append("int"); break;
      case DataValue::DOUBLE_VALUE: (void)type_b->Append("double"); break;
      case DataValue::STRING_VALUE: (void)type_b->Append("string"); break;
      case DataValue::INT_LIST: (void)type_b->Append("int_list"); break;
      case DataValue::DOUBLE_LIST: (void)type_b->Append("double_list"); break;
      case DataValue::STRING_LIST: (void)type_b->Append("string_list"); break;
      default: (void)type_b->Append("string"); break;
    }
  }
}

// ==================== JSON helpers for DataProcessing metadata ====================

/// Unescape a JSON string (reverse of ParquetFile::jsonEscape()).
std::string unescapeJsonString(const std::string& s)
{
  std::string result;
  result.reserve(s.size());
  for (size_t i = 0; i < s.size(); ++i)
  {
    if (s[i] == '\\' && i + 1 < s.size())
    {
      ++i;
      switch (s[i])
      {
        case '"':  result += '"'; break;
        case '\\': result += '\\'; break;
        case 'n':  result += '\n'; break;
        case 'r':  result += '\r'; break;
        case 't':  result += '\t'; break;
        case 'u':
          if (i + 4 < s.size())
          {
            unsigned int cp = 0;
            for (int j = 0; j < 4; ++j) // parse 4 hex digits
            {
              char ch = s[i + 1 + j];
              cp <<= 4;
              if (ch >= '0' && ch <= '9') cp |= (ch - '0');
              else if (ch >= 'a' && ch <= 'f') cp |= (ch - 'a' + 10);
              else if (ch >= 'A' && ch <= 'F') cp |= (ch - 'A' + 10);
            }
            // Encode as UTF-8
            if (cp < 0x80)
            {
              result += static_cast<char>(cp);
            }
            else if (cp < 0x800)
            {
              result += static_cast<char>(0xC0 | (cp >> 6));
              result += static_cast<char>(0x80 | (cp & 0x3F));
            }
            else
            {
              result += static_cast<char>(0xE0 | (cp >> 12));
              result += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
              result += static_cast<char>(0x80 | (cp & 0x3F));
            }
            i += 4;
          }
          break;
        default: result += s[i]; break;
      }
    }
    else
    {
      result += s[i];
    }
  }
  return result;
}

/// Serialize a vector of DataProcessing to a JSON string.
std::string serializeDataProcessing(const std::vector<DataProcessing>& dps)
{
  std::string json = "[";
  for (size_t i = 0; i < dps.size(); ++i)
  {
    const auto& dp = dps[i];
    if (i > 0) json += ",";
    json += "{";

    // software_name
    json += "\"software_name\":\"" + ParquetFile::jsonEscape(dp.getSoftware().getName()) + "\"";

    // software_version
    json += ",\"software_version\":\"" + ParquetFile::jsonEscape(dp.getSoftware().getVersion()) + "\"";

    // completion_time
    std::string time_str;
    if (dp.getCompletionTime().isValid())
    {
      time_str = dp.getCompletionTime().toString("yyyy-MM-ddThh:mm:ss");
    }
    json += ",\"completion_time\":\"" + ParquetFile::jsonEscape(time_str) + "\"";

    // actions
    json += ",\"actions\":[";
    const auto& actions = dp.getProcessingActions();
    bool first_action = true;
    for (const auto& action : actions)
    {
      if (!first_action) json += ",";
      json += "\"" + ParquetFile::jsonEscape(DataProcessing::processingActionToString(action)) + "\"";
      first_action = false;
    }
    json += "]";

    // metavalues
    json += ",\"metavalues\":" + serializeMetaValues(dp);

    json += "}";
  }
  json += "]";
  return json;
}

/// Skip whitespace in a JSON string starting at position pos.
void skipWhitespace(const std::string& s, size_t& pos)
{
  while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t' || s[pos] == '\n' || s[pos] == '\r'))
  {
    ++pos;
  }
}

/// Parse a JSON string value starting at position pos (must point to opening quote).
/// Returns the unescaped string value and advances pos past the closing quote.
std::string parseJsonString(const std::string& s, size_t& pos)
{
  if (pos >= s.size() || s[pos] != '"') return "";
  ++pos; // skip opening quote
  std::string raw;
  while (pos < s.size() && s[pos] != '"')
  {
    if (s[pos] == '\\' && pos + 1 < s.size())
    {
      raw += s[pos];
      raw += s[pos + 1];
      if (s[pos + 1] == 'u' && pos + 5 < s.size())
      {
        raw += StringUtils::substr(s, pos + 2, 4);
        pos += 6;
      }
      else
      {
        pos += 2;
      }
    }
    else
    {
      raw += s[pos];
      ++pos;
    }
  }
  if (pos < s.size()) ++pos; // skip closing quote
  return unescapeJsonString(raw);
}

void skipJsonValue(const std::string& s, size_t& pos)
{
  skipWhitespace(s, pos);
  if (pos >= s.size()) return;
  if (s[pos] == '"')
  {
    parseJsonString(s, pos);
    return;
  }
  if (s[pos] == '[' || s[pos] == '{')
  {
    // skip the nested array/object, honouring strings (which may contain brackets)
    int depth = 0;
    while (pos < s.size())
    {
      const char c = s[pos];
      if (c == '"') { parseJsonString(s, pos); continue; }
      if (c == '[' || c == '{') ++depth;
      else if (c == ']' || c == '}')
      {
        if (--depth == 0) { ++pos; return; }
      }
      ++pos;
    }
    return;
  }
  // number, true, false, null: up to the next separator
  while (pos < s.size() && s[pos] != ',' && s[pos] != '}' && s[pos] != ']') ++pos;
}

std::string parseJsonStringOrSkip(const std::string& s, size_t& pos)
{
  skipWhitespace(s, pos);
  if (pos < s.size() && s[pos] == '"') return parseJsonString(s, pos);
  skipJsonValue(s, pos); // not a string (number, literal, array, object): ignore, keep parsing
  return {};
}

/// Deserialize a JSON string into a vector of DataProcessing.
std::vector<DataProcessing> deserializeDataProcessing(const std::string& json, const char* context)
{
  std::vector<DataProcessing> result;
  if (json.empty()) return result;

  size_t pos = 0;
  skipWhitespace(json, pos);
  if (pos >= json.size() || json[pos] != '[') return result;
  ++pos; // skip '['

  while (pos < json.size())
  {
    skipWhitespace(json, pos);
    if (pos >= json.size() || json[pos] == ']') break;
    if (json[pos] == ',') { ++pos; continue; }
    if (json[pos] != '{') break;
    ++pos; // skip '{'

    DataProcessing dp;
    std::string software_name, software_version, completion_time;
    std::vector<std::string> action_names;

    // Parse key-value pairs in the object
    while (pos < json.size())
    {
      skipWhitespace(json, pos);
      if (pos >= json.size() || json[pos] == '}') { ++pos; break; }
      if (json[pos] == ',') { ++pos; continue; }

      // Parse key
      const size_t key_start = pos;
      std::string key = parseJsonString(json, pos);
      if (pos == key_start) { pos = json.size(); break; } // no string key: malformed, give up
      skipWhitespace(json, pos);
      if (pos < json.size() && json[pos] == ':') ++pos; // skip ':'
      skipWhitespace(json, pos);

      if (key == "software_name")
      {
        software_name = parseJsonStringOrSkip(json, pos);
      }
      else if (key == "software_version")
      {
        software_version = parseJsonStringOrSkip(json, pos);
      }
      else if (key == "completion_time")
      {
        completion_time = parseJsonStringOrSkip(json, pos);
      }
      else if (key == "actions")
      {
        // Parse array of strings
        if (pos < json.size() && json[pos] == '[')
        {
          ++pos; // skip '['
          while (pos < json.size())
          {
            skipWhitespace(json, pos);
            if (pos >= json.size() || json[pos] == ']') { ++pos; break; }
            if (json[pos] == ',') { ++pos; continue; }
            if (json[pos] == '"') { action_names.push_back(parseJsonString(json, pos)); continue; }
            const size_t value_start = pos;
            skipJsonValue(json, pos); // not a string: ignore
            if (pos == value_start) { pos = json.size(); break; } // stray '}' or similar: malformed, give up
          }
        }
        else
        {
          skipJsonValue(json, pos); // not an array: ignore, keep parsing
        }
      }
      else if (key == "metavalues")
      {
        parseMetaValueArray(json, pos, dp);
      }
      else
      {
        skipJsonValue(json, pos);
      }
    }

    // Apply parsed fields to DataProcessing
    dp.getSoftware().setName(software_name);
    dp.getSoftware().setVersion(software_version);
    if (!completion_time.empty())
    {
      dp.setCompletionTime(DateTime::fromString(completion_time, "yyyy-MM-ddThh:mm:ss"));
    }
    for (const auto& action_name : action_names)
    {
      try
      {
        dp.getProcessingActions().insert(DataProcessing::toProcessingAction(action_name));
      }
      catch (...)
      {
        OPENMS_LOG_WARN << context << ": Unknown processing action: " << action_name << std::endl;
      }
    }

    result.push_back(std::move(dp));
  }

  return result;
}

// ==================== MetaInfoInterface JSON helpers ====================

/// Serialize all MetaValues of a MetaInfoInterface to a JSON array string
/// (shape: `[{"name":..., "value":..., "type":...}, ...]`).
std::string serializeMetaValues(const MetaInfoInterface& mii)
{
  std::string json = "[";
  std::vector<std::string> keys;
  mii.getKeys(keys);
  bool first = true;
  for (const auto& key : keys)
  {
    if (!first) json += ",";
    const DataValue& val = mii.getMetaValue(key);
    std::string type_str;
    switch (val.valueType())
    {
      case DataValue::INT_VALUE: type_str = "int"; break;
      case DataValue::DOUBLE_VALUE: type_str = "double"; break;
      case DataValue::STRING_VALUE: type_str = "string"; break;
      case DataValue::INT_LIST: type_str = "int_list"; break;
      case DataValue::DOUBLE_LIST: type_str = "double_list"; break;
      case DataValue::STRING_LIST: type_str = "string_list"; break;
      default: type_str = "string"; break;
    }
    json += "{\"name\":\"" + ParquetFile::jsonEscape(std::string(key))
          + "\",\"value\":\"" + ParquetFile::jsonEscape(ArrowIOHelpers::metaValueString(val))
          + "\",\"type\":\"" + type_str + "\"}";
    first = false;
  }
  json += "]";
  return json;
}

/// Parse the JSON array of {name, value, type} objects that starts at @p pos (at its '[') into
/// @p target and advance @p pos past it. Per-entry type mismatches fall back to the raw string
/// value; unknown keys, a non-array value and entries that are not objects are skipped;
/// malformed input stops the parse.
void parseMetaValueArray(const std::string& json, size_t& pos, MetaInfoInterface& target)
{
  skipWhitespace(json, pos);
  if (pos >= json.size()) return;
  if (json[pos] != '[') { skipJsonValue(json, pos); return; } // not an array: ignore, keep parsing
  ++pos;

  while (pos < json.size())
  {
    skipWhitespace(json, pos);
    if (pos >= json.size() || json[pos] == ']') { ++pos; break; }
    if (json[pos] == ',') { ++pos; continue; }
    if (json[pos] != '{')
    {
      const size_t value_start = pos;
      skipJsonValue(json, pos); // an entry that is not an object: ignore it
      if (pos == value_start) { pos = json.size(); break; } // stray '}' or similar: malformed, give up
      continue;
    }
    ++pos;

    std::string mv_name, mv_value, mv_type;
    bool value_is_string = false; // an entry without a JSON-string value is dropped, not stored as ""
    while (pos < json.size())
    {
      skipWhitespace(json, pos);
      if (pos >= json.size() || json[pos] == '}') { ++pos; break; }
      if (json[pos] == ',') { ++pos; continue; }

      const size_t key_start = pos;
      std::string mk = parseJsonString(json, pos);
      if (pos == key_start) { pos = json.size(); break; } // no string key: malformed, give up
      skipWhitespace(json, pos);
      if (pos < json.size() && json[pos] == ':') ++pos;
      skipWhitespace(json, pos);

      if (mk == "name") mv_name = parseJsonStringOrSkip(json, pos);
      else if (mk == "value")
      {
        value_is_string = pos < json.size() && json[pos] == '"';
        mv_value = parseJsonStringOrSkip(json, pos);
      }
      else if (mk == "type") mv_type = parseJsonStringOrSkip(json, pos);
      else skipJsonValue(json, pos);
    }

    if (!mv_name.empty() && value_is_string)
    {
      if (mv_type == "int")
      {
        try { target.setMetaValue(mv_name, DataValue(std::stoll(mv_value))); }
        catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
      }
      else if (mv_type == "double" || mv_type == "float")
      {
        try { target.setMetaValue(mv_name, DataValue(StringUtils::toDouble(mv_value))); }
        catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
      }
      else if (mv_type == "int_list")
      {
        try
        {
          std::string s(mv_value);
          if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
          target.setMetaValue(mv_name, DataValue(ListUtils::create<Int>(s)));
        }
        catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
      }
      else if (mv_type == "double_list")
      {
        try
        {
          std::string s(mv_value);
          if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
          target.setMetaValue(mv_name, DataValue(ListUtils::create<double>(s)));
        }
        catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
      }
      else if (mv_type == "string_list")
      {
        try
        {
          std::string s(mv_value);
          if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
          auto sl = ListUtils::create<std::string>(s);
          for (auto& e : sl) { StringUtils::trim(e); }
          target.setMetaValue(mv_name, DataValue(sl));
        }
        catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
      }
      else
      {
        target.setMetaValue(mv_name, DataValue(mv_value));
      }
    }
  }
}

void deserializeMetaValues(const std::string& json, MetaInfoInterface& target)
{
  size_t pos = 0;
  parseMetaValueArray(json, pos, target);
}


/// Write an Arrow table to a Parquet file with QPX-style metadata.
bool writeParquetTable(
  std::shared_ptr<arrow::Table> table,
  const std::string& filename,
  const std::string& file_type,
  const ParquetWriteConfig& config,
  const char* context,
  const std::unordered_map<std::string, std::string>& extra_metadata)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << context << ": null table passed to writeParquetTable" << std::endl;
    return false;
  }

  // Generate RFC 4122 version-4 UUID
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint32_t> dist;
  uint8_t bytes[16];
  for (int i = 0; i < 4; ++i)
  {
    uint32_t r = dist(gen);
    std::memcpy(bytes + i * 4, &r, 4);
  }
  bytes[6] = (bytes[6] & 0x0F) | 0x40; // version 4
  bytes[8] = (bytes[8] & 0x3F) | 0x80; // variant 1
  char buf[37];
  std::snprintf(buf, sizeof(buf),
    "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
    bytes[0], bytes[1], bytes[2], bytes[3],
    bytes[4], bytes[5], bytes[6], bytes[7],
    bytes[8], bytes[9], bytes[10], bytes[11],
    bytes[12], bytes[13], bytes[14], bytes[15]);
  std::string uuid_str(buf);

  std::vector<std::pair<std::string, std::string>> md_pairs = {
    {"creator", "OpenMS"},
    {"file_type", file_type},
    {"creation_date", DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ")},
    {"uuid", uuid_str},
    {"software_provider", "OpenMS"}
  };
  for (const auto& kv : extra_metadata)
  {
    md_pairs.emplace_back(kv.first, kv.second);
  }
  auto metadata = std::make_shared<arrow::KeyValueMetadata>();
  for (const auto& p : md_pairs)
  {
    (void)metadata->Append(p.first, p.second);
  }
  table = table->ReplaceSchemaMetadata(metadata);

  // Open output file
  auto result = arrow::io::FileOutputStream::Open(filename);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << context << ": Failed to open file: " << filename << std::endl;
    return false;
  }
  const auto& outfile = *result;

  // Configure Parquet writer
  auto builder = parquet::WriterProperties::Builder();

  switch (config.compression)
  {
    case ParquetWriteConfig::Compression::NONE:
      builder.compression(arrow::Compression::UNCOMPRESSED);
      break;
    case ParquetWriteConfig::Compression::SNAPPY:
      builder.compression(arrow::Compression::SNAPPY);
      break;
    case ParquetWriteConfig::Compression::GZIP:
      builder.compression(arrow::Compression::GZIP);
      builder.compression_level(config.compression_level);
      break;
    case ParquetWriteConfig::Compression::LZ4:
      builder.compression(arrow::Compression::LZ4);
      break;
    case ParquetWriteConfig::Compression::ZSTD:
      builder.compression(arrow::Compression::ZSTD);
      builder.compression_level(config.compression_level);
      break;
  }

  builder.data_pagesize(config.data_page_size);

  if (config.write_statistics)
  {
    builder.enable_statistics();
  }
  else
  {
    builder.disable_statistics();
  }

  auto writer_props = builder.build();
  auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();

  // Write
  auto write_status = parquet::arrow::WriteTable(
    *table, arrow::default_memory_pool(), outfile,
    config.row_group_size, writer_props, arrow_props);

  // FileOutputStream::Open above already created (and truncated) the file, so any failure from
  // here on leaves a partial .parquet behind -- and a truncated Parquet file has no footer, so a
  // reader reports it as corrupt rather than as the smaller table it looks like. Close before
  // removing: on Windows an open handle blocks the unlink.
  const auto abandon = [&](const std::string& what)
  {
    OPENMS_LOG_ERROR << context << ": " << what << std::endl;
    (void)outfile->Close();
    if (!File::remove(filename))
    {
      OPENMS_LOG_ERROR << context << ": Failed to remove incomplete output "
                       << filename << std::endl;
    }
    return false;
  };

  if (!write_status.ok())
  {
    return abandon("Failed to write Parquet: " + write_status.ToString());
  }

  // Close explicitly rather than leaving it to the destructor, which swallows the error: the
  // final flush is where a full disk surfaces, and reporting success there would hand back a
  // truncated file.
  auto close_status = outfile->Close();
  if (!close_status.ok())
  {
    return abandon("Failed to close " + filename + ": " + close_status.ToString());
  }

  return true;
}

// ==================== Import helpers ====================

/// Read a single Parquet file into an Arrow table.
std::shared_ptr<arrow::Table> readParquetTable(const std::string& filename, const char* context)
{
  auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
  if (!infile_result.ok())
  {
    OPENMS_LOG_ERROR << context << ": Failed to open file: " << filename << std::endl;
    return nullptr;
  }
  const auto& infile = *infile_result;

  auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
  if (!reader_result.ok())
  {
    OPENMS_LOG_ERROR << context << ": Failed to create Parquet reader for: " << filename << std::endl;
    return nullptr;
  }
  auto reader = std::move(reader_result.ValueOrDie());

  std::shared_ptr<arrow::Table> table;
  auto read_status = reader->ReadTable(&table);
  if (!read_status.ok())
  {
    OPENMS_LOG_ERROR << context << ": Failed to read table: " << read_status.ToString() << std::endl;
    return nullptr;
  }

  auto combined = table->CombineChunks(arrow::default_memory_pool());
  if (!combined.ok())
  {
    OPENMS_LOG_ERROR << context << ": Failed to combine chunks" << std::endl;
    return nullptr;
  }

  return *combined;
}
} // namespace OpenMS::Internal::ArrowFileIO
