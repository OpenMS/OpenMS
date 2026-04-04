// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureMapArrowIO.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/CHEMISTRY/ProForma.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>
#include <parquet/properties.h>

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <random>
#include <unordered_set>

namespace OpenMS
{

namespace // anonymous
{
  // ==================== Export helpers ====================

  /// Append all MetaValues (excluding specified keys) from a MetaInfoInterface to the struct builders.
  void appendMetaValues_(
    const MetaInfoInterface& mii,
    std::shared_ptr<arrow::StringBuilder>& name_b,
    std::shared_ptr<arrow::StringBuilder>& value_b,
    std::shared_ptr<arrow::StringBuilder>& type_b,
    std::shared_ptr<arrow::StructBuilder>& struct_b,
    const std::unordered_set<std::string>& excluded_keys)
  {
    std::vector<String> keys;
    mii.getKeys(keys);
    for (const auto& key : keys)
    {
      if (excluded_keys.count(key)) continue;
      const DataValue& val = mii.getMetaValue(key);
      (void)struct_b->Append();
      (void)name_b->Append(key);
      (void)value_b->Append(val.toString());
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

  /// Escape a string for JSON (handles quotes, backslashes, and control characters).
  std::string escapeJsonString_(const std::string& s)
  {
    std::string result;
    result.reserve(s.size() + 8);
    for (char c : s)
    {
      switch (c)
      {
        case '"':  result += "\\\""; break;
        case '\\': result += "\\\\"; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default:
          if (static_cast<unsigned char>(c) < 0x20)
          {
            char buf[8];
            std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned char>(c));
            result += buf;
          }
          else
          {
            result += c;
          }
      }
    }
    return result;
  }

  /// Unescape a JSON string (reverse of escapeJsonString_).
  std::string unescapeJsonString_(const std::string& s)
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
  std::string serializeDataProcessing_(const std::vector<DataProcessing>& dps)
  {
    std::string json = "[";
    for (size_t i = 0; i < dps.size(); ++i)
    {
      const auto& dp = dps[i];
      if (i > 0) json += ",";
      json += "{";

      // software_name
      json += "\"software_name\":\"" + escapeJsonString_(dp.getSoftware().getName()) + "\"";

      // software_version
      json += ",\"software_version\":\"" + escapeJsonString_(dp.getSoftware().getVersion()) + "\"";

      // completion_time
      std::string time_str;
      if (dp.getCompletionTime().isValid())
      {
        time_str = dp.getCompletionTime().toString("yyyy-MM-ddThh:mm:ss");
      }
      json += ",\"completion_time\":\"" + escapeJsonString_(time_str) + "\"";

      // actions
      json += ",\"actions\":[";
      const auto& actions = dp.getProcessingActions();
      bool first_action = true;
      for (const auto& action : actions)
      {
        if (!first_action) json += ",";
        json += "\"" + escapeJsonString_(DataProcessing::processingActionToString(action)) + "\"";
        first_action = false;
      }
      json += "]";

      // metavalues
      json += ",\"metavalues\":[";
      std::vector<String> keys;
      dp.getKeys(keys);
      bool first_mv = true;
      for (const auto& key : keys)
      {
        if (!first_mv) json += ",";
        const DataValue& val = dp.getMetaValue(key);
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
        json += "{\"name\":\"" + escapeJsonString_(std::string(key))
              + "\",\"value\":\"" + escapeJsonString_(val.toString())
              + "\",\"type\":\"" + type_str + "\"}";
        first_mv = false;
      }
      json += "]";

      json += "}";
    }
    json += "]";
    return json;
  }

  /// Skip whitespace in a JSON string starting at position pos.
  void skipWhitespace_(const std::string& s, size_t& pos)
  {
    while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t' || s[pos] == '\n' || s[pos] == '\r'))
    {
      ++pos;
    }
  }

  /// Parse a JSON string value starting at position pos (must point to opening quote).
  /// Returns the unescaped string value and advances pos past the closing quote.
  std::string parseJsonString_(const std::string& s, size_t& pos)
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
          raw += s.substr(pos + 2, 4);
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
    return unescapeJsonString_(raw);
  }

  /// Deserialize a JSON string into a vector of DataProcessing.
  std::vector<DataProcessing> deserializeDataProcessing_(const std::string& json)
  {
    std::vector<DataProcessing> result;
    if (json.empty()) return result;

    size_t pos = 0;
    skipWhitespace_(json, pos);
    if (pos >= json.size() || json[pos] != '[') return result;
    ++pos; // skip '['

    while (pos < json.size())
    {
      skipWhitespace_(json, pos);
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
        skipWhitespace_(json, pos);
        if (pos >= json.size() || json[pos] == '}') { ++pos; break; }
        if (json[pos] == ',') { ++pos; continue; }

        // Parse key
        std::string key = parseJsonString_(json, pos);
        skipWhitespace_(json, pos);
        if (pos < json.size() && json[pos] == ':') ++pos; // skip ':'
        skipWhitespace_(json, pos);

        if (key == "software_name")
        {
          software_name = parseJsonString_(json, pos);
        }
        else if (key == "software_version")
        {
          software_version = parseJsonString_(json, pos);
        }
        else if (key == "completion_time")
        {
          completion_time = parseJsonString_(json, pos);
        }
        else if (key == "actions")
        {
          // Parse array of strings
          if (pos < json.size() && json[pos] == '[')
          {
            ++pos; // skip '['
            while (pos < json.size())
            {
              skipWhitespace_(json, pos);
              if (pos >= json.size() || json[pos] == ']') { ++pos; break; }
              if (json[pos] == ',') { ++pos; continue; }
              action_names.push_back(parseJsonString_(json, pos));
            }
          }
        }
        else if (key == "metavalues")
        {
          // Parse array of {name, value, type} objects
          if (pos < json.size() && json[pos] == '[')
          {
            ++pos; // skip '['
            while (pos < json.size())
            {
              skipWhitespace_(json, pos);
              if (pos >= json.size() || json[pos] == ']') { ++pos; break; }
              if (json[pos] == ',') { ++pos; continue; }
              if (json[pos] != '{') break;
              ++pos; // skip '{'

              std::string mv_name, mv_value, mv_type;
              while (pos < json.size())
              {
                skipWhitespace_(json, pos);
                if (pos >= json.size() || json[pos] == '}') { ++pos; break; }
                if (json[pos] == ',') { ++pos; continue; }

                std::string mk = parseJsonString_(json, pos);
                skipWhitespace_(json, pos);
                if (pos < json.size() && json[pos] == ':') ++pos;
                skipWhitespace_(json, pos);

                if (mk == "name") mv_name = parseJsonString_(json, pos);
                else if (mk == "value") mv_value = parseJsonString_(json, pos);
                else if (mk == "type") mv_type = parseJsonString_(json, pos);
                else parseJsonString_(json, pos); // skip unknown value
              }

              if (!mv_name.empty())
              {
                if (mv_type == "int")
                {
                  try { dp.setMetaValue(mv_name, DataValue(std::stoi(mv_value))); }
                  catch (...) { dp.setMetaValue(mv_name, DataValue(mv_value)); }
                }
                else if (mv_type == "double" || mv_type == "float")
                {
                  try { dp.setMetaValue(mv_name, DataValue(std::stod(mv_value))); }
                  catch (...) { dp.setMetaValue(mv_name, DataValue(mv_value)); }
                }
                else if (mv_type == "int_list")
                {
                  try
                  {
                    String s(mv_value);
                    if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
                    dp.setMetaValue(mv_name, DataValue(ListUtils::create<Int>(s)));
                  }
                  catch (...) { dp.setMetaValue(mv_name, DataValue(mv_value)); }
                }
                else if (mv_type == "double_list")
                {
                  try
                  {
                    String s(mv_value);
                    if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
                    dp.setMetaValue(mv_name, DataValue(ListUtils::create<double>(s)));
                  }
                  catch (...) { dp.setMetaValue(mv_name, DataValue(mv_value)); }
                }
                else if (mv_type == "string_list")
                {
                  try
                  {
                    String s(mv_value);
                    if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
                    auto sl = ListUtils::create<String>(s);
                    for (auto& e : sl) { e = e.trim(); }
                    dp.setMetaValue(mv_name, DataValue(sl));
                  }
                  catch (...) { dp.setMetaValue(mv_name, DataValue(mv_value)); }
                }
                else
                {
                  dp.setMetaValue(mv_name, DataValue(mv_value));
                }
              }
            }
          }
        }
        else
        {
          // Skip unknown value (string)
          parseJsonString_(json, pos);
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
          OPENMS_LOG_WARN << "FeatureMapArrowIO: Unknown processing action: " << action_name << std::endl;
        }
      }

      result.push_back(std::move(dp));
    }

    return result;
  }

  /// Write an Arrow table to a Parquet file with QPX-style metadata.
  bool writeArrowTableToParquet_(
    std::shared_ptr<arrow::Table> table,
    const String& filename,
    const std::string& file_type,
    const ParquetWriteConfig& config,
    const std::unordered_map<std::string, std::string>& extra_metadata = {})
  {
    if (!table)
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: null table passed to writeArrowTableToParquet_" << std::endl;
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
      {"qpx_version", "1.0"},
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
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to open file: " << filename << std::endl;
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

    if (!write_status.ok())
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to write Parquet: "
                       << write_status.ToString() << std::endl;
      return false;
    }

    return true;
  }

  /// Count total features recursively (top-level + all subordinates).
  Size countFeaturesRecursive_(const std::vector<Feature>& features)
  {
    Size count = 0;
    for (const auto& f : features)
    {
      ++count;
      count += countFeaturesRecursive_(f.getSubordinates());
    }
    return count;
  }

  // ==================== Import helpers ====================

  /// Read a single Parquet file into an Arrow table.
  std::shared_ptr<arrow::Table> readParquetTable_(const String& filename)
  {
    auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
    if (!infile_result.ok())
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to open file: " << filename << std::endl;
      return nullptr;
    }
    const auto& infile = *infile_result;

    auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    if (!reader_result.ok())
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to create Parquet reader for: " << filename << std::endl;
      return nullptr;
    }
    auto reader = std::move(reader_result.ValueOrDie());

    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    if (!read_status.ok())
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to read table: " << read_status.ToString() << std::endl;
      return nullptr;
    }

    auto combined = table->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok())
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to combine chunks" << std::endl;
      return nullptr;
    }

    return *combined;
  }

  /// Fetch a named column from a table, combining chunks if needed.
  /// Returns nullptr if column not found (logs error if required).
  std::shared_ptr<arrow::Array> getColumn_(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name,
    bool required = true)
  {
    auto column = table->GetColumnByName(name);
    if (!column)
    {
      if (required)
      {
        OPENMS_LOG_ERROR << "FeatureMapArrowIO: Missing required column '" << name << "'" << std::endl;
      }
      return nullptr;
    }
    if (column->num_chunks() == 0)
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Column '" << name << "' has no chunks" << std::endl;
      return nullptr;
    }
    if (column->num_chunks() == 1)
    {
      return column->chunk(0);
    }
    auto combined = arrow::Concatenate(column->chunks(), arrow::default_memory_pool());
    if (!combined.ok())
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to combine chunks for column '" << name << "'" << std::endl;
      return nullptr;
    }
    return *combined;
  }

  /// Get string value at a row, returning empty string if null.
  String getStringValue_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    if (!array || array->IsNull(row)) return "";
    return std::static_pointer_cast<arrow::StringArray>(array)->GetString(row);
  }

  /// Get double value at a row, returning default_val if null.
  double getDoubleValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, double default_val = 0.0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::DoubleArray>(array)->Value(row);
  }

  /// Get float value at a row, returning default_val if null.
  float getFloatValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, float default_val = 0.0f)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::FloatArray>(array)->Value(row);
  }

  /// Get int64 value at a row, returning default_val if null.
  int64_t getInt64Value_(const std::shared_ptr<arrow::Array>& array, int64_t row, int64_t default_val = 0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::Int64Array>(array)->Value(row);
  }

  /// Get int32 value at a row, returning default_val if null.
  int32_t getInt32Value_(const std::shared_ptr<arrow::Array>& array, int64_t row, int32_t default_val = 0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::Int32Array>(array)->Value(row);
  }

  /// Get boolean value at a row, returning default_val if null.
  bool getBoolValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, bool default_val = false)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::BooleanArray>(array)->Value(row);
  }

  /// Check if value at row is null.
  bool isNull_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    return !array || array->IsNull(row);
  }

  /// Read metavalues from a list<struct{name,value,value_type}> column at a given row.
  /// Sets them on the target MetaInfoInterface, excluding specified keys.
  void readMetaValues_(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    MetaInfoInterface& target,
    const std::unordered_set<std::string>& excluded_keys = {})
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
    auto value_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(1));
    auto type_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(2));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      std::string name = name_arr->GetString(i);
      if (excluded_keys.count(name)) continue;

      std::string value_str = value_arr->GetString(i);
      std::string type_str = type_arr->GetString(i);

      if (type_str == "int")
      {
        try { target.setMetaValue(name, static_cast<int>(std::stol(value_str))); }
        catch (...) { target.setMetaValue(name, value_str); }
      }
      else if (type_str == "double" || type_str == "float")
      {
        try { target.setMetaValue(name, std::stod(value_str)); }
        catch (...) { target.setMetaValue(name, value_str); }
      }
      else if (type_str == "int_list")
      {
        try
        {
          String s(value_str);
          if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
          target.setMetaValue(name, DataValue(ListUtils::create<Int>(s)));
        }
        catch (...) { target.setMetaValue(name, value_str); }
      }
      else if (type_str == "double_list")
      {
        try
        {
          String s(value_str);
          if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
          target.setMetaValue(name, DataValue(ListUtils::create<double>(s)));
        }
        catch (...) { target.setMetaValue(name, value_str); }
      }
      else if (type_str == "string_list")
      {
        try
        {
          String s(value_str);
          if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
          auto sl = ListUtils::create<String>(s);
          for (auto& e : sl) { e = e.trim(); }
          target.setMetaValue(name, DataValue(sl));
        }
        catch (...) { target.setMetaValue(name, value_str); }
      }
      else
      {
        target.setMetaValue(name, value_str);
      }
    }
  }

  /// Read convex hulls from a list<struct{hull_index, points: list<struct{x, y}>}> column at a given row.
  void readConvexHulls_(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    Feature& feature)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto hull_struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!hull_struct_arr || hull_struct_arr->length() == 0) return;

    auto points_list_arr = std::static_pointer_cast<arrow::ListArray>(hull_struct_arr->field(1));

    for (int64_t h = 0; h < hull_struct_arr->length(); ++h)
    {
      ConvexHull2D hull;
      ConvexHull2D::PointArrayType hull_points;

      auto point_struct_arr = std::static_pointer_cast<arrow::StructArray>(points_list_arr->value_slice(h));
      if (point_struct_arr && point_struct_arr->length() > 0)
      {
        auto x_arr = std::static_pointer_cast<arrow::DoubleArray>(point_struct_arr->field(0));
        auto y_arr = std::static_pointer_cast<arrow::DoubleArray>(point_struct_arr->field(1));

        hull_points.reserve(point_struct_arr->length());
        for (int64_t p = 0; p < point_struct_arr->length(); ++p)
        {
          hull_points.emplace_back(x_arr->Value(p), y_arr->Value(p));
        }
      }

      hull.setHullPoints(hull_points);
      feature.getConvexHulls().push_back(hull);
    }
  }

} // anonymous namespace


std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportFeaturesToArrow(
  const FeatureMap& feature_map)
{
  // -- Scalar column builders --
  arrow::Int64Builder unique_id_builder;
  arrow::Int64Builder parent_feature_id_builder;
  arrow::Int32Builder depth_builder;
  arrow::DoubleBuilder rt_builder, mz_builder;
  arrow::FloatBuilder intensity_builder;
  arrow::Int32Builder charge_builder;
  arrow::FloatBuilder overall_quality_builder, quality_rt_builder, quality_mz_builder;
  arrow::FloatBuilder width_builder;
  arrow::DoubleBuilder rt_bb_min_builder, rt_bb_max_builder, mz_bb_min_builder, mz_bb_max_builder;

  // -- convex_hulls: list<struct{hull_index: int32, points: list<struct{x: float64, y: float64}>}> --
  auto point_x_b = std::make_shared<arrow::DoubleBuilder>();
  auto point_y_b = std::make_shared<arrow::DoubleBuilder>();
  auto point_struct_type = arrow::struct_({
    arrow::field("x", arrow::float64()),
    arrow::field("y", arrow::float64())
  });
  auto point_struct_b = std::make_shared<arrow::StructBuilder>(
    point_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{point_x_b, point_y_b});
  auto points_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), point_struct_b);

  auto hull_index_b = std::make_shared<arrow::Int32Builder>();
  auto hull_struct_type = arrow::struct_({
    arrow::field("hull_index", arrow::int32()),
    arrow::field("points", arrow::list(point_struct_type))
  });
  auto hull_struct_b = std::make_shared<arrow::StructBuilder>(
    hull_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{hull_index_b, points_list_b});
  arrow::ListBuilder convex_hulls_builder(arrow::default_memory_pool(), hull_struct_b);

  // -- metavalues: list<struct{name: utf8, value: utf8, value_type: utf8}> --
  auto mv_name_b = std::make_shared<arrow::StringBuilder>();
  auto mv_value_b = std::make_shared<arrow::StringBuilder>();
  auto mv_type_b = std::make_shared<arrow::StringBuilder>();
  auto mv_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("value", arrow::utf8()),
    arrow::field("value_type", arrow::utf8())
  });
  auto mv_struct_b = std::make_shared<arrow::StructBuilder>(
    mv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{mv_name_b, mv_value_b, mv_type_b});
  arrow::ListBuilder metavalues_builder(arrow::default_memory_pool(), mv_struct_b);

  // Count total features (top-level + subordinates) for capacity reservation
  Size num_rows = countFeaturesRecursive_(feature_map.getData());

  // Reserve capacity for scalar builders
  arrow::Status status;

  #define RESERVE_OR_RETURN(builder_name) \
    status = builder_name.Reserve(num_rows); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "FeatureMapArrowIO: " #builder_name " Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  RESERVE_OR_RETURN(unique_id_builder)
  RESERVE_OR_RETURN(parent_feature_id_builder)
  RESERVE_OR_RETURN(depth_builder)
  RESERVE_OR_RETURN(rt_builder)
  RESERVE_OR_RETURN(mz_builder)
  RESERVE_OR_RETURN(intensity_builder)
  RESERVE_OR_RETURN(charge_builder)
  RESERVE_OR_RETURN(overall_quality_builder)
  RESERVE_OR_RETURN(quality_rt_builder)
  RESERVE_OR_RETURN(quality_mz_builder)
  RESERVE_OR_RETURN(width_builder)
  RESERVE_OR_RETURN(rt_bb_min_builder)
  RESERVE_OR_RETURN(rt_bb_max_builder)
  RESERVE_OR_RETURN(mz_bb_min_builder)
  RESERVE_OR_RETURN(mz_bb_max_builder)

  #undef RESERVE_OR_RETURN

  // Exclude FWHM from metavalues (it is already stored as a dedicated column or derived from convex hulls)
  static const std::unordered_set<std::string> excluded_mvs = {"FWHM"};

  // Recursive lambda to flatten features depth-first
  std::function<void(const Feature&, int64_t, int32_t)> appendFeature =
    [&](const Feature& feature, int64_t parent_id, int32_t depth)
  {
    // === unique_id (not nullable) ===
    (void)unique_id_builder.Append(static_cast<int64_t>(feature.getUniqueId()));

    // === parent_feature_id (nullable: null for top-level) ===
    if (depth == 0)
    {
      (void)parent_feature_id_builder.AppendNull();
    }
    else
    {
      (void)parent_feature_id_builder.Append(parent_id);
    }

    // === depth (not nullable) ===
    (void)depth_builder.Append(depth);

    // === rt (not nullable) ===
    (void)rt_builder.Append(feature.getRT());

    // === mz (not nullable) ===
    (void)mz_builder.Append(feature.getMZ());

    // === intensity (not nullable) ===
    (void)intensity_builder.Append(feature.getIntensity());

    // === charge (not nullable) ===
    (void)charge_builder.Append(static_cast<int32_t>(feature.getCharge()));

    // === overall_quality (not nullable) ===
    (void)overall_quality_builder.Append(feature.getOverallQuality());

    // === quality_rt (not nullable) ===
    (void)quality_rt_builder.Append(feature.getQuality(0));

    // === quality_mz (not nullable) ===
    (void)quality_mz_builder.Append(feature.getQuality(1));

    // === width (nullable: null if 0.0, since 0.0 is the unset default; round-trip correct) ===
    float w = feature.getWidth();
    if (w == 0.0f)
    {
      (void)width_builder.AppendNull();
    }
    else
    {
      (void)width_builder.Append(w);
    }

    // === bounding box (nullable: null if no convex hulls) ===
    const auto& hulls = feature.getConvexHulls();
    if (hulls.empty())
    {
      (void)rt_bb_min_builder.AppendNull();
      (void)rt_bb_max_builder.AppendNull();
      (void)mz_bb_min_builder.AppendNull();
      (void)mz_bb_max_builder.AppendNull();
    }
    else
    {
      auto bb = feature.getConvexHull().getBoundingBox();
      (void)rt_bb_min_builder.Append(bb.minPosition()[0]);
      (void)rt_bb_max_builder.Append(bb.maxPosition()[0]);
      (void)mz_bb_min_builder.Append(bb.minPosition()[1]);
      (void)mz_bb_max_builder.Append(bb.maxPosition()[1]);
    }

    // === convex_hulls (not nullable) ===
    (void)convex_hulls_builder.Append(); // begin list for this feature
    for (int32_t hull_idx = 0; hull_idx < static_cast<int32_t>(hulls.size()); ++hull_idx)
    {
      const auto& hull = hulls[hull_idx];
      const auto& points = hull.getHullPoints();

      (void)hull_struct_b->Append(); // begin struct for this hull
      (void)hull_index_b->Append(hull_idx);
      (void)points_list_b->Append(); // begin list of points
      for (const auto& pt : points)
      {
        (void)point_struct_b->Append();
        (void)point_x_b->Append(pt[0]); // RT
        (void)point_y_b->Append(pt[1]); // MZ
      }
    }

    // === metavalues (not nullable) ===
    (void)metavalues_builder.Append();
    appendMetaValues_(feature, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, excluded_mvs);

    // Recurse into subordinates
    int64_t this_id = static_cast<int64_t>(feature.getUniqueId());
    for (const auto& sub : feature.getSubordinates())
    {
      appendFeature(sub, this_id, depth + 1);
    }
  };

  // Walk all top-level features
  for (const auto& feature : feature_map)
  {
    appendFeature(feature, 0, 0);
  }

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_unique_id, arr_parent_id, arr_depth;
  std::shared_ptr<arrow::Array> arr_rt, arr_mz, arr_intensity;
  std::shared_ptr<arrow::Array> arr_charge, arr_overall_quality;
  std::shared_ptr<arrow::Array> arr_quality_rt, arr_quality_mz;
  std::shared_ptr<arrow::Array> arr_width;
  std::shared_ptr<arrow::Array> arr_rt_bb_min, arr_rt_bb_max, arr_mz_bb_min, arr_mz_bb_max;
  std::shared_ptr<arrow::Array> arr_convex_hulls, arr_metavalues;

  #define FINISH_OR_RETURN(builder_name, arr_name) \
    status = builder_name.Finish(&arr_name); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "FeatureMapArrowIO: " #builder_name " Finish failed: " << status.ToString() << std::endl; return nullptr; }

  FINISH_OR_RETURN(unique_id_builder, arr_unique_id)
  FINISH_OR_RETURN(parent_feature_id_builder, arr_parent_id)
  FINISH_OR_RETURN(depth_builder, arr_depth)
  FINISH_OR_RETURN(rt_builder, arr_rt)
  FINISH_OR_RETURN(mz_builder, arr_mz)
  FINISH_OR_RETURN(intensity_builder, arr_intensity)
  FINISH_OR_RETURN(charge_builder, arr_charge)
  FINISH_OR_RETURN(overall_quality_builder, arr_overall_quality)
  FINISH_OR_RETURN(quality_rt_builder, arr_quality_rt)
  FINISH_OR_RETURN(quality_mz_builder, arr_quality_mz)
  FINISH_OR_RETURN(width_builder, arr_width)
  FINISH_OR_RETURN(rt_bb_min_builder, arr_rt_bb_min)
  FINISH_OR_RETURN(rt_bb_max_builder, arr_rt_bb_max)
  FINISH_OR_RETURN(mz_bb_min_builder, arr_mz_bb_min)
  FINISH_OR_RETURN(mz_bb_max_builder, arr_mz_bb_max)
  FINISH_OR_RETURN(convex_hulls_builder, arr_convex_hulls)
  FINISH_OR_RETURN(metavalues_builder, arr_metavalues)

  #undef FINISH_OR_RETURN

  // Build schema from registry (17 columns)
  auto schema = FeatureSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_unique_id, arr_parent_id, arr_depth,
    arr_rt, arr_mz, arr_intensity,
    arr_charge, arr_overall_quality,
    arr_quality_rt, arr_quality_mz,
    arr_width,
    arr_rt_bb_min, arr_rt_bb_max, arr_mz_bb_min, arr_mz_bb_max,
    arr_convex_hulls, arr_metavalues
  });

  // Validate table against registry schema
  auto validation = ArrowSchemaValidation::validate(table, FeatureSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}

std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportPSMsToArrow(
  const FeatureMap& feature_map)
{
  // 1. Collect all PeptideIdentifications with their associated feature_ids.
  //    For each PeptideIdentification we store (feature_unique_id, is_null).
  PeptideIdentificationList all_pep_ids;
  std::vector<std::pair<int64_t, bool>> feature_ids_per_pep_id; // (feature_id, is_null)

  // Helper: recursively collect PeptideIdentifications from a Feature and its subordinates.
  std::function<void(const Feature&)> collectFromFeature = [&](const Feature& f)
  {
    for (const auto& pep_id : f.getPeptideIdentifications())
    {
      all_pep_ids.push_back(pep_id);
      feature_ids_per_pep_id.push_back({static_cast<int64_t>(f.getUniqueId()), false});
    }
    for (const auto& sub : f.getSubordinates())
    {
      collectFromFeature(sub);
    }
  };

  for (const auto& feature : feature_map)
  {
    collectFromFeature(feature);
  }

  // Unassigned peptide identifications get feature_id = null
  for (const auto& pep_id : feature_map.getUnassignedPeptideIdentifications())
  {
    all_pep_ids.push_back(pep_id);
    feature_ids_per_pep_id.push_back({0, true}); // null
  }

  // 2. Call QPXFile::exportToArrow to produce the base PSM table.
  //    We use export_all_psms=true so all hits per PeptideIdentification are exported.
  auto base_table = QPXFile::exportToArrow(
    feature_map.getProteinIdentifications(), all_pep_ids, true);
  if (!base_table) { return nullptr; }

  // 3. Build the feature_id column.
  //    QPXFile::exportToArrow (with export_all_psms=true) produces one row per
  //    PeptideHit across all PeptideIdentifications. PeptideIdentifications
  //    with empty hits are skipped. We emit one feature_id per PeptideHit.
  arrow::Int64Builder feature_id_builder;

  for (size_t i = 0; i < all_pep_ids.size(); ++i)
  {
    // Skip PeptideIdentifications with no hits (QPXFile skips them too)
    if (all_pep_ids[i].getHits().empty())
    {
      continue;
    }

    // Emit one feature_id per PeptideHit in this PeptideIdentification
    for (size_t j = 0; j < all_pep_ids[i].getHits().size(); ++j)
    {
      if (feature_ids_per_pep_id[i].second) // is_null
      {
        (void)feature_id_builder.AppendNull();
      }
      else
      {
        (void)feature_id_builder.Append(feature_ids_per_pep_id[i].first);
      }
    }
  }

  std::shared_ptr<arrow::Array> feature_id_array;
  auto status = feature_id_builder.Finish(&feature_id_array);
  if (!status.ok())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: feature_id_builder Finish failed: " << status.ToString() << std::endl;
    return nullptr;
  }

  // Verify row count matches
  if (feature_id_array->length() != base_table->num_rows())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: feature_id column length (" << feature_id_array->length()
                     << ") does not match table row count (" << base_table->num_rows() << ")" << std::endl;
    return nullptr;
  }

  // 4. Add feature_id as the first column in the table.
  auto chunked_feature_id = std::make_shared<arrow::ChunkedArray>(feature_id_array);
  auto result = base_table->AddColumn(0, arrow::field("feature_unique_id", arrow::int64()), chunked_feature_id);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: AddColumn failed: " << result.status().ToString() << std::endl;
    return nullptr;
  }
  return *result;
}

bool FeatureMapArrowIO::exportToParquet(
  const FeatureMap& feature_map,
  const String& directory,
  const ParquetWriteConfig& config)
{
  // 1. Create output directory
  try
  {
    std::filesystem::create_directories(std::string(directory));
  }
  catch (const std::filesystem::filesystem_error& e)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to create directory: " << e.what() << std::endl;
    return false;
  }

  // 2. Export features table (with FeatureMap-level metadata)
  auto features_table = exportFeaturesToArrow(feature_map);
  if (!features_table)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to create features Arrow table" << std::endl;
    return false;
  }

  // Collect FeatureMap-level metadata (DocumentIdentifier + DataProcessing)
  std::unordered_map<std::string, std::string> feature_map_metadata;
  feature_map_metadata["document_id"] = feature_map.getIdentifier();
  feature_map_metadata["loaded_file_path"] = feature_map.getLoadedFilePath();
  feature_map_metadata["loaded_file_type"] = FileTypes::typeToName(feature_map.getLoadedFileType());
  feature_map_metadata["data_processing"] = serializeDataProcessing_(feature_map.getDataProcessing());

  if (!writeArrowTableToParquet_(features_table, directory + "/features.parquet", "features", config, feature_map_metadata))
  {
    return false;
  }

  // 3. Export PSMs table
  auto psms_table = exportPSMsToArrow(feature_map);
  if (!psms_table)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to create PSMs Arrow table" << std::endl;
    return false;
  }
  if (!writeArrowTableToParquet_(psms_table, directory + "/psms.parquet", "psms", config))
  {
    return false;
  }

  // 4. Delegate protein data to ProteinIdentificationArrowIO
  const auto& prot_ids = feature_map.getProteinIdentifications();
  if (!ProteinIdentificationArrowIO::exportProteinsToParquet(
          prot_ids, directory + "/proteins.parquet", config))
  {
    return false;
  }
  if (!ProteinIdentificationArrowIO::exportProteinGroupsToParquet(
          prot_ids, directory + "/protein_groups.parquet", config))
  {
    return false;
  }
  if (!ProteinIdentificationArrowIO::exportSearchParamsToParquet(
          prot_ids, directory + "/search_params.parquet", config))
  {
    return false;
  }

  return true;
}

bool FeatureMapArrowIO::importFeaturesFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  FeatureMap& feature_map)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: null table passed to importFeaturesFromArrow" << std::endl;
    return false;
  }

  // Combine chunks for reliable single-chunk access
  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to combine chunks" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;

  int64_t num_rows = tbl->num_rows();
  if (num_rows == 0)
  {
    return true;
  }

  // Validate table schema against registry (subset mode — file may have extra columns)
  auto validation = ArrowSchemaValidation::validate(tbl, FeatureSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible schema: " << validation.toString() << "\n";
    return false;
  }

  // Get all columns
  auto col_unique_id = getColumn_(tbl, FeatureSchema::UNIQUE_ID);
  auto col_parent_id = getColumn_(tbl, FeatureSchema::PARENT_FEATURE_ID);
  auto col_depth = getColumn_(tbl, FeatureSchema::DEPTH);
  auto col_rt = getColumn_(tbl, FeatureSchema::RT);
  auto col_mz = getColumn_(tbl, FeatureSchema::MZ);
  auto col_intensity = getColumn_(tbl, FeatureSchema::INTENSITY);
  auto col_charge = getColumn_(tbl, FeatureSchema::CHARGE);
  auto col_overall_quality = getColumn_(tbl, FeatureSchema::QUALITY);
  auto col_quality_rt = getColumn_(tbl, FeatureSchema::QUALITY_RT);
  auto col_quality_mz = getColumn_(tbl, FeatureSchema::QUALITY_MZ);
  auto col_width = getColumn_(tbl, FeatureSchema::WIDTH, /*required=*/false);
  auto col_convex_hulls = getColumn_(tbl, FeatureSchema::CONVEX_HULLS, /*required=*/false);
  auto col_metavalues = getColumn_(tbl, FeatureSchema::METAVALUES, /*required=*/false);

  if (!col_unique_id || !col_parent_id || !col_depth || !col_rt || !col_mz ||
      !col_intensity || !col_charge || !col_overall_quality ||
      !col_quality_rt || !col_quality_mz)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Missing one or more required columns" << std::endl;
    return false;
  }

  // Build a vector of {depth, row_index, unique_id, parent_id, Feature} entries
  struct FeatureEntry
  {
    int32_t depth;
    int64_t row_index;
    int64_t unique_id;
    int64_t parent_id;  // -1 if top-level
    Feature feature;
  };

  std::vector<FeatureEntry> entries;
  entries.reserve(num_rows);

  // No metavalue keys to exclude for features
  static const std::unordered_set<std::string> excluded_mvs = {};

  for (int64_t i = 0; i < num_rows; ++i)
  {
    FeatureEntry entry;
    entry.row_index = i;
    entry.depth = getInt32Value_(col_depth, i, 0);
    entry.unique_id = getInt64Value_(col_unique_id, i, 0);
    entry.parent_id = isNull_(col_parent_id, i) ? -1 : getInt64Value_(col_parent_id, i, 0);

    Feature& f = entry.feature;
    f.setUniqueId(static_cast<UInt64>(entry.unique_id));
    f.setRT(getDoubleValue_(col_rt, i));
    f.setMZ(getDoubleValue_(col_mz, i));
    f.setIntensity(getFloatValue_(col_intensity, i));
    f.setCharge(static_cast<int>(getInt32Value_(col_charge, i)));
    f.setOverallQuality(getFloatValue_(col_overall_quality, i));
    f.setQuality(0, getFloatValue_(col_quality_rt, i));
    f.setQuality(1, getFloatValue_(col_quality_mz, i));

    // width: null means unset (default 0.0); only set if non-null and non-zero
    float w = getFloatValue_(col_width, i, 0.0f);
    if (w != 0.0f)
    {
      f.setWidth(w);
    }

    // Read convex hulls
    if (col_convex_hulls)
    {
      readConvexHulls_(col_convex_hulls, i, f);
    }

    // Read metavalues
    if (col_metavalues)
    {
      readMetaValues_(col_metavalues, i, f, excluded_mvs);
    }

    entries.push_back(std::move(entry));
  }

  // Sort by depth descending so children are processed before parents.
  // This allows us to build the tree bottom-up: attach children to parents
  // in the entries vector first, then add fully assembled top-level features
  // to the FeatureMap.
  std::stable_sort(entries.begin(), entries.end(),
    [](const FeatureEntry& a, const FeatureEntry& b) { return a.depth > b.depth; });

  // Build unique_id -> index-in-entries lookup map
  std::unordered_map<int64_t, size_t> id_to_index;
  for (size_t i = 0; i < entries.size(); ++i)
  {
    id_to_index[entries[i].unique_id] = i;
  }

  // Bottom-up assembly: for each entry (deepest first), attach it to its parent
  for (auto& entry : entries)
  {
    if (entry.parent_id != -1)
    {
      auto it = id_to_index.find(entry.parent_id);
      if (it != id_to_index.end())
      {
        entries[it->second].feature.getSubordinates().push_back(std::move(entry.feature));
      }
      else
      {
        OPENMS_LOG_WARN << "FeatureMapArrowIO: Could not find parent feature with id "
                        << entry.parent_id << " for feature " << entry.unique_id
                        << ". Will be added as top-level." << std::endl;
        // Mark as top-level by clearing parent_id
        entry.parent_id = -1;
      }
    }
  }

  // Add top-level features to the FeatureMap (preserving original order by row_index)
  // First, collect top-level entries sorted by row_index
  std::vector<size_t> top_level_indices;
  for (size_t i = 0; i < entries.size(); ++i)
  {
    if (entries[i].parent_id == -1)
    {
      top_level_indices.push_back(i);
    }
  }
  std::sort(top_level_indices.begin(), top_level_indices.end(),
    [&](size_t a, size_t b) { return entries[a].row_index < entries[b].row_index; });

  for (size_t idx : top_level_indices)
  {
    feature_map.push_back(std::move(entries[idx].feature));
  }

  return true;
}

bool FeatureMapArrowIO::importPSMsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  FeatureMap& feature_map)
{
  if (!table || table->num_rows() == 0) { return true; }

  // Combine chunks for reliable single-chunk access
  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to combine chunks in importPSMsFromArrow" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;
  int64_t num_rows = tbl->num_rows();

  auto psm_validation = ArrowSchemaValidation::validate(tbl, PSMSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!psm_validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible PSM schema: " << psm_validation.toString() << "\n";
    return false;
  }

  // Build feature lookup: unique_id -> Feature* (recursively includes subordinates)
  std::unordered_map<int64_t, Feature*> feature_lookup;
  std::function<void(Feature&)> buildLookup = [&](Feature& f)
  {
    feature_lookup[static_cast<int64_t>(f.getUniqueId())] = &f;
    for (auto& sub : f.getSubordinates())
    {
      buildLookup(sub);
    }
  };
  for (auto& f : feature_map)
  {
    buildLookup(f);
  }

  // Read columns
  auto col_feature_id = getColumn_(tbl, "feature_unique_id");
  auto col_p_id = getColumn_(tbl, PSMSchema::PEPTIDE_IDENTIFICATION_INDEX);
  auto col_peptidoform = getColumn_(tbl, PSMSchema::PEPTIDOFORM, /*required=*/false);
  auto col_sequence = getColumn_(tbl, PSMSchema::SEQUENCE, /*required=*/false);
  auto col_charge = getColumn_(tbl, PSMSchema::PRECURSOR_CHARGE);
  auto col_score = getColumn_(tbl, PSMSchema::SCORE);
  auto col_score_type = getColumn_(tbl, PSMSchema::SCORE_TYPE);
  auto col_rank = getColumn_(tbl, PSMSchema::RANK, /*required=*/false);
  auto col_rt = getColumn_(tbl, PSMSchema::RT, /*required=*/false);
  auto col_mz = getColumn_(tbl, PSMSchema::OBSERVED_MZ, /*required=*/false);
  auto col_spec_ref = getColumn_(tbl, PSMSchema::SPECTRUM_REFERENCE, /*required=*/false);
  auto col_run_id = getColumn_(tbl, PSMSchema::RUN_IDENTIFIER, /*required=*/false);
  auto col_is_decoy = getColumn_(tbl, PSMSchema::IS_DECOY, /*required=*/false);
  auto col_protein_accs = getColumn_(tbl, PSMSchema::PROTEIN_ACCESSIONS, /*required=*/false);
  auto col_additional_scores = getColumn_(tbl, PSMSchema::ADDITIONAL_SCORES, /*required=*/false);
  auto col_psm_metavalues = getColumn_(tbl, PSMSchema::PSM_METAVALUES, /*required=*/false);
  auto col_spectrum_metavalues = getColumn_(tbl, PSMSchema::SPECTRUM_METAVALUES, /*required=*/false);
  auto col_predicted_rt = getColumn_(tbl, PSMSchema::PREDICTED_RT, /*required=*/false);
  auto col_ion_mobility = getColumn_(tbl, PSMSchema::ION_MOBILITY, /*required=*/false);
  auto col_hsb = getColumn_(tbl, PSMSchema::HIGHER_SCORE_BETTER, /*required=*/false);
  auto col_scan = getColumn_(tbl, PSMSchema::SCAN, /*required=*/false);
  auto col_ref_file = getColumn_(tbl, PSMSchema::REFERENCE_FILE_NAME, /*required=*/false);

  if (!col_feature_id || !col_p_id || !col_charge || !col_score || !col_score_type)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Missing required columns for PSM import" << std::endl;
    return false;
  }

  // Group rows by P_ID to reconstruct PeptideIdentifications.
  // Each unique P_ID value corresponds to one PeptideIdentification.
  // Rows with the same P_ID share the same PeptideIdentification-level data
  // (rt, mz, score_type, spectrum_reference, run_identifier).
  //
  // We also track the feature_id for each group so we know where to attach it.
  // All rows in the same P_ID group should have the same feature_id.

  // Build a lookup for higher_score_better from ProteinIdentifications
  std::unordered_map<std::string, bool> higher_score_better_lookup;
  for (const auto& prot_id : feature_map.getProteinIdentifications())
  {
    higher_score_better_lookup[prot_id.getIdentifier()] = prot_id.isHigherScoreBetter();
  }

  // Process rows in order, grouping consecutive rows with the same P_ID.
  struct PepIdGroup
  {
    int64_t feature_id;
    bool feature_id_is_null;
    PeptideIdentification pep_id;
  };

  std::vector<PepIdGroup> groups;
  int32_t current_p_id = -1;

  for (int64_t row = 0; row < num_rows; ++row)
  {
    int32_t p_id = getInt32Value_(col_p_id, row, -1);

    // Start a new group if P_ID changes.
    // Note: This assumes rows with the same P_ID are contiguous in the table,
    // which is guaranteed by QPXFile::exportToArrow. If the Parquet file has been
    // externally sorted/filtered, non-contiguous rows with the same P_ID will be
    // split into separate PeptideIdentification objects.
    if (groups.empty() || p_id != current_p_id)
    {
      PepIdGroup group;
      group.feature_id_is_null = isNull_(col_feature_id, row);
      group.feature_id = group.feature_id_is_null ? 0 : getInt64Value_(col_feature_id, row, 0);

      // Set PeptideIdentification-level fields
      PeptideIdentification& pid = group.pep_id;
      pid.setScoreType(getStringValue_(col_score_type, row));
      // Run identifier (links to ProteinIdentification)
      if (col_run_id && !isNull_(col_run_id, row))
      {
        pid.setIdentifier(getStringValue_(col_run_id, row));
      }

      // higher_score_better: prefer per-PSM column, fall back to ProteinIdentification lookup
      if (col_hsb && !isNull_(col_hsb, row))
      {
        pid.setHigherScoreBetter(getBoolValue_(col_hsb, row, true));
      }
      else if (col_run_id && !isNull_(col_run_id, row))
      {
        auto hsb_it = higher_score_better_lookup.find(pid.getIdentifier());
        if (hsb_it != higher_score_better_lookup.end())
        {
          pid.setHigherScoreBetter(hsb_it->second);
        }
        else
        {
          pid.setHigherScoreBetter(true);
        }
      }
      else
      {
        pid.setHigherScoreBetter(true);
      }

      // RT
      if (col_rt && !isNull_(col_rt, row))
      {
        pid.setRT(getDoubleValue_(col_rt, row));
      }

      // MZ
      if (col_mz && !isNull_(col_mz, row))
      {
        pid.setMZ(getDoubleValue_(col_mz, row));
      }

      // Spectrum reference
      if (col_spec_ref && !isNull_(col_spec_ref, row))
      {
        pid.setSpectrumReference(getStringValue_(col_spec_ref, row));
      }

      // spectrum_metavalues -> PeptideIdentification metavalues
      if (col_spectrum_metavalues)
      {
        readMetaValues_(col_spectrum_metavalues, row, pid);
      }

      groups.push_back(std::move(group));
      current_p_id = p_id;
    }

    // Add a PeptideHit to the current group
    PeptideHit hit;

    // Reconstruct sequence from peptidoform (ProForma) if available, else from sequence column
    if (col_peptidoform && !isNull_(col_peptidoform, row))
    {
      String peptidoform_str = getStringValue_(col_peptidoform, row);
      if (!peptidoform_str.empty())
      {
        try
        {
          auto pf = ProForma::parse(peptidoform_str);
          hit.setSequence(ProForma::toAASequence(pf, ProForma::ConversionPolicy::BEST_EFFORT));
        }
        catch (...)
        {
          // Fallback: try unmodified sequence
          if (col_sequence && !isNull_(col_sequence, row))
          {
            hit.setSequence(AASequence::fromString(getStringValue_(col_sequence, row)));
          }
        }
      }
    }
    else if (col_sequence && !isNull_(col_sequence, row))
    {
      hit.setSequence(AASequence::fromString(getStringValue_(col_sequence, row)));
    }

    hit.setCharge(static_cast<Int>(getInt32Value_(col_charge, row, 0)));
    hit.setScore(getDoubleValue_(col_score, row, 0.0));

    if (col_rank && !isNull_(col_rank, row))
    {
      hit.setRank(static_cast<UInt>(getInt32Value_(col_rank, row, 0)));
    }

    // is_decoy -> target_decoy metavalue
    if (col_is_decoy && !isNull_(col_is_decoy, row))
    {
      bool is_decoy = getBoolValue_(col_is_decoy, row, false);
      hit.setMetaValue("target_decoy", is_decoy ? "decoy" : "target");
    }

    // protein_accessions -> PeptideEvidence
    if (col_protein_accs && !isNull_(col_protein_accs, row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_protein_accs);
      auto values = std::static_pointer_cast<arrow::StringArray>(list_arr->values());
      int64_t start = list_arr->value_offset(row);
      int64_t end = start + list_arr->value_length(row);
      for (int64_t k = start; k < end; ++k)
      {
        PeptideEvidence ev;
        ev.setProteinAccession(values->GetString(k));
        hit.addPeptideEvidence(ev);
      }
    }

    // additional_scores -> metavalues on PeptideHit
    if (col_additional_scores && !isNull_(col_additional_scores, row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_additional_scores);
      auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->values());
      auto names_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("score_name"));
      auto values_arr = std::static_pointer_cast<arrow::DoubleArray>(struct_arr->GetFieldByName("score_value"));

      int64_t start = list_arr->value_offset(row);
      int64_t end = start + list_arr->value_length(row);
      for (int64_t k = start; k < end; ++k)
      {
        String name = names_arr->GetString(k);
        double value = values_arr->Value(k);
        hit.setMetaValue(name, value);
      }
    }

    // predicted_rt -> PeptideHit metavalue
    if (col_predicted_rt && !isNull_(col_predicted_rt, row))
    {
      hit.setMetaValue("predicted_RT", getDoubleValue_(col_predicted_rt, row));
    }

    // ion_mobility -> PeptideHit metavalue
    if (col_ion_mobility && !isNull_(col_ion_mobility, row))
    {
      hit.setMetaValue("ion_mobility", getDoubleValue_(col_ion_mobility, row));
    }

    // scan -> PeptideHit metavalue
    if (col_scan && !isNull_(col_scan, row))
    {
      hit.setMetaValue("scan", static_cast<int>(getInt32Value_(col_scan, row)));
    }

    // reference_file_name -> PeptideHit metavalue
    if (col_ref_file && !isNull_(col_ref_file, row))
    {
      hit.setMetaValue("reference_file_name", getStringValue_(col_ref_file, row));
    }

    // psm_metavalues -> PeptideHit metavalues (exclude keys that have dedicated columns)
    if (col_psm_metavalues)
    {
      static const std::unordered_set<std::string> psm_excluded_mvs =
        {"target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
         "scan", "reference_file_name"};
      readMetaValues_(col_psm_metavalues, row, hit, psm_excluded_mvs);
    }

    groups.back().pep_id.getHits().push_back(std::move(hit));
  }

  // Assign PeptideIdentifications to features or as unassigned
  for (auto& group : groups)
  {
    if (group.feature_id_is_null)
    {
      // Unassigned
      feature_map.getUnassignedPeptideIdentifications().push_back(std::move(group.pep_id));
    }
    else
    {
      auto it = feature_lookup.find(group.feature_id);
      if (it != feature_lookup.end())
      {
        it->second->getPeptideIdentifications().push_back(std::move(group.pep_id));
      }
      else
      {
        OPENMS_LOG_WARN << "FeatureMapArrowIO: Could not find feature with id "
                        << group.feature_id << " for PSM. Adding as unassigned." << std::endl;
        feature_map.getUnassignedPeptideIdentifications().push_back(std::move(group.pep_id));
      }
    }
  }

  return true;
}

bool FeatureMapArrowIO::importFromParquet(
  const String& directory,
  FeatureMap& feature_map)
{
  feature_map = FeatureMap{};

  // 1. Import protein identification data
  std::vector<ProteinIdentification> prot_ids;
  if (!ProteinIdentificationArrowIO::importFromParquet(
          directory + "/proteins.parquet",
          directory + "/protein_groups.parquet",
          directory + "/search_params.parquet",
          prot_ids))
  {
    return false;
  }
  feature_map.setProteinIdentifications(prot_ids);

  // 2. Import features and FeatureMap-level metadata
  auto features_table = readParquetTable_(directory + "/features.parquet");
  if (!features_table) { return false; }

  // Extract FeatureMap-level metadata from file-level key-value metadata
  auto schema_md = features_table->schema()->metadata();
  if (schema_md)
  {
    // DocumentIdentifier fields
    int idx = schema_md->FindKey("document_id");
    if (idx >= 0) feature_map.setIdentifier(schema_md->value(idx));

    idx = schema_md->FindKey("loaded_file_path");
    if (idx >= 0) feature_map.setLoadedFilePath(schema_md->value(idx));

    // loaded_file_type: stored as type name string (e.g. "featureXML").
    // Note: DocumentIdentifier::setLoadedFileType(String) expects a file path and
    // determines type by content/extension, so we cannot directly restore the type
    // from a type name. The type can be re-derived from loaded_file_path if needed.

    // DataProcessing
    idx = schema_md->FindKey("data_processing");
    if (idx >= 0)
    {
      feature_map.setDataProcessing(deserializeDataProcessing_(schema_md->value(idx)));
    }
  }

  if (!importFeaturesFromArrow(features_table, feature_map))
  {
    return false;
  }

  // 3. Import PSMs (links to features already in the map)
  auto psms_table = readParquetTable_(directory + "/psms.parquet");
  if (!psms_table) { return false; }
  if (!importPSMsFromArrow(psms_table, feature_map))
  {
    return false;
  }

  return true;
}

} // namespace OpenMS
