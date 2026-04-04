// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusMapArrowIO.h>

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

#include <cctype>
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
              for (int j = 0; j < 4; ++j)
              {
                char ch = s[i + 1 + j];
                cp <<= 4;
                if (ch >= '0' && ch <= '9') cp |= (ch - '0');
                else if (ch >= 'a' && ch <= 'f') cp |= (ch - 'a' + 10);
                else if (ch >= 'A' && ch <= 'F') cp |= (ch - 'A' + 10);
              }
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

      json += "\"software_name\":\"" + escapeJsonString_(dp.getSoftware().getName()) + "\"";
      json += ",\"software_version\":\"" + escapeJsonString_(dp.getSoftware().getVersion()) + "\"";

      std::string time_str;
      if (dp.getCompletionTime().isValid())
      {
        time_str = dp.getCompletionTime().toString("yyyy-MM-ddThh:mm:ss");
      }
      json += ",\"completion_time\":\"" + escapeJsonString_(time_str) + "\"";

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
  std::string parseJsonString_(const std::string& s, size_t& pos)
  {
    if (pos >= s.size() || s[pos] != '"') return "";
    ++pos;
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
    if (pos < s.size()) ++pos;
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
    ++pos;

    while (pos < json.size())
    {
      skipWhitespace_(json, pos);
      if (pos >= json.size() || json[pos] == ']') break;
      if (json[pos] == ',') { ++pos; continue; }
      if (json[pos] != '{') break;
      ++pos;

      DataProcessing dp;
      std::string software_name, software_version, completion_time;
      std::vector<std::string> action_names;

      while (pos < json.size())
      {
        skipWhitespace_(json, pos);
        if (pos >= json.size() || json[pos] == '}') { ++pos; break; }
        if (json[pos] == ',') { ++pos; continue; }

        std::string key = parseJsonString_(json, pos);
        skipWhitespace_(json, pos);
        if (pos < json.size() && json[pos] == ':') ++pos;
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
          if (pos < json.size() && json[pos] == '[')
          {
            ++pos;
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
          if (pos < json.size() && json[pos] == '[')
          {
            ++pos;
            while (pos < json.size())
            {
              skipWhitespace_(json, pos);
              if (pos >= json.size() || json[pos] == ']') { ++pos; break; }
              if (json[pos] == ',') { ++pos; continue; }
              if (json[pos] != '{') break;
              ++pos;

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
                else parseJsonString_(json, pos);
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
          parseJsonString_(json, pos);
        }
      }

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
          OPENMS_LOG_WARN << "ConsensusMapArrowIO: Unknown processing action: " << action_name << std::endl;
        }
      }

      result.push_back(std::move(dp));
    }

    return result;
  }

  // ==================== MetaInfoInterface JSON helpers ====================

  /// Serialize MetaInfoInterface metavalues to a JSON array string.
  std::string serializeMetaValues_(const MetaInfoInterface& mii)
  {
    std::string json = "[";
    std::vector<String> keys;
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
      json += "{\"name\":\"" + escapeJsonString_(std::string(key))
            + "\",\"value\":\"" + escapeJsonString_(val.toString())
            + "\",\"type\":\"" + type_str + "\"}";
      first = false;
    }
    json += "]";
    return json;
  }

  /// Deserialize a JSON array of {name, value, type} objects into a MetaInfoInterface.
  void deserializeMetaValues_(const std::string& json, MetaInfoInterface& target)
  {
    if (json.empty()) return;

    size_t pos = 0;
    skipWhitespace_(json, pos);
    if (pos >= json.size() || json[pos] != '[') return;
    ++pos;

    while (pos < json.size())
    {
      skipWhitespace_(json, pos);
      if (pos >= json.size() || json[pos] == ']') break;
      if (json[pos] == ',') { ++pos; continue; }
      if (json[pos] != '{') break;
      ++pos;

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
        else parseJsonString_(json, pos);
      }

      if (!mv_name.empty())
      {
        if (mv_type == "int")
        {
          try { target.setMetaValue(mv_name, DataValue(std::stoi(mv_value))); }
          catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
        }
        else if (mv_type == "double" || mv_type == "float")
        {
          try { target.setMetaValue(mv_name, DataValue(std::stod(mv_value))); }
          catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
        }
        else if (mv_type == "int_list")
        {
          try
          {
            String s(mv_value);
            if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
            target.setMetaValue(mv_name, DataValue(ListUtils::create<Int>(s)));
          }
          catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
        }
        else if (mv_type == "double_list")
        {
          try
          {
            String s(mv_value);
            if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
            target.setMetaValue(mv_name, DataValue(ListUtils::create<double>(s)));
          }
          catch (...) { target.setMetaValue(mv_name, DataValue(mv_value)); }
        }
        else if (mv_type == "string_list")
        {
          try
          {
            String s(mv_value);
            if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
            auto sl = ListUtils::create<String>(s);
            for (auto& e : sl) { e = e.trim(); }
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

  // ==================== Column header JSON helpers ====================

  /// Serialize ConsensusMap column headers to a JSON string (including metavalues).
  std::string serializeColumnHeaders_(const ConsensusMap::ColumnHeaders& headers)
  {
    std::string json = "[";
    bool first = true;
    for (const auto& [map_index, header] : headers)
    {
      if (!first) json += ",";
      json += "{\"map_index\":" + std::to_string(map_index)
            + ",\"filename\":\"" + escapeJsonString_(header.filename) + "\""
            + ",\"label\":\"" + escapeJsonString_(header.label) + "\""
            + ",\"size\":" + std::to_string(header.size)
            + ",\"unique_id\":" + std::to_string(header.unique_id)
            + ",\"metavalues\":" + serializeMetaValues_(header)
            + "}";
      first = false;
    }
    json += "]";
    return json;
  }

  /// Deserialize a JSON string into ConsensusMap column headers.
  ConsensusMap::ColumnHeaders deserializeColumnHeaders_(const std::string& json)
  {
    ConsensusMap::ColumnHeaders result;
    if (json.empty()) return result;

    size_t pos = 0;
    skipWhitespace_(json, pos);
    if (pos >= json.size() || json[pos] != '[') return result;
    ++pos;

    while (pos < json.size())
    {
      skipWhitespace_(json, pos);
      if (pos >= json.size() || json[pos] == ']') break;
      if (json[pos] == ',') { ++pos; continue; }
      if (json[pos] != '{') break;
      ++pos;

      ConsensusMap::ColumnHeader header;
      UInt64 map_index = 0;

      while (pos < json.size())
      {
        skipWhitespace_(json, pos);
        if (pos >= json.size() || json[pos] == '}') { ++pos; break; }
        if (json[pos] == ',') { ++pos; continue; }

        std::string key = parseJsonString_(json, pos);
        skipWhitespace_(json, pos);
        if (pos < json.size() && json[pos] == ':') ++pos;
        skipWhitespace_(json, pos);

        if (key == "map_index" || key == "size" || key == "unique_id")
        {
          std::string num_str;
          while (pos < json.size() && (std::isdigit(json[pos]) || json[pos] == '-'))
          {
            num_str += json[pos++];
          }
          if (key == "map_index") map_index = std::stoull(num_str);
          else if (key == "size") header.size = std::stoull(num_str);
          else if (key == "unique_id") header.unique_id = std::stoull(num_str);
        }
        else if (key == "metavalues")
        {
          // Parse the nested metavalues array inline
          // Find the matching ']' to extract the substring for deserializeMetaValues_
          if (pos < json.size() && json[pos] == '[')
          {
            size_t start = pos;
            int depth = 0;
            while (pos < json.size())
            {
              if (json[pos] == '[') ++depth;
              else if (json[pos] == ']') { --depth; if (depth == 0) { ++pos; break; } }
              else if (json[pos] == '"') { ++pos; while (pos < json.size() && json[pos] != '"') { if (json[pos] == '\\') ++pos; ++pos; } }
              ++pos;
            }
            deserializeMetaValues_(json.substr(start, pos - start), header);
          }
        }
        else
        {
          std::string val = parseJsonString_(json, pos);
          if (key == "filename") header.filename = val;
          else if (key == "label") header.label = val;
        }
      }
      result[map_index] = header;
    }
    return result;
  }

  // ==================== Parquet I/O helpers ====================

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
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: null table passed to writeArrowTableToParquet_" << std::endl;
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
    bytes[6] = (bytes[6] & 0x0F) | 0x40;
    bytes[8] = (bytes[8] & 0x3F) | 0x80;
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

    auto result = arrow::io::FileOutputStream::Open(filename);
    if (!result.ok())
    {
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to open file: " << filename << std::endl;
      return false;
    }
    const auto& outfile = *result;

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

    auto write_status = parquet::arrow::WriteTable(
      *table, arrow::default_memory_pool(), outfile,
      config.row_group_size, writer_props, arrow_props);

    if (!write_status.ok())
    {
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to write Parquet: "
                       << write_status.ToString() << std::endl;
      return false;
    }

    return true;
  }

  // ==================== Import helpers ====================

  /// Read a single Parquet file into an Arrow table.
  std::shared_ptr<arrow::Table> readParquetTable_(const String& filename)
  {
    auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
    if (!infile_result.ok())
    {
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to open file: " << filename << std::endl;
      return nullptr;
    }
    const auto& infile = *infile_result;

    auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    if (!reader_result.ok())
    {
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to create Parquet reader for: " << filename << std::endl;
      return nullptr;
    }
    auto reader = std::move(reader_result.ValueOrDie());

    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    if (!read_status.ok())
    {
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to read table: " << read_status.ToString() << std::endl;
      return nullptr;
    }

    auto combined = table->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok())
    {
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to combine chunks" << std::endl;
      return nullptr;
    }

    return *combined;
  }

  /// Fetch a named column from a table, combining chunks if needed.
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
        OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Missing required column '" << name << "'" << std::endl;
      }
      return nullptr;
    }
    if (column->num_chunks() == 0)
    {
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Column '" << name << "' has no chunks" << std::endl;
      return nullptr;
    }
    if (column->num_chunks() == 1)
    {
      return column->chunk(0);
    }
    auto combined = arrow::Concatenate(column->chunks(), arrow::default_memory_pool());
    if (!combined.ok())
    {
      OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to combine chunks for column '" << name << "'" << std::endl;
      return nullptr;
    }
    return *combined;
  }

  String getStringValue_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    if (!array || array->IsNull(row)) return "";
    return std::static_pointer_cast<arrow::StringArray>(array)->GetString(row);
  }

  double getDoubleValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, double default_val = 0.0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::DoubleArray>(array)->Value(row);
  }

  float getFloatValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, float default_val = 0.0f)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::FloatArray>(array)->Value(row);
  }

  int64_t getInt64Value_(const std::shared_ptr<arrow::Array>& array, int64_t row, int64_t default_val = 0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::Int64Array>(array)->Value(row);
  }

  int32_t getInt32Value_(const std::shared_ptr<arrow::Array>& array, int64_t row, int32_t default_val = 0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::Int32Array>(array)->Value(row);
  }

  bool getBoolValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, bool default_val = false)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::BooleanArray>(array)->Value(row);
  }

  bool isNull_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    return !array || array->IsNull(row);
  }

  /// Read metavalues from a list<struct{name,value,value_type}> column at a given row.
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

  /// Read FeatureHandles from a list<struct{...}> column at a given row into a ConsensusFeature.
  void readHandles_(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    ConsensusFeature& cf)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto map_index_arr = std::static_pointer_cast<arrow::Int64Array>(struct_arr->field(0));
    auto unique_id_arr = std::static_pointer_cast<arrow::Int64Array>(struct_arr->field(1));
    auto rt_arr = std::static_pointer_cast<arrow::DoubleArray>(struct_arr->field(2));
    auto mz_arr = std::static_pointer_cast<arrow::DoubleArray>(struct_arr->field(3));
    auto intensity_arr = std::static_pointer_cast<arrow::FloatArray>(struct_arr->field(4));
    auto charge_arr = std::static_pointer_cast<arrow::Int32Array>(struct_arr->field(5));
    auto width_arr = std::static_pointer_cast<arrow::FloatArray>(struct_arr->field(6));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      FeatureHandle handle;
      handle.setMapIndex(static_cast<UInt64>(map_index_arr->Value(i)));
      handle.setUniqueId(static_cast<UInt64>(unique_id_arr->Value(i)));
      handle.setRT(rt_arr->Value(i));
      handle.setMZ(mz_arr->Value(i));
      handle.setIntensity(intensity_arr->Value(i));
      handle.setCharge(static_cast<Int>(charge_arr->Value(i)));
      handle.setWidth(width_arr->Value(i));
      cf.insert(handle);
    }
  }

} // anonymous namespace


// ==================== Export methods ====================

std::shared_ptr<arrow::Table> ConsensusMapArrowIO::exportFeaturesToArrow(
  const ConsensusMap& cmap)
{
  // -- Scalar column builders --
  arrow::Int64Builder unique_id_builder;
  arrow::DoubleBuilder rt_builder, mz_builder;
  arrow::FloatBuilder intensity_builder;
  arrow::Int32Builder charge_builder;
  arrow::FloatBuilder quality_builder;
  arrow::FloatBuilder width_builder;

  // -- handles: list<struct{map_index, unique_id, rt, mz, intensity, charge, width}> --
  auto handle_map_index_b = std::make_shared<arrow::Int64Builder>();
  auto handle_unique_id_b = std::make_shared<arrow::Int64Builder>();
  auto handle_rt_b = std::make_shared<arrow::DoubleBuilder>();
  auto handle_mz_b = std::make_shared<arrow::DoubleBuilder>();
  auto handle_intensity_b = std::make_shared<arrow::FloatBuilder>();
  auto handle_charge_b = std::make_shared<arrow::Int32Builder>();
  auto handle_width_b = std::make_shared<arrow::FloatBuilder>();

  auto handle_struct_type = arrow::struct_({
    arrow::field("map_index", arrow::int64()),
    arrow::field("unique_id", arrow::int64()),
    arrow::field("rt", arrow::float64()),
    arrow::field("mz", arrow::float64()),
    arrow::field("intensity", arrow::float32()),
    arrow::field("charge", arrow::int32()),
    arrow::field("width", arrow::float32())
  });

  auto handle_struct_b = std::make_shared<arrow::StructBuilder>(
    handle_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{
      handle_map_index_b, handle_unique_id_b, handle_rt_b, handle_mz_b,
      handle_intensity_b, handle_charge_b, handle_width_b});
  arrow::ListBuilder handles_builder(arrow::default_memory_pool(), handle_struct_b);

  // -- metavalues: list<struct{name, value, value_type}> --
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

  Size num_rows = cmap.size();

  arrow::Status status;

  #define RESERVE_OR_RETURN(builder_name) \
    status = builder_name.Reserve(num_rows); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowIO: " #builder_name " Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  RESERVE_OR_RETURN(unique_id_builder)
  RESERVE_OR_RETURN(rt_builder)
  RESERVE_OR_RETURN(mz_builder)
  RESERVE_OR_RETURN(intensity_builder)
  RESERVE_OR_RETURN(charge_builder)
  RESERVE_OR_RETURN(quality_builder)
  RESERVE_OR_RETURN(width_builder)

  #undef RESERVE_OR_RETURN

  // Exclude FWHM from metavalues (it is already stored as the dedicated 'width' column)
  static const std::unordered_set<std::string> excluded_mvs = {"FWHM"};

  for (const auto& cf : cmap)
  {
    (void)unique_id_builder.Append(static_cast<int64_t>(cf.getUniqueId()));
    (void)rt_builder.Append(cf.getRT());
    (void)mz_builder.Append(cf.getMZ());
    (void)intensity_builder.Append(cf.getIntensity());
    (void)charge_builder.Append(static_cast<int32_t>(cf.getCharge()));
    (void)quality_builder.Append(cf.getQuality());

    // Width: 0.0 is the unset default, stored as null to distinguish from explicitly set values.
    // On import, null is left as default (0.0), so round-trip is numerically correct.
    float w = cf.getWidth();
    if (w == 0.0f)
    {
      (void)width_builder.AppendNull();
    }
    else
    {
      (void)width_builder.Append(w);
    }

    // handles
    (void)handles_builder.Append();
    for (const auto& handle : cf.getFeatures())
    {
      (void)handle_struct_b->Append();
      (void)handle_map_index_b->Append(static_cast<int64_t>(handle.getMapIndex()));
      (void)handle_unique_id_b->Append(static_cast<int64_t>(handle.getUniqueId()));
      (void)handle_rt_b->Append(handle.getRT());
      (void)handle_mz_b->Append(handle.getMZ());
      (void)handle_intensity_b->Append(handle.getIntensity());
      (void)handle_charge_b->Append(static_cast<int32_t>(handle.getCharge()));
      (void)handle_width_b->Append(handle.getWidth());
    }

    // metavalues
    (void)metavalues_builder.Append();
    appendMetaValues_(cf, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, excluded_mvs);
  }

  // Finalize arrays
  std::shared_ptr<arrow::Array> arr_unique_id, arr_rt, arr_mz, arr_intensity;
  std::shared_ptr<arrow::Array> arr_charge, arr_quality, arr_width;
  std::shared_ptr<arrow::Array> arr_handles, arr_metavalues;

  #define FINISH_OR_RETURN(builder_name, arr_name) \
    status = builder_name.Finish(&arr_name); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowIO: " #builder_name " Finish failed: " << status.ToString() << std::endl; return nullptr; }

  FINISH_OR_RETURN(unique_id_builder, arr_unique_id)
  FINISH_OR_RETURN(rt_builder, arr_rt)
  FINISH_OR_RETURN(mz_builder, arr_mz)
  FINISH_OR_RETURN(intensity_builder, arr_intensity)
  FINISH_OR_RETURN(charge_builder, arr_charge)
  FINISH_OR_RETURN(quality_builder, arr_quality)
  FINISH_OR_RETURN(width_builder, arr_width)
  FINISH_OR_RETURN(handles_builder, arr_handles)
  FINISH_OR_RETURN(metavalues_builder, arr_metavalues)

  #undef FINISH_OR_RETURN

  // Build schema from registry
  auto schema = ConsensusFeatureSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_unique_id, arr_rt, arr_mz, arr_intensity,
    arr_charge, arr_quality, arr_width,
    arr_handles, arr_metavalues
  });

  // Validate table against registry schema
  auto validation = ArrowSchemaValidation::validate(table, ConsensusFeatureSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}

std::shared_ptr<arrow::Table> ConsensusMapArrowIO::exportPSMsToArrow(
  const ConsensusMap& cmap)
{
  // Collect all PeptideIdentifications with their associated consensus_feature_unique_ids.
  PeptideIdentificationList all_pep_ids;
  std::vector<std::pair<int64_t, bool>> feature_ids_per_pep_id; // (id, is_null)

  // Collect from consensus features (flat - no subordinates)
  for (const auto& cf : cmap)
  {
    for (const auto& pep_id : cf.getPeptideIdentifications())
    {
      all_pep_ids.push_back(pep_id);
      feature_ids_per_pep_id.push_back({static_cast<int64_t>(cf.getUniqueId()), false});
    }
  }

  // Unassigned peptide identifications get consensus_feature_unique_id = null
  for (const auto& pep_id : cmap.getUnassignedPeptideIdentifications())
  {
    all_pep_ids.push_back(pep_id);
    feature_ids_per_pep_id.push_back({0, true});
  }

  // Delegate to QPXFile to produce the base PSM table
  auto base_table = QPXFile::exportToArrow(
    cmap.getProteinIdentifications(), all_pep_ids, true);
  if (!base_table) { return nullptr; }

  // Build the consensus_feature_unique_id column (one value per PeptideHit)
  arrow::Int64Builder feature_id_builder;

  for (size_t i = 0; i < all_pep_ids.size(); ++i)
  {
    if (all_pep_ids[i].getHits().empty())
    {
      continue;
    }

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
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: feature_id_builder Finish failed: " << status.ToString() << std::endl;
    return nullptr;
  }

  if (feature_id_array->length() != base_table->num_rows())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: consensus_feature_unique_id column length (" << feature_id_array->length()
                     << ") does not match table row count (" << base_table->num_rows() << ")" << std::endl;
    return nullptr;
  }

  auto chunked_feature_id = std::make_shared<arrow::ChunkedArray>(feature_id_array);
  auto result = base_table->AddColumn(0, arrow::field("consensus_feature_unique_id", arrow::int64()), chunked_feature_id);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: AddColumn failed: " << result.status().ToString() << std::endl;
    return nullptr;
  }
  return *result;
}

bool ConsensusMapArrowIO::exportToParquet(
  const ConsensusMap& cmap,
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
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to create directory: " << e.what() << std::endl;
    return false;
  }

  // 2. Export consensus features table (with ConsensusMap-level metadata)
  auto features_table = exportFeaturesToArrow(cmap);
  if (!features_table)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to create consensus features Arrow table" << std::endl;
    return false;
  }

  std::unordered_map<std::string, std::string> cmap_metadata;
  cmap_metadata["column_headers"] = serializeColumnHeaders_(cmap.getColumnHeaders());
  cmap_metadata["experiment_type"] = cmap.getExperimentType();
  cmap_metadata["document_id"] = cmap.getIdentifier();
  cmap_metadata["loaded_file_path"] = cmap.getLoadedFilePath();
  cmap_metadata["loaded_file_type"] = FileTypes::typeToName(cmap.getLoadedFileType());
  cmap_metadata["data_processing"] = serializeDataProcessing_(cmap.getDataProcessing());
  cmap_metadata["unique_id"] = std::to_string(cmap.getUniqueId());
  cmap_metadata["cmap_metavalues"] = serializeMetaValues_(cmap);

  if (!writeArrowTableToParquet_(features_table, directory + "/consensus_features.parquet",
                                  "consensus_features", config, cmap_metadata))
  {
    return false;
  }

  // 3. Export PSMs table
  auto psms_table = exportPSMsToArrow(cmap);
  if (!psms_table)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to create PSMs Arrow table" << std::endl;
    return false;
  }
  if (!writeArrowTableToParquet_(psms_table, directory + "/psms.parquet", "psms", config))
  {
    return false;
  }

  // 4. Delegate protein data to ProteinIdentificationArrowIO
  const auto& prot_ids = cmap.getProteinIdentifications();
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

// ==================== Import methods ====================

bool ConsensusMapArrowIO::importFeaturesFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  ConsensusMap& cmap)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: null table passed to importFeaturesFromArrow" << std::endl;
    return false;
  }

  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to combine chunks" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;

  int64_t num_rows = tbl->num_rows();
  if (num_rows == 0)
  {
    return true;
  }

  // Validate table schema against registry (subset mode — file may have extra columns)
  auto validation = ArrowSchemaValidation::validate(tbl, ConsensusFeatureSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible schema: " << validation.toString() << "\n";
    return false;
  }

  auto col_unique_id = getColumn_(tbl, ConsensusFeatureSchema::UNIQUE_ID);
  auto col_rt = getColumn_(tbl, ConsensusFeatureSchema::RT);
  auto col_mz = getColumn_(tbl, ConsensusFeatureSchema::MZ);
  auto col_intensity = getColumn_(tbl, ConsensusFeatureSchema::INTENSITY);
  auto col_charge = getColumn_(tbl, ConsensusFeatureSchema::CHARGE);
  auto col_quality = getColumn_(tbl, ConsensusFeatureSchema::QUALITY);
  auto col_width = getColumn_(tbl, ConsensusFeatureSchema::WIDTH, /*required=*/false);
  auto col_handles = getColumn_(tbl, ConsensusFeatureSchema::HANDLES, /*required=*/false);
  auto col_metavalues = getColumn_(tbl, ConsensusFeatureSchema::METAVALUES, /*required=*/false);

  if (!col_unique_id || !col_rt || !col_mz || !col_intensity || !col_charge || !col_quality)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Missing one or more required columns" << std::endl;
    return false;
  }

  for (int64_t i = 0; i < num_rows; ++i)
  {
    ConsensusFeature cf;
    cf.setUniqueId(static_cast<UInt64>(getInt64Value_(col_unique_id, i, 0)));
    cf.setRT(getDoubleValue_(col_rt, i));
    cf.setMZ(getDoubleValue_(col_mz, i));
    cf.setIntensity(getFloatValue_(col_intensity, i));
    cf.setCharge(static_cast<Int>(getInt32Value_(col_charge, i)));
    cf.setQuality(getFloatValue_(col_quality, i));

    // Width: null means unset (default 0.0), so we only call setWidth for non-null values.
    if (col_width && !isNull_(col_width, i))
    {
      cf.setWidth(getFloatValue_(col_width, i));
    }

    if (col_handles)
    {
      readHandles_(col_handles, i, cf);
    }

    if (col_metavalues)
    {
      readMetaValues_(col_metavalues, i, cf);
    }

    cmap.push_back(std::move(cf));
  }

  return true;
}

bool ConsensusMapArrowIO::importPSMsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  ConsensusMap& cmap)
{
  if (!table || table->num_rows() == 0) { return true; }

  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to combine chunks in importPSMsFromArrow" << std::endl;
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

  // Build consensus feature lookup: unique_id -> ConsensusFeature*
  std::unordered_map<int64_t, ConsensusFeature*> feature_lookup;
  for (auto& cf : cmap)
  {
    feature_lookup[static_cast<int64_t>(cf.getUniqueId())] = &cf;
  }

  // Read columns
  auto col_feature_id = getColumn_(tbl, "consensus_feature_unique_id");
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
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Missing required columns for PSM import" << std::endl;
    return false;
  }

  // Build a lookup for higher_score_better from ProteinIdentifications
  std::unordered_map<std::string, bool> higher_score_better_lookup;
  for (const auto& prot_id : cmap.getProteinIdentifications())
  {
    higher_score_better_lookup[prot_id.getIdentifier()] = prot_id.isHigherScoreBetter();
  }

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

    if (groups.empty() || p_id != current_p_id)
    {
      PepIdGroup group;
      group.feature_id_is_null = isNull_(col_feature_id, row);
      group.feature_id = group.feature_id_is_null ? 0 : getInt64Value_(col_feature_id, row, 0);

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

      if (col_rt && !isNull_(col_rt, row))
      {
        pid.setRT(getDoubleValue_(col_rt, row));
      }

      if (col_mz && !isNull_(col_mz, row))
      {
        pid.setMZ(getDoubleValue_(col_mz, row));
      }

      if (col_spec_ref && !isNull_(col_spec_ref, row))
      {
        pid.setSpectrumReference(getStringValue_(col_spec_ref, row));
      }

      if (col_spectrum_metavalues)
      {
        readMetaValues_(col_spectrum_metavalues, row, pid);
      }

      groups.push_back(std::move(group));
      current_p_id = p_id;
    }

    // Add a PeptideHit to the current group
    PeptideHit hit;

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

    if (col_is_decoy && !isNull_(col_is_decoy, row))
    {
      bool is_decoy = getBoolValue_(col_is_decoy, row, false);
      hit.setMetaValue("target_decoy", is_decoy ? "decoy" : "target");
    }

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

    if (col_predicted_rt && !isNull_(col_predicted_rt, row))
    {
      hit.setMetaValue("predicted_RT", getDoubleValue_(col_predicted_rt, row));
    }

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

    if (col_psm_metavalues)
    {
      static const std::unordered_set<std::string> psm_excluded_mvs =
        {"target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
         "scan", "reference_file_name"};
      readMetaValues_(col_psm_metavalues, row, hit, psm_excluded_mvs);
    }

    groups.back().pep_id.getHits().push_back(std::move(hit));
  }

  // Assign PeptideIdentifications to consensus features or as unassigned
  for (auto& group : groups)
  {
    if (group.feature_id_is_null)
    {
      cmap.getUnassignedPeptideIdentifications().push_back(std::move(group.pep_id));
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
        OPENMS_LOG_WARN << "ConsensusMapArrowIO: Could not find consensus feature with id "
                        << group.feature_id << " for PSM. Adding as unassigned." << std::endl;
        cmap.getUnassignedPeptideIdentifications().push_back(std::move(group.pep_id));
      }
    }
  }

  return true;
}

bool ConsensusMapArrowIO::importFromParquet(
  const String& directory,
  ConsensusMap& cmap)
{
  cmap = ConsensusMap{};

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
  cmap.setProteinIdentifications(prot_ids);

  // 2. Import consensus features and ConsensusMap-level metadata
  auto features_table = readParquetTable_(directory + "/consensus_features.parquet");
  if (!features_table) { return false; }

  // Extract ConsensusMap-level metadata from file-level key-value metadata
  auto schema_md = features_table->schema()->metadata();
  if (schema_md)
  {
    int idx = schema_md->FindKey("column_headers");
    if (idx >= 0) cmap.setColumnHeaders(deserializeColumnHeaders_(schema_md->value(idx)));

    idx = schema_md->FindKey("experiment_type");
    if (idx >= 0) cmap.setExperimentType(schema_md->value(idx));

    idx = schema_md->FindKey("document_id");
    if (idx >= 0) cmap.setIdentifier(schema_md->value(idx));

    idx = schema_md->FindKey("loaded_file_path");
    if (idx >= 0) cmap.setLoadedFilePath(schema_md->value(idx));

    idx = schema_md->FindKey("data_processing");
    if (idx >= 0)
    {
      cmap.setDataProcessing(deserializeDataProcessing_(schema_md->value(idx)));
    }

    idx = schema_md->FindKey("unique_id");
    if (idx >= 0)
    {
      try { cmap.setUniqueId(std::stoull(schema_md->value(idx))); }
      catch (...) {}
    }

    idx = schema_md->FindKey("cmap_metavalues");
    if (idx >= 0)
    {
      deserializeMetaValues_(schema_md->value(idx), cmap);
    }
  }

  if (!importFeaturesFromArrow(features_table, cmap))
  {
    return false;
  }

  // 3. Import PSMs (links to consensus features already in the map)
  auto psms_table = readParquetTable_(directory + "/psms.parquet");
  if (!psms_table) { return false; }
  if (!importPSMsFromArrow(psms_table, cmap))
  {
    return false;
  }

  return true;
}

} // namespace OpenMS
