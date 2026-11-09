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

#include <charconv>
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

// ==================== JSON writing ====================

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

// ==================== JSON reading ====================
//
// The metadata JSON of these files is small, so a strict parser is not worth its failure mode:
// nlohmann::json rejects the whole document for a single byte that is not valid UTF-8, which
// ParquetFile::jsonEscape() passes through unchanged. This reader is lenient instead: values of
// an unexpected JSON type are skipped (and counted), unknown keys are ignored, and malformed
// input stops the parse instead of hanging. Callers log one warning per document.

namespace
{
  void appendUtf8(std::string& out, unsigned int cp)
  {
    if (cp < 0x80)
    {
      out += static_cast<char>(cp);
    }
    else if (cp < 0x800)
    {
      out += static_cast<char>(0xC0 | (cp >> 6));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    }
    else if (cp < 0x10000)
    {
      out += static_cast<char>(0xE0 | (cp >> 12));
      out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    }
    else
    {
      out += static_cast<char>(0xF0 | (cp >> 18));
      out += static_cast<char>(0x80 | ((cp >> 12) & 0x3F));
      out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    }
  }

  /// The 4 hex digits at @p s[i..i+3] as a code unit, or false if they are not 4 hex digits.
  bool parseHex4(const std::string& s, size_t i, unsigned int& cu)
  {
    if (i + 4 > s.size()) return false;
    cu = 0;
    for (size_t j = i; j < i + 4; ++j)
    {
      const char ch = s[j];
      cu <<= 4;
      if (ch >= '0' && ch <= '9') cu |= static_cast<unsigned int>(ch - '0');
      else if (ch >= 'a' && ch <= 'f') cu |= static_cast<unsigned int>(ch - 'a' + 10);
      else if (ch >= 'A' && ch <= 'F') cu |= static_cast<unsigned int>(ch - 'A' + 10);
      else return false;
    }
    return true;
  }

  /// Resolve the escapes of a JSON string body. \uXXXX surrogate pairs become one code point; a
  /// lone surrogate becomes U+FFFD; an escape with invalid hex digits is kept as written.
  std::string unescapeJsonString(const std::string& s)
  {
    std::string result;
    result.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i)
    {
      if (s[i] != '\\' || i + 1 >= s.size())
      {
        result += s[i];
        continue;
      }
      const char c = s[++i];
      switch (c)
      {
        case '"':  result += '"'; break;
        case '\\': result += '\\'; break;
        case '/':  result += '/'; break;
        case 'b':  result += '\b'; break;
        case 'f':  result += '\f'; break;
        case 'n':  result += '\n'; break;
        case 'r':  result += '\r'; break;
        case 't':  result += '\t'; break;
        case 'u':
        {
          unsigned int cu = 0;
          if (!parseHex4(s, i + 1, cu))
          {
            result += "\\u"; // not an escape we understand: keep it as written
            break;
          }
          i += 4;
          if (cu >= 0xD800 && cu <= 0xDBFF) // high surrogate: needs a following \uDC00-\uDFFF
          {
            unsigned int low = 0;
            if (i + 2 < s.size() && s[i + 1] == '\\' && s[i + 2] == 'u' && parseHex4(s, i + 3, low)
                && low >= 0xDC00 && low <= 0xDFFF)
            {
              appendUtf8(result, 0x10000 + ((cu - 0xD800) << 10) + (low - 0xDC00));
              i += 6;
            }
            else
            {
              appendUtf8(result, 0xFFFD);
            }
          }
          else if (cu >= 0xDC00 && cu <= 0xDFFF) // lone low surrogate
          {
            appendUtf8(result, 0xFFFD);
          }
          else
          {
            appendUtf8(result, cu);
          }
          break;
        }
        default: result += c; break; // not valid JSON; keep the character
      }
    }
    return result;
  }

  /// Lenient reader over one JSON document (see the comment above).
  struct JsonReader
  {
    const std::string& s;
    size_t pos = 0;
    size_t skipped = 0;     ///< values of an unexpected type, and entries dropped because of them
    bool malformed = false; ///< the parse stopped early
    size_t stopped_at = 0;  ///< where it stopped (for the warning)

    void skipWhitespace()
    {
      while (pos < s.size() && (s[pos] == ' ' || s[pos] == '\t' || s[pos] == '\n' || s[pos] == '\r')) ++pos;
    }

    bool at(char c)
    {
      skipWhitespace();
      return pos < s.size() && s[pos] == c;
    }

    /// Stop reading: mark the document malformed at the current position.
    void stop()
    {
      markMalformed();
      pos = s.size();
    }

    void markMalformed()
    {
      if (!malformed) stopped_at = pos;
      malformed = true;
    }

    /// The string at pos (which must be at its opening quote), advancing past it.
    std::string string()
    {
      const size_t start = pos;
      ++pos; // opening quote
      std::string raw;
      while (pos < s.size() && s[pos] != '"')
      {
        if (s[pos] == '\\' && pos + 1 < s.size())
        {
          raw += s[pos++];
        }
        raw += s[pos++];
      }
      if (pos >= s.size())
      {
        pos = start;
        stop(); // unterminated string: report where it starts
        return unescapeJsonString(raw);
      }
      ++pos; // closing quote
      return unescapeJsonString(raw);
    }

    /// Skip the value of any kind at pos. Does not advance at ',', ']' or '}'.
    void skipValue()
    {
      skipWhitespace();
      if (pos >= s.size()) return;
      if (s[pos] == '"')
      {
        string();
        return;
      }
      if (s[pos] == '[' || s[pos] == '{')
      {
        const size_t start = pos;
        int depth = 0;
        while (pos < s.size())
        {
          const char c = s[pos];
          if (c == '"') { string(); continue; }
          if (c == '[' || c == '{') ++depth;
          else if (c == ']' || c == '}')
          {
            if (--depth == 0) { ++pos; return; }
          }
          ++pos;
        }
        pos = start;
        stop(); // unbalanced: report where the value starts
        return;
      }
      // number, true, false, null: up to the next separator
      while (pos < s.size() && s[pos] != ',' && s[pos] != '}' && s[pos] != ']') ++pos;
    }

    /// The string value at pos into @p out; for any other kind of value, skip it and return false.
    bool stringOrSkip(std::string& out)
    {
      if (at('"'))
      {
        out = string();
        return true;
      }
      skipValue();
      return false;
    }

    /// A non-negative integer value at pos into @p out; for anything else (negative, fractional,
    /// out of range, not a number), skip it and return false.
    bool unsignedOrSkip(UInt64& out)
    {
      skipWhitespace();
      if (pos < s.size() && s[pos] >= '0' && s[pos] <= '9')
      {
        const char* first = s.data() + pos;
        const char* last = s.data() + s.size();
        const auto [ptr, ec] = std::from_chars(first, last, out);
        pos += static_cast<size_t>(ptr - first);
        const bool integer_ends = pos >= s.size() || (s[pos] != '.' && s[pos] != 'e' && s[pos] != 'E');
        if (ec == std::errc() && integer_ends) return true;
      }
      skipValue();
      return false;
    }

    /// Call @p member(key) for each member of the object at pos (which must be at its '{'). A
    /// value that @p member leaves unread is skipped.
    template <typename F>
    void forEachMember(F member)
    {
      ++pos; // '{'
      while (true)
      {
        skipWhitespace();
        if (pos >= s.size()) { markMalformed(); return; } // unterminated object
        if (s[pos] == '}') { ++pos; return; }
        if (s[pos] == ',') { ++pos; continue; }
        if (s[pos] != '"') { stop(); return; } // key is not a string
        const std::string key = string();
        if (!at(':')) { stop(); return; }
        ++pos;
        skipWhitespace();
        const size_t value_start = pos;
        member(key);
        if (malformed) return;
        if (pos == value_start) skipValue();
        if (pos == value_start) { stop(); return; } // e.g. a stray '}' or ']' as the value
      }
    }

    /// Call @p element() for each element of the array at pos (which must be at its '['). An
    /// element that @p element leaves unread is skipped.
    template <typename F>
    void forEachElement(F element)
    {
      ++pos; // '['
      while (true)
      {
        skipWhitespace();
        if (pos >= s.size()) { markMalformed(); return; } // unterminated array
        if (s[pos] == ']') { ++pos; return; }
        if (s[pos] == ',') { ++pos; continue; }
        const size_t value_start = pos;
        element();
        if (malformed) return;
        if (pos == value_start) skipValue();
        if (pos == value_start) { stop(); return; } // e.g. a stray '}'
      }
    }

    /// Log one warning for this document if anything was skipped or the parse stopped early.
    void report(const char* context, const char* key) const
    {
      if (skipped == 0 && !malformed) return;
      OPENMS_LOG_WARN << context << ": '" << key << "' metadata JSON: ";
      if (skipped > 0)
      {
        OPENMS_LOG_WARN << skipped << " value(s) of an unexpected type (and the entries that needed them) ignored";
      }
      if (malformed)
      {
        OPENMS_LOG_WARN << (skipped > 0 ? "; " : "") << "malformed JSON, reading stopped at character "
                        << stopped_at << ", later entries not read";
      }
      OPENMS_LOG_WARN << std::endl;
    }
  };

  /// Store metavalue @p name = @p value of declared @p type in @p target. A value that does not
  /// parse as its declared type is stored as the raw string.
  void setMetaValueFromStrings(const std::string& name, const std::string& value, const std::string& type,
                               MetaInfoInterface& target)
  {
    const auto unbracketed = [&value]()
    {
      std::string s(value);
      if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
      return s;
    };
    try
    {
      if (type == "int")
      {
        target.setMetaValue(name, DataValue(std::stoll(value)));
      }
      else if (type == "double" || type == "float")
      {
        target.setMetaValue(name, DataValue(StringUtils::toDouble(value)));
      }
      else if (type == "int_list")
      {
        target.setMetaValue(name, DataValue(ListUtils::create<Int>(unbracketed())));
      }
      else if (type == "double_list")
      {
        target.setMetaValue(name, DataValue(ListUtils::create<double>(unbracketed())));
      }
      else if (type == "string_list")
      {
        auto sl = ListUtils::create<std::string>(unbracketed());
        for (auto& e : sl) { StringUtils::trim(e); }
        target.setMetaValue(name, DataValue(sl));
      }
      else
      {
        target.setMetaValue(name, DataValue(value));
      }
    }
    catch (...)
    {
      target.setMetaValue(name, DataValue(value));
    }
  }

  /// Read the metavalue array at pos (as written by serializeMetaValues()) into @p target. An entry
  /// without a string "name" and a string "value" is dropped.
  void readMetaValueArray(JsonReader& r, MetaInfoInterface& target)
  {
    if (!r.at('['))
    {
      r.skipValue();
      ++r.skipped;
      return;
    }
    r.forEachElement([&]()
    {
      if (!r.at('{'))
      {
        r.skipValue();
        ++r.skipped;
        return;
      }
      std::string name, value, type;
      bool name_ok = false, value_ok = false;
      r.forEachMember([&](const std::string& key)
      {
        if (key == "name") name_ok = r.stringOrSkip(name);
        else if (key == "value") value_ok = r.stringOrSkip(value);
        else if (key == "type" && !r.stringOrSkip(type)) ++r.skipped;
      });
      if (!name_ok || name.empty() || !value_ok)
      {
        ++r.skipped;
        return;
      }
      setMetaValueFromStrings(name, value, type, target);
    });
  }
} // namespace

std::vector<DataProcessing> deserializeDataProcessing(const std::string& json, const char* context)
{
  std::vector<DataProcessing> result;
  JsonReader r{json};
  r.skipWhitespace();
  if (r.pos >= json.size()) return result; // empty
  if (!r.at('['))
  {
    r.stop();
    r.report(context, "data_processing");
    return result;
  }

  r.forEachElement([&]()
  {
    if (!r.at('{'))
    {
      r.skipValue();
      ++r.skipped;
      return;
    }
    DataProcessing dp;
    std::string software_name, software_version, completion_time;
    std::vector<std::string> action_names;
    r.forEachMember([&](const std::string& key)
    {
      if (key == "software_name")
      {
        if (!r.stringOrSkip(software_name)) ++r.skipped;
      }
      else if (key == "software_version")
      {
        if (!r.stringOrSkip(software_version)) ++r.skipped;
      }
      else if (key == "completion_time")
      {
        if (!r.stringOrSkip(completion_time)) ++r.skipped;
      }
      else if (key == "actions")
      {
        if (!r.at('['))
        {
          r.skipValue();
          ++r.skipped;
          return;
        }
        r.forEachElement([&]()
        {
          std::string name;
          if (r.stringOrSkip(name)) action_names.push_back(std::move(name));
          else ++r.skipped;
        });
      }
      else if (key == "metavalues")
      {
        readMetaValueArray(r, dp);
      }
    });

    // a record cut short by malformed input keeps what was read before that spot
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
  });

  r.report(context, "data_processing");
  return result;
}

void deserializeMetaValues(const std::string& json, MetaInfoInterface& target, const char* context, const char* key)
{
  JsonReader r{json};
  r.skipWhitespace();
  if (r.pos >= json.size()) return; // empty
  readMetaValueArray(r, target);
  r.report(context, key);
}

ConsensusMap::ColumnHeaders deserializeColumnHeaders(const std::string& json, const char* context)
{
  ConsensusMap::ColumnHeaders result;
  JsonReader r{json};
  r.skipWhitespace();
  if (r.pos >= json.size()) return result; // empty
  if (!r.at('['))
  {
    r.stop();
    r.report(context, "column_headers");
    return result;
  }

  r.forEachElement([&]()
  {
    if (!r.at('{'))
    {
      r.skipValue();
      ++r.skipped;
      return;
    }
    ConsensusMap::ColumnHeader header;
    UInt64 map_index = 0;
    bool has_map_index = false;
    r.forEachMember([&](const std::string& key)
    {
      if (key == "map_index")
      {
        has_map_index = r.unsignedOrSkip(map_index);
        if (!has_map_index) ++r.skipped;
      }
      else if (key == "size")
      {
        UInt64 v = 0;
        if (r.unsignedOrSkip(v)) header.size = v;
        else ++r.skipped;
      }
      else if (key == "unique_id")
      {
        UInt64 v = 0;
        if (r.unsignedOrSkip(v)) header.unique_id = v;
        else ++r.skipped;
      }
      else if (key == "filename")
      {
        if (!r.stringOrSkip(header.filename)) ++r.skipped;
      }
      else if (key == "label")
      {
        if (!r.stringOrSkip(header.label)) ++r.skipped;
      }
      else if (key == "metavalues")
      {
        readMetaValueArray(r, header);
      }
    });
    if (!has_map_index)
    {
      ++r.skipped; // without its index, the header cannot be placed: drop it
      return;
    }
    result[map_index] = std::move(header);
  });

  r.report(context, "column_headers");
  return result;
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
