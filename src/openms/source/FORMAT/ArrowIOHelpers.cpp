// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/QPXValueValidation.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <map>
#include <set>
#include <unordered_set>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <arrow/table.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

namespace OpenMS
{
namespace ArrowIOHelpers
{

std::string generateUuidV4()
{
  return UniqueIdGenerator::getUUID();
}

Size tableRowCount(const std::shared_ptr<arrow::Table>& table)
{
  return table ? static_cast<Size>(table->num_rows()) : 0;
}

namespace
{
  arrow::Compression::type toArrowCompression(ParquetWriteConfig::Compression c)
  {
    switch (c)
    {
      case ParquetWriteConfig::Compression::NONE:   return arrow::Compression::UNCOMPRESSED;
      case ParquetWriteConfig::Compression::SNAPPY: return arrow::Compression::SNAPPY;
      case ParquetWriteConfig::Compression::GZIP:   return arrow::Compression::GZIP;
      case ParquetWriteConfig::Compression::LZ4:    return arrow::Compression::LZ4;
      case ParquetWriteConfig::Compression::ZSTD:   return arrow::Compression::ZSTD;
    }
    return arrow::Compression::ZSTD;
  }

  /// QPX `compression_format` token for @p c, or "" if QPX does not define one.
  /// The spec's vocabulary is zstd|snappy|gzip|lzo|none — LZ4 has no token.
  std::string qpxCompressionName(ParquetWriteConfig::Compression c)
  {
    switch (c)
    {
      case ParquetWriteConfig::Compression::NONE:   return "none";
      case ParquetWriteConfig::Compression::SNAPPY: return "snappy";
      case ParquetWriteConfig::Compression::GZIP:   return "gzip";
      case ParquetWriteConfig::Compression::ZSTD:   return "zstd";
      case ParquetWriteConfig::Compression::LZ4:    return "";
    }
    return "";
  }

  /// Validate an OpenMS-supported QPX table before a generic Parquet helper opens its output.
  /// Non-QPX tables, and QPX sidecars not produced by these exporters, pass through unchanged.
  bool validateQPXValuesForWrite(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& filename)
  {
    const auto metadata = table->schema()->metadata();
    if (!metadata || !metadata->Contains("qpx_version")) { return true; }

    const auto file_type_result = metadata->Get("file_type");
    if (!file_type_result.ok())
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: QPX table for " << filename
                       << " has no file_type metadata." << std::endl;
      return false;
    }

    const std::string file_type = file_type_result.ValueOrDie();
    QPXValueValidation::View view;
    if (file_type == "psm_file")
    {
      view = QPXValueValidation::View::PSM;
    }
    else if (file_type == "feature_file")
    {
      view = QPXValueValidation::View::FEATURE;
    }
    else if (file_type == "pg_file")
    {
      view = QPXValueValidation::View::PROTEIN_GROUP;
    }
    else
    {
      return true;
    }

    QPXValueValidation validator(view);
    const auto validation = validator.validate(table);
    if (!validation.valid)
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: refusing invalid QPX " << file_type << " table for "
                       << filename << ": " << validation.toString() << std::endl;
      return false;
    }
    return true;
  }
}

std::shared_ptr<const arrow::KeyValueMetadata> qpxFileMetadata(
  const std::string& file_type,
  const ParquetWriteConfig& config,
  const std::map<std::string, std::string>& extra)
{
  const std::string compression = qpxCompressionName(config.compression);
  if (compression.empty())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: QPX does not define a compression_format token for LZ4; "
                        "use ZSTD (default), SNAPPY, GZIP or NONE for QPX output." << std::endl;
    return nullptr;
  }

  std::vector<std::string> keys{
    "qpx_version", "file_type", "creator", "software_provider", "creation_date",
    "compression_format", "uuid"};
  std::vector<std::string> values{
    // QPX 1.1 (bigbio/qpx#220): the pg view is re-keyed from a scalar run_file_name onto
    // grouped_runs (list<string>). Breaking, but shipped as a minor under the spec's pre-2.0
    // stabilisation rule, so a reader must consult the version rather than assume 1.x is additive.
    // The psm and feature views are unchanged by 1.1 and carry the same version key.
    "1.1",
    file_type,
    "OpenMS",
    "OpenMS " + VersionInfo::getVersion(),
    DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ"),
    compression,
    generateUuidV4()};

  for (const auto& [k, v] : extra)
  {
    keys.push_back(k);
    values.push_back(v);
  }

  return std::make_shared<arrow::KeyValueMetadata>(std::move(keys), std::move(values));
}

std::string qpxScanFormat(const std::string& native_id)
{
  if (StringUtils::hasPrefix(native_id, "index=")) { return "index"; }
  if (SpectrumNativeIDParser::isNativeID(native_id)) { return "scan"; }
  return "";
}

std::string qpxIntensityLabel(const std::string& column_label, const std::string& channel_name)
{
  // IsobaricChannelExtractor writes both fields, but older/synthetic ConsensusMaps may carry the
  // complete "<method>_<channel>" identity only in ColumnHeader::label. Recover the channel from
  // that label after validating the method instead of rejecting otherwise unambiguous input.
  const auto sep = column_label.rfind('_');
  const std::string method_name = (sep == std::string::npos) ? "" : column_label.substr(0, sep);
  const auto method = IsobaricQuantitationMethod::methodTypeFromName(method_name);
  std::string resolved_channel_name = channel_name;
  if (resolved_channel_name.empty()
      && method != IsobaricQuantitationMethod::MethodType::UNKNOWN
      && sep + 1 < column_label.size())
  {
    resolved_channel_name = column_label.substr(sep + 1);
  }

  if (!resolved_channel_name.empty()
      && method != IsobaricQuantitationMethod::MethodType::UNKNOWN)
  {
    if (StringUtils::hasPrefix(method_name, "tmt"))
    {
      return "TMT" + resolved_channel_name;
    }
    if (StringUtils::hasPrefix(method_name, "itraq"))
    {
      return "ITRAQ" + resolved_channel_name;
    }
  }

  // No channel => not an isobaric map.
  if (channel_name.empty())
  {
    // ProteomicsLFQ stamps "label-free" on every header; an unset label means the same.
    if (column_label.empty() || column_label == "label-free") { return "LFQ"; }

    std::string lower = column_label;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == "lfq") { return "LFQ"; }
    if (lower == "no_label" || lower == "silac light" || lower == "light")
    {
      return "SILAC light";
    }
    if (lower == "silac medium" || lower == "medium"
        || ((lower.find("arg6") != std::string::npos
             || lower.find("lys4") != std::string::npos
             || lower.find("lys6") != std::string::npos)
            && lower.find("arg10") == std::string::npos
            && lower.find("lys8") == std::string::npos))
    {
      return "SILAC medium";
    }
    if (lower == "silac heavy" || lower == "heavy"
        || lower.find("arg10") != std::string::npos
        || lower.find("lys8") != std::string::npos)
    {
      return "SILAC heavy";
    }

    // FeatureFinderMultiplex uses title case while the shared SDRF/QPX vocabulary is uppercase.
    if (StringUtils::hasPrefix(lower, "dimethyl"))
    {
      const std::string suffix = column_label.substr(std::string("Dimethyl").size());
      const std::string canonical = "DIMETHYL" + suffix;
      if (qpxIsCanonicalIntensityLabel(canonical)) { return canonical; }
    }

    // Already-canonical values can occur after loading a QPX-enriched intermediate.
    if (qpxIsCanonicalIntensityLabel(column_label)) { return column_label; }

    OPENMS_LOG_ERROR << "ArrowIOHelpers: column header label '" << column_label
                     << "' is not in the canonical SDRF/QPX intensity-label vocabulary."
                     << std::endl;
    return "";
  }

  OPENMS_LOG_ERROR << "ArrowIOHelpers: cannot derive a QPX intensity label for channel '"
                   << channel_name << "' — column header label '" << column_label
                   << "' does not name a known isobaric quantitation method." << std::endl;
  return "";
}

namespace
{
  enum class SILACRole
  {
    UNKNOWN,
    LIGHT,
    MEDIUM,
    HEAVY
  };

  SILACRole silacRole(const std::string& label)
  {
    std::string lower = label;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower.empty() || lower == "no_label" || lower == "silac light" || lower == "light")
    {
      return SILACRole::LIGHT;
    }
    const bool medium = lower == "silac medium" || lower == "medium"
                        || lower.find("arg6") != std::string::npos
                        || lower.find("lys4") != std::string::npos
                        || lower.find("lys6") != std::string::npos;
    const bool heavy = lower == "silac heavy" || lower == "heavy"
                       || lower.find("arg10") != std::string::npos
                       || lower.find("lys8") != std::string::npos;
    if (medium == heavy) { return SILACRole::UNKNOWN; }
    return medium ? SILACRole::MEDIUM : SILACRole::HEAVY;
  }

  bool hasSILACMarker(const std::string& label)
  {
    const SILACRole role = silacRole(label);
    if (role == SILACRole::MEDIUM || role == SILACRole::HEAVY) { return true; }

    std::string lower = label;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    // An explicit SILAC label with an unknown role must still enter the map-level validation
    // branch and be rejected. Ordinary labels containing "arg" or "lys" are not SILAC markers.
    return lower.find("silac") != std::string::npos;
  }

  std::string canonicalSILACLabel(SILACRole role)
  {
    switch (role)
    {
      case SILACRole::LIGHT:  return "SILAC light";
      case SILACRole::MEDIUM: return "SILAC medium";
      case SILACRole::HEAVY:  return "SILAC heavy";
      case SILACRole::UNKNOWN: break;
    }
    return "";
  }
}

std::map<std::uint64_t, std::string> qpxIntensityLabels(const ConsensusMap& cmap)
{
  using HeaderRef = std::pair<std::uint64_t, const ConsensusMap::ColumnHeader*>;
  std::map<std::string, std::vector<HeaderRef>> by_source;
  for (const auto& [map_index, header] : cmap.getColumnHeaders())
  {
    by_source[header.filename].emplace_back(map_index, &header);
  }

  std::map<std::uint64_t, std::string> result;
  for (const auto& [source, headers] : by_source)
  {
    bool has_isobaric_channel = false;
    bool has_silac_marker = false;
    const bool supported_silac_plex = headers.size() == 2 || headers.size() == 3;
    bool all_silac_roles_known = true;
    std::vector<SILACRole> roles;
    roles.reserve(headers.size());
    for (const auto& [map_index, header] : headers)
    {
      (void)map_index;
      has_isobaric_channel = has_isobaric_channel
                             || (header->metaValueExists("channel_name")
                                 && !header->getMetaValue("channel_name").toString().empty());
      has_silac_marker = has_silac_marker || hasSILACMarker(header->label);
      const auto role = silacRole(header->label);
      roles.push_back(role);
      all_silac_roles_known = all_silac_roles_known && role != SILACRole::UNKNOWN;
    }

    const bool silac_shaped = !has_isobaric_channel && has_silac_marker;
    if (silac_shaped)
    {
      if (!supported_silac_plex || !all_silac_roles_known)
      {
        OPENMS_LOG_ERROR << "ArrowIOHelpers: source run '" << source
                         << "' has SILAC labels but is not a supported, unambiguous two- or "
                            "three-plex."
                         << std::endl;
        for (const auto& [map_index, header] : headers)
        {
          (void)header;
          result[map_index] = "";
        }
        continue;
      }

      std::vector<SILACRole> resolved_roles = roles;
      bool valid_plex = true;
      if (headers.size() == 2)
      {
        // QPX names the non-light channel of a two-plex "heavy", including Arg6/Lys4
        // experiments where the same mass class would be "medium" in a three-plex.
        const Size light_count = static_cast<Size>(std::count(
          roles.begin(), roles.end(), SILACRole::LIGHT));
        valid_plex = light_count == 1;
        if (valid_plex)
        {
          for (auto& role : resolved_roles)
          {
            if (role != SILACRole::LIGHT) { role = SILACRole::HEAVY; }
          }
        }
      }
      else
      {
        valid_plex = std::count(roles.begin(), roles.end(), SILACRole::LIGHT) == 1
                     && std::count(roles.begin(), roles.end(), SILACRole::MEDIUM) == 1
                     && std::count(roles.begin(), roles.end(), SILACRole::HEAVY) == 1;
      }

      if (!valid_plex)
      {
        OPENMS_LOG_ERROR << "ArrowIOHelpers: source run '" << source
                         << "' has a SILAC-shaped column set that cannot be represented as a "
                         << headers.size() << "-plex with unique channel roles." << std::endl;
        for (const auto& [map_index, header] : headers)
        {
          (void)header;
          result[map_index] = "";
        }
        continue;
      }

      for (Size i = 0; i < headers.size(); ++i)
      {
        result[headers[i].first] = canonicalSILACLabel(resolved_roles[i]);
      }
      continue;
    }

    for (const auto& [map_index, header] : headers)
    {
      const std::string channel_name = header->metaValueExists("channel_name")
                                     ? header->getMetaValue("channel_name").toString() : "";
      result[map_index] = qpxIntensityLabel(header->label, channel_name);
    }
  }
  return result;
}

bool qpxIsCanonicalIntensityLabel(const std::string& label)
{
  static const std::unordered_set<std::string> canonical = []
  {
    std::unordered_set<std::string> labels{
      "LFQ", "SILAC light", "SILAC medium", "SILAC heavy",
      "MTRAQ0", "MTRAQ4", "MTRAQ8",
      "DIMETHYL0", "DIMETHYL2", "DIMETHYL4", "DIMETHYL6", "DIMETHYL8"};

    using MethodType = IsobaricQuantitationMethod::MethodType;
    for (int value = static_cast<int>(MethodType::TMT_6PLEX);
         value < static_cast<int>(MethodType::SIZE_OF_METHODTYPE); ++value)
    {
      const auto type = static_cast<MethodType>(value);
      auto method = IsobaricQuantitationMethod::create(type);
      if (!method) { continue; }
      const std::string method_name(method->getMethodName());
      std::string prefix;
      if (StringUtils::hasPrefix(method_name, "tmt")) { prefix = "TMT"; }
      else if (StringUtils::hasPrefix(method_name, "itraq")) { prefix = "ITRAQ"; }
      else { continue; }
      for (const auto& channel : method->getChannelInformation())
      {
        labels.insert(prefix + channel.name);
      }
    }
    return labels;
  }();
  return canonical.contains(label);
}

std::string qpxRunFileName(const std::string& ms_run_path)
{
  // File::stemName() already maps "" -> "".
  return File::stemName(ms_run_path);
}

bool qpxWarnOnRunNameCollisions(const std::string& context,
                                const std::vector<std::string>& ms_run_paths)
{
  std::map<std::string, std::set<std::string>> stem_to_paths;
  for (const auto& path : ms_run_paths)
  {
    if (path.empty()) { continue; }
    stem_to_paths[qpxRunFileName(path)].insert(path);
  }

  bool unique = true;
  for (const auto& [stem, paths] : stem_to_paths)
  {
    if (paths.size() < 2) { continue; }
    unique = false;
    OPENMS_LOG_WARN << context << ": several MS runs share the QPX run name '" << stem
                    << "' after stripping path and extension: "
                    << ListUtils::concatenate(StringList(paths.begin(), paths.end()), ", ")
                    << ". They cannot be told apart in the exported QPX tables." << std::endl;
  }
  return unique;
}

std::string qpxScanFormat(const std::vector<std::string>& native_ids)
{
  std::string resolved;
  for (const auto& id : native_ids)
  {
    const std::string fmt = qpxScanFormat(id);
    if (fmt.empty()) { continue; } // unrecognized IDs carry no evidence either way
    if (resolved.empty()) { resolved = fmt; }
    else if (resolved != fmt)
    {
      // Two conventions in one export: neither token would describe the whole file, so
      // omit scan_format rather than assert a wrong one.
      OPENMS_LOG_WARN << "ArrowIOHelpers: spectrum native IDs mix '" << resolved << "' and '"
                      << fmt << "' conventions; omitting QPX scan_format." << std::endl;
      return "";
    }
  }
  return resolved;
}

bool writeTableToParquet(
  const std::shared_ptr<arrow::Table>& table,
  const std::string& filename,
  const ParquetWriteConfig& config)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: null table passed to writeTableToParquet (" << filename << ")" << std::endl;
    return false;
  }
  if (!validateQPXValuesForWrite(table, filename)) { return false; }

  auto file_result = arrow::io::FileOutputStream::Open(filename);
  if (!file_result.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to open " << filename
                     << " for writing: " << file_result.status().ToString() << std::endl;
    return false;
  }
  const auto& outfile = file_result.ValueOrDie();

  auto builder = parquet::WriterProperties::Builder();
  builder.compression(toArrowCompression(config.compression));
  if (config.compression == ParquetWriteConfig::Compression::ZSTD ||
      config.compression == ParquetWriteConfig::Compression::GZIP)
  {
    builder.compression_level(config.compression_level);
  }
  builder.data_pagesize(config.data_page_size);
  if (config.write_statistics) { builder.enable_statistics(); }
  else                         { builder.disable_statistics(); }

  auto writer_properties = builder.build();
  auto arrow_properties = parquet::ArrowWriterProperties::Builder().store_schema()->build();

  auto status = parquet::arrow::WriteTable(
    *table,
    arrow::default_memory_pool(),
    outfile,
    config.row_group_size,
    writer_properties,
    arrow_properties);

  // FileOutputStream::Open above already created (and truncated) the file, so any failure from
  // here on leaves a partial .parquet behind -- and a truncated Parquet file has no footer, so a
  // reader reports it as corrupt rather than as the smaller table it looks like. Close before
  // removing: on Windows an open handle blocks the unlink.
  const auto abandon = [&](const std::string& what)
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: " << what << std::endl;
    (void)outfile->Close();
    if (!File::remove(filename))
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to remove incomplete output " << filename
                       << std::endl;
    }
    return false;
  };

  if (!status.ok())
  {
    return abandon("Failed to write " + filename + ": " + status.ToString());
  }

  auto close_status = outfile->Close();
  if (!close_status.ok())
  {
    return abandon("Failed to close " + filename + ": " + close_status.ToString());
  }

  return true;
}

bool concatenateAndWriteToParquet(
  const std::vector<std::shared_ptr<arrow::Table>>& tables,
  const std::string& filename,
  const ParquetWriteConfig& config,
  const std::shared_ptr<const arrow::KeyValueMetadata>& metadata)
{
  if (tables.empty()) { return true; }

  auto concat_result = arrow::ConcatenateTables(tables);
  if (!concat_result.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to concatenate tables for " << filename
                     << ": " << concat_result.status().ToString() << std::endl;
    return false;
  }

  auto table = concat_result.ValueOrDie();
  // The concatenated table inherits the first input's schema metadata, which carries that
  // input's uuid/creation_date. Stamp one identity for the merged file instead.
  if (metadata) { table = table->ReplaceSchemaMetadata(metadata); }
  return writeTableToParquet(table, filename, config);
}

// ---------------------------------------------------------------------------
// Read helpers
// ---------------------------------------------------------------------------

std::shared_ptr<arrow::Array> getColumn(
  const std::shared_ptr<arrow::Table>& table,
  const std::string& name,
  bool required)
{
  auto column = table->GetColumnByName(name);
  if (!column)
  {
    if (required)
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: Missing required column '" << name << "'" << std::endl;
    }
    return nullptr;
  }
  if (column->num_chunks() == 0)
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Column '" << name << "' has no chunks" << std::endl;
    return nullptr;
  }
  if (column->num_chunks() == 1)
  {
    return column->chunk(0);
  }
  auto combined = arrow::Concatenate(column->chunks(), arrow::default_memory_pool());
  if (!combined.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to combine chunks for column '" << name << "'" << std::endl;
    return nullptr;
  }
  return *combined;
}

std::string getStringValue(const std::shared_ptr<arrow::Array>& array, int64_t row)
{
  if (!array || array->IsNull(row)) return "";
  return std::static_pointer_cast<arrow::StringArray>(array)->GetString(row);
}

double getDoubleValue(const std::shared_ptr<arrow::Array>& array, int64_t row, double default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::DoubleArray>(array)->Value(row);
}

float getFloatValue(const std::shared_ptr<arrow::Array>& array, int64_t row, float default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::FloatArray>(array)->Value(row);
}

int32_t getInt32Value(const std::shared_ptr<arrow::Array>& array, int64_t row, int32_t default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::Int32Array>(array)->Value(row);
}

int64_t getInt64Value(const std::shared_ptr<arrow::Array>& array, int64_t row, int64_t default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::Int64Array>(array)->Value(row);
}

bool getBoolValue(const std::shared_ptr<arrow::Array>& array, int64_t row, bool default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::BooleanArray>(array)->Value(row);
}

bool isNull(const std::shared_ptr<arrow::Array>& array, int64_t row)
{
  return !array || array->IsNull(row);
}

void readMetaValues(
  const std::shared_ptr<arrow::Array>& array,
  int64_t row,
  MetaInfoInterface& target,
  const std::unordered_set<std::string>& excluded_keys)
{
  if (!array || array->IsNull(row)) return;
  auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
  const int64_t off = list_arr->value_offset(row);
  const int64_t len = list_arr->value_length(row);
  if (len == 0) return;

  // Index the shared child struct via value_offset/length instead of value_slice(row),
  // which allocates a sliced array on every call (once per metavalue list per row).
  auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->values());
  auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
  auto value_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(1));
  auto type_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(2));

  for (int64_t i = off; i < off + len; ++i)
  {
    std::string name = name_arr->GetString(i);
    if (excluded_keys.contains(name)) continue;

    std::string value_str = value_arr->GetString(i);
    std::string type_str = type_arr->GetString(i);

    if (type_str == "int")
    {
      try { target.setMetaValue(name, static_cast<SignedSize>(std::stoll(value_str))); }
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
        std::string s(value_str);
        if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
        target.setMetaValue(name, DataValue(ListUtils::create<Int>(s)));
      }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else if (type_str == "double_list")
    {
      try
      {
        std::string s(value_str);
        if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
        target.setMetaValue(name, DataValue(ListUtils::create<double>(s)));
      }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else if (type_str == "string_list")
    {
      try
      {
        std::string s(value_str);
        if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
        auto sl = ListUtils::create<std::string>(s);
        for (auto& e : sl) { StringUtils::trim(e); }
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

} // namespace ArrowIOHelpers
} // namespace OpenMS
