// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/QPXIdentity.h>
#include <OpenMS/FORMAT/QPXValueValidation.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/MetaInfo.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/MS1LabelState.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include <variant>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <arrow/ipc/writer.h>
#include <arrow/util/base64.h>
#include <arrow/table.h>
#include <arrow/util/thread_pool.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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
    const std::string& filename,
    bool extra_columns,
    int threads)
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
    validator.setExtraColumnsAllowed(extra_columns);
    validator.setThreads(threads);
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
  // The identity declaration belongs to the view, so it is stamped here rather than by each
  // exporter: a producer that forgot it would emit ids no reader could re-derive. qpxc reads
  // identity_composite back to re-derive the ids on conversion (its _source_identity_composite),
  // which is what lets an OpenMS collection survive a round trip with its cross-references intact.
  std::string primary_key;
  std::string identity_composite;
  if      (file_type == "feature_file") { primary_key = QPXFeatureSchema::FEATURE_ID; identity_composite = QPXIdentity::FEATURE_COMPOSITE; }
  else if (file_type == "psm_file")     { primary_key = QPXPSMSchema::PSM_ID;         identity_composite = QPXIdentity::PSM_COMPOSITE; }
  else if (file_type == "pg_file")      { primary_key = QPXPgSchema::PG_ID;           identity_composite = QPXIdentity::PG_COMPOSITE; }
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

  if (!primary_key.empty())
  {
    keys.push_back("primary_key");
    values.push_back(primary_key);
    keys.push_back("identity_composite");
    values.push_back(identity_composite);
  }

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

std::vector<std::pair<std::string, std::string>> qpxCvParams(const MetaInfoInterface& hit)
{
  return qpxCvParams(hit, MS1LabelState::Keys());
}

std::vector<std::pair<std::string, std::string>> qpxCvParams(const MetaInfoInterface& hit,
                                                              const MS1LabelState::Keys& keys)
{
  // The label state written by MS1LabeledWorkflow; see MS1LabelState for the meaning of each key.
  // Looked up by index: this runs once per exported row, and a lookup by name takes the registry lock.
  std::vector<std::pair<std::string, std::string>> params;
  const auto add = [&](const std::string& name, UInt index)
  {
    if (index != static_cast<UInt>(-1) && hit.metaValueExists(index))
    {
      params.emplace_back(name, hit.getMetaValue(index).toString());
    }
  };
  add(MS1LabelState::LABELED_SEQUENCE, keys.labeled_sequence);
  add(MS1LabelState::REMOVED_LABELS, keys.removed_labels);
  add(MS1LabelState::CHANNEL, keys.channel);
  return params;
}

std::string qpxRunFileName(const std::string& ms_run_path)
{
  // File::stemName() already maps "" -> "".
  return File::stemName(ms_run_path);
}

std::vector<Int32> qpxScanComponents(const std::string& spectrum_reference)
{
  if (spectrum_reference.empty()) { return {}; }
  const std::string regex_str = SpectrumNativeIDParser::getRegExFromNativeID(spectrum_reference);
  if (regex_str.empty()) { return {}; }

  // getRegExFromNativeID() maps every native-ID convention onto one of a handful of fixed
  // patterns, but compiling a RegularExpression parses the pattern each time -- and this runs once per
  // exported row. Cached per thread rather than shared: the exporters call it from inside an
  // OpenMP region, where one shared cache would need a lock and reintroduce the cost it saves.
  thread_local std::map<std::string, RegularExpression> compiled;
  auto entry = compiled.find(regex_str);
  if (entry == compiled.end()) { entry = compiled.emplace(regex_str, RegularExpression(regex_str)).first; }

  const Int scan = SpectrumNativeIDParser::extractScanNumber(spectrum_reference, entry->second, true);
  if (scan < 0) { return {}; }
  return {static_cast<Int32>(scan)};
}

QPXRunFileNameKeys::QPXRunFileNameKeys() :
  reference_file_name(MetaInfo::registry().getIndex("reference_file_name")),
  run_file_name(MetaInfo::registry().getIndex("run_file_name"))
{
}

std::string qpxPsmRunFileName(const MetaInfoInterface& hit,
                              const MetaInfoInterface& identification,
                              const std::string& resolved_run_file,
                              const QPXRunFileNameKeys& keys)
{
  // UInt(-1) is guarded explicitly: the index overload of metaValueExists() does not special-case
  // the "not registered" sentinel the way the string overload does.
  const auto has = [](const MetaInfoInterface& meta, UInt index)
  { return index != static_cast<UInt>(-1) && meta.metaValueExists(index); };

  std::string run_file;
  if (has(hit, keys.reference_file_name))
  {
    run_file = hit.getMetaValue(keys.reference_file_name).toString();
  }
  else if (has(hit, keys.run_file_name))
  {
    run_file = hit.getMetaValue(keys.run_file_name).toString();
  }
  else if (has(identification, keys.reference_file_name))
  {
    run_file = identification.getMetaValue(keys.reference_file_name).toString();
  }
  if (run_file.empty()) { run_file = resolved_run_file; }

  // Stem every source of the value, not just the fallback, so the column is a usable join key
  // across the psm, feature and pg tables no matter which branch above supplied it.
  return qpxRunFileName(run_file);
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
  return qpxScanFormat(native_ids, 1);
}

std::string qpxScanFormat(const std::vector<std::string>& native_ids, int threads)
{
  // Per range of IDs (one per thread): the first convention seen and the first one that differs from it. Merged in
  // order below, the result and the warning (the first pair of conventions in ID order) are those of one serial pass.
  struct Range { std::string resolved, conflict; };
  int ranges = 1;
#ifdef _OPENMP
  if (threads > 1 && !omp_in_parallel()) { ranges = static_cast<int>(std::min<size_t>(static_cast<size_t>(threads), std::max<size_t>(1, native_ids.size() / 4096))); }
#endif
  std::vector<Range> parts(ranges);
  std::exception_ptr error; // (exceptions must not leave the parallel region)
  #pragma omp parallel for num_threads(ranges) schedule(static) if(ranges > 1)
  for (int r = 0; r < ranges; ++r)
  {
    try
    {
      Range& part = parts[r];
      const size_t begin = native_ids.size() * r / ranges, end = native_ids.size() * (r + 1) / ranges;
      for (size_t i = begin; i < end; ++i)
      {
        const std::string fmt = qpxScanFormat(native_ids[i]);
        if (fmt.empty()) { continue; } // unrecognized IDs carry no evidence either way
        if (part.resolved.empty()) { part.resolved = fmt; }
        else if (part.resolved != fmt) { part.conflict = fmt; break; }
      }
    }
    catch (...)
    {
      #pragma omp critical (ArrowIOHelpers_qpxScanFormat)
      if (!error) { error = std::current_exception(); }
    }
  }
  if (error) { std::rethrow_exception(error); }
  std::string resolved;
  for (const Range& part : parts)
  {
    std::string conflict;
    if (resolved.empty()) { resolved = part.resolved; conflict = part.conflict; }
    else if (!part.resolved.empty() && part.resolved != resolved) { conflict = part.resolved; }
    else { conflict = part.conflict; }
    if (!conflict.empty())
    {
      // Two conventions in one export: neither token would describe the whole file, so
      // omit scan_format rather than assert a wrong one.
      OPENMS_LOG_WARN << "ArrowIOHelpers: spectrum native IDs mix '" << resolved << "' and '"
                      << conflict << "' conventions; omitting QPX scan_format." << std::endl;
      return "";
    }
  }
  return resolved;
}

namespace
{
  /// The rows of @p c (the size of the vector of its type).
  size_t parquetColumnRows(const ParquetColumn& c)
  {
    switch (c.type)
    {
      case ParquetColumn::Type::INT16:
      case ParquetColumn::Type::INT32:
      case ParquetColumn::Type::BOOL:
      case ParquetColumn::Type::INT64:                return c.int64_values.size();
      case ParquetColumn::Type::DOUBLE:
      case ParquetColumn::Type::FLOAT:                return c.double_values.size();
      case ParquetColumn::Type::INT32_LIST:           return c.int32_lists.size();
      case ParquetColumn::Type::STRUCT_LIST:          return c.struct_lists.size();
      case ParquetColumn::Type::UTF8:                 return c.utf8_values.size();
      case ParquetColumn::Type::FLOAT_LIST:           return c.float_lists.size();
      case ParquetColumn::Type::UTF8_LIST:            return c.utf8_lists.size();
      case ParquetColumn::Type::SCORE_LIST:           return c.score_lists.size();
    }
    return 0;
  }

  /// The Arrow type of a STRUCT_LIST of @p fields (null for a field of another type)
  std::shared_ptr<arrow::DataType> structListType(const std::vector<ParquetColumn::StructField>& fields)
  {
    std::vector<std::shared_ptr<arrow::Field>> out;
    for (const ParquetColumn::StructField& f : fields)
    {
      std::shared_ptr<arrow::DataType> type;
      switch (f.type)
      {
        case ParquetColumn::Type::UTF8:        type = arrow::utf8(); break;
        case ParquetColumn::Type::INT32:       type = arrow::int32(); break;
        case ParquetColumn::Type::BOOL:        type = arrow::boolean(); break;
        case ParquetColumn::Type::DOUBLE:      type = arrow::float64(); break;
        case ParquetColumn::Type::FLOAT:       type = arrow::float32(); break;
        case ParquetColumn::Type::STRUCT_LIST: type = structListType(f.fields); break;
        default: break;
      }
      if (!type) { return nullptr; }
      out.push_back(arrow::field(f.name, type, f.nullable));
    }
    return arrow::list(arrow::struct_(out));
  }

  /// Appends a list of @p entries of the structs of @p fields to @p b (whose type is structListType(@p fields))
  arrow::Status appendStructEntries(arrow::ListBuilder& b, const std::vector<ParquetColumn::StructField>& fields,
                                    const std::vector<ParquetColumn::StructEntry>& entries, const std::string& column)
  {
    ARROW_RETURN_NOT_OK(b.Append());
    auto& entry_builder = static_cast<arrow::StructBuilder&>(*b.value_builder());
    for (const ParquetColumn::StructEntry& entry : entries)
    {
      if (entry.size() != fields.size())
      {
        return arrow::Status::Invalid("column '", column, "': an entry of ", entry.size(), " values for ", fields.size(), " fields");
      }
      ARROW_RETURN_NOT_OK(entry_builder.Append());
      for (size_t k = 0; k < fields.size(); ++k)
      {
        const ParquetColumn::StructField& f = fields[k];
        const auto& v = entry[k].value;
        arrow::ArrayBuilder& fb = *entry_builder.field_builder(static_cast<int>(k));
        if (std::holds_alternative<std::monostate>(v))
        {
          if (!f.nullable) { return arrow::Status::Invalid("column '", column, "': the field '", f.name, "' is non-null"); }
          ARROW_RETURN_NOT_OK(fb.AppendNull());
          continue;
        }
        const auto mismatch = [&]() { return arrow::Status::Invalid("column '", column, "': a value of another type for the field '", f.name, "'"); };
        switch (f.type)
        {
          case ParquetColumn::Type::UTF8:
            if (!std::holds_alternative<std::string>(v)) { return mismatch(); }
            ARROW_RETURN_NOT_OK(static_cast<arrow::StringBuilder&>(fb).Append(std::get<std::string>(v)));
            break;
          case ParquetColumn::Type::INT32:
            if (!std::holds_alternative<Int32>(v)) { return mismatch(); }
            ARROW_RETURN_NOT_OK(static_cast<arrow::Int32Builder&>(fb).Append(std::get<Int32>(v)));
            break;
          case ParquetColumn::Type::BOOL:
            if (!std::holds_alternative<bool>(v)) { return mismatch(); }
            ARROW_RETURN_NOT_OK(static_cast<arrow::BooleanBuilder&>(fb).Append(std::get<bool>(v)));
            break;
          case ParquetColumn::Type::DOUBLE:
            if (!std::holds_alternative<double>(v)) { return mismatch(); }
            ARROW_RETURN_NOT_OK(static_cast<arrow::DoubleBuilder&>(fb).Append(std::get<double>(v)));
            break;
          case ParquetColumn::Type::FLOAT:
            if (!std::holds_alternative<double>(v)) { return mismatch(); }
            ARROW_RETURN_NOT_OK(static_cast<arrow::FloatBuilder&>(fb).Append(static_cast<float>(std::get<double>(v))));
            break;
          case ParquetColumn::Type::STRUCT_LIST:
            if (!std::holds_alternative<std::vector<ParquetColumn::StructEntry>>(v)) { return mismatch(); }
            ARROW_RETURN_NOT_OK(appendStructEntries(static_cast<arrow::ListBuilder&>(fb), f.fields,
                                                    std::get<std::vector<ParquetColumn::StructEntry>>(v), column));
            break;
          default:
            return mismatch();
        }
      }
    }
    return arrow::Status::OK();
  }

  /// The Arrow array and field of @p c (@p rows rows); a builder error is returned as its status.
  /// The number of values of all lists of @p lists (for the exact reserve of a list column's values)
  template <typename List>
  int64_t totalSize(const std::vector<List>& lists)
  {
    int64_t n = 0;
    for (const List& list : lists) { n += static_cast<int64_t>(list.size()); }
    return n;
  }

  arrow::Status buildParquetColumn(const ParquetColumn& c, size_t rows, std::shared_ptr<arrow::Field>& field,
                                   std::shared_ptr<arrow::Array>& array)
  {
    auto isNull = [&c](size_t row) { return !c.null_rows.empty() && c.null_rows[row]; };
    // One row of each type: append the value of @p row (not null) to @p b.
    auto build = [&](arrow::ArrayBuilder& b, const std::shared_ptr<arrow::DataType>& type, auto&& append) -> arrow::Status
    {
      ARROW_RETURN_NOT_OK(b.Reserve(static_cast<int64_t>(rows)));
      for (size_t row = 0; row < rows; ++row)
      {
        if (isNull(row) && !c.nullable)
        {
          return arrow::Status::Invalid("column '", c.name, "' is non-null, row ", row, " is null");
        }
        ARROW_RETURN_NOT_OK(isNull(row) ? b.AppendNull() : append(row));
      }
      ARROW_RETURN_NOT_OK(b.Finish(&array));
      field = arrow::field(c.name, type, c.nullable);
      return arrow::Status::OK();
    };
    auto int32Of = [&c](int64_t v, size_t row, int32_t& out) -> arrow::Status
    {
      if (v < std::numeric_limits<int32_t>::min() || v > std::numeric_limits<int32_t>::max())
      {
        return arrow::Status::Invalid("column '", c.name, "' row ", row, ": ", v, " is outside the range of int32");
      }
      out = static_cast<int32_t>(v);
      return arrow::Status::OK();
    };
    switch (c.type)
    {
      case ParquetColumn::Type::INT16:
      {
        arrow::Int16Builder b;
        return build(b, arrow::int16(), [&](size_t row)
        {
          const int64_t v = c.int64_values[row];
          if (v < std::numeric_limits<int16_t>::min() || v > std::numeric_limits<int16_t>::max())
          {
            return arrow::Status::Invalid("column '", c.name, "' row ", row, ": ", v, " is outside the range of int16");
          }
          return b.Append(static_cast<int16_t>(v));
        });
      }
      case ParquetColumn::Type::INT64:
      {
        arrow::Int64Builder b;
        return build(b, arrow::int64(), [&](size_t row) { return b.Append(c.int64_values[row]); });
      }
      case ParquetColumn::Type::INT32:
      {
        arrow::Int32Builder b;
        return build(b, arrow::int32(), [&](size_t row)
        {
          int32_t v = 0;
          ARROW_RETURN_NOT_OK(int32Of(c.int64_values[row], row, v));
          return b.Append(v);
        });
      }
      case ParquetColumn::Type::BOOL:
      {
        arrow::BooleanBuilder b;
        return build(b, arrow::boolean(), [&](size_t row) { return b.Append(c.int64_values[row] != 0); });
      }
      case ParquetColumn::Type::DOUBLE:
      {
        arrow::DoubleBuilder b;
        return build(b, arrow::float64(), [&](size_t row) { return b.Append(c.double_values[row]); });
      }
      case ParquetColumn::Type::FLOAT:
      {
        arrow::FloatBuilder b;
        return build(b, arrow::float32(), [&](size_t row) { return b.Append(static_cast<float>(c.double_values[row])); });
      }
      case ParquetColumn::Type::INT32_LIST:
      {
        auto values = std::make_shared<arrow::Int32Builder>();
        arrow::ListBuilder b(arrow::default_memory_pool(), values);
        ARROW_RETURN_NOT_OK(values->Reserve(totalSize(c.int32_lists)));
        return build(b, arrow::list(arrow::int32()), [&](size_t row)
        {
          ARROW_RETURN_NOT_OK(b.Append());
          return values->AppendValues(c.int32_lists[row]);
        });
      }
      case ParquetColumn::Type::STRUCT_LIST:
      {
        const auto type = structListType(c.struct_fields);
        if (!type) { return arrow::Status::Invalid("column '", c.name, "': a struct field of an unsupported type"); }
        std::unique_ptr<arrow::ArrayBuilder> made;
        ARROW_RETURN_NOT_OK(arrow::MakeBuilder(arrow::default_memory_pool(), type, &made));
        auto& b = static_cast<arrow::ListBuilder&>(*made);
        return build(b, type, [&](size_t row) { return appendStructEntries(b, c.struct_fields, c.struct_lists[row], c.name); });
      }
      case ParquetColumn::Type::UTF8:
      {
        arrow::StringBuilder b;
        int64_t bytes = 0;
        for (const std::string& v : c.utf8_values) { bytes += static_cast<int64_t>(v.size()); }
        ARROW_RETURN_NOT_OK(b.ReserveData(bytes));
        return build(b, arrow::utf8(), [&](size_t row) { return b.Append(c.utf8_values[row]); });
      }
      case ParquetColumn::Type::FLOAT_LIST:
      {
        auto values = std::make_shared<arrow::FloatBuilder>();
        arrow::ListBuilder b(arrow::default_memory_pool(), values);
        ARROW_RETURN_NOT_OK(values->Reserve(totalSize(c.float_lists)));
        return build(b, arrow::list(arrow::float32()), [&](size_t row)
        {
          ARROW_RETURN_NOT_OK(b.Append());
          return values->AppendValues(c.float_lists[row]);
        });
      }
      case ParquetColumn::Type::UTF8_LIST:
      {
        auto values = std::make_shared<arrow::StringBuilder>();
        arrow::ListBuilder b(arrow::default_memory_pool(), values);
        ARROW_RETURN_NOT_OK(values->Reserve(totalSize(c.utf8_lists)));
        int64_t bytes = 0;
        for (const auto& list : c.utf8_lists) { for (const std::string& v : list) { bytes += static_cast<int64_t>(v.size()); } }
        ARROW_RETURN_NOT_OK(values->ReserveData(bytes));
        return build(b, arrow::list(arrow::utf8()), [&](size_t row)
        {
          ARROW_RETURN_NOT_OK(b.Append());
          return values->AppendValues(c.utf8_lists[row]);
        });
      }
      case ParquetColumn::Type::SCORE_LIST:
      {
        // the type of the QPX additional_scores, so that the column can replace that of a QPX view
        const auto type = QPXPSMSchema::additionalScoresType();
        auto name = std::make_shared<arrow::StringBuilder>();
        auto value = std::make_shared<arrow::DoubleBuilder>();
        auto higher_better = std::make_shared<arrow::BooleanBuilder>();
        auto entry = std::make_shared<arrow::StructBuilder>(std::static_pointer_cast<arrow::ListType>(type)->value_type(),
          arrow::default_memory_pool(), std::vector<std::shared_ptr<arrow::ArrayBuilder>>{name, value, higher_better});
        arrow::ListBuilder b(arrow::default_memory_pool(), entry, type);
        return build(b, type, [&](size_t row)
        {
          ARROW_RETURN_NOT_OK(b.Append());
          for (const ParquetColumn::Score& s : c.score_lists[row])
          {
            ARROW_RETURN_NOT_OK(entry->Append());
            ARROW_RETURN_NOT_OK(name->Append(s.name));
            ARROW_RETURN_NOT_OK(value->Append(s.value));
            ARROW_RETURN_NOT_OK(s.higher_better ? higher_better->Append(*s.higher_better) : higher_better->AppendNull());
          }
          return arrow::Status::OK();
        });
      }
    }
    return arrow::Status::Invalid("unknown column type");
  }
} // namespace

std::shared_ptr<arrow::Table> withColumns(
  const std::shared_ptr<arrow::Table>& table,
  const std::vector<ParquetColumn>& columns)
{
  if (!table)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "withColumns: null table");
  }
  const size_t rows = static_cast<size_t>(table->num_rows());
  for (const ParquetColumn& c : columns)
  {
    if (parquetColumnRows(c) != rows || (!c.null_rows.empty() && c.null_rows.size() != rows))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Column '" + c.name + "' has " + std::to_string(parquetColumnRows(c)) + " rows, the table " + std::to_string(rows));
    }
  }

  std::shared_ptr<arrow::Table> out = table;
  for (const ParquetColumn& c : columns)
  {
    std::shared_ptr<arrow::Field> field;
    std::shared_ptr<arrow::Array> array;
    arrow::Status status = buildParquetColumn(c, rows, field, array);
    if (status.ok())
    {
      const auto chunked = std::make_shared<arrow::ChunkedArray>(array);
      const int existing = out->schema()->GetFieldIndex(c.name);
      auto result = existing >= 0 ? out->SetColumn(existing, field, chunked) : out->AddColumn(out->num_columns(), field, chunked);
      status = result.status();
      if (status.ok()) { out = result.ValueOrDie(); }
    }
    if (!status.ok())
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: adding column '" << c.name << "' failed: " << status.ToString() << std::endl;
      return nullptr;
    }
  }
  return out;
}

std::shared_ptr<arrow::Table> tableOfColumns(
  const std::vector<ParquetColumn>& columns,
  Size rows)
{
  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;
  for (const ParquetColumn& c : columns)
  {
    if (parquetColumnRows(c) != rows || (!c.null_rows.empty() && c.null_rows.size() != rows))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Column '" + c.name + "' has " + std::to_string(parquetColumnRows(c)) + " rows, the table " + std::to_string(rows));
    }
    std::shared_ptr<arrow::Field> field;
    std::shared_ptr<arrow::Array> array;
    const arrow::Status status = buildParquetColumn(c, rows, field, array);
    if (!status.ok())
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: building column '" << c.name << "' failed: " << status.ToString() << std::endl;
      return nullptr;
    }
    fields.push_back(field);
    arrays.push_back(array);
  }
  return arrow::Table::Make(arrow::schema(fields), arrays, static_cast<int64_t>(rows));
}

namespace
{
  /// The writer properties of @p config for the fields of @p schema: compression, page size, statistics, and no
  /// dictionary encoding for the leaves of the plain_columns and for every floating-point leaf (m/z, intensities,
  /// scores: mostly distinct values, whose dictionary page each row group would otherwise try and drop again)
  std::shared_ptr<parquet::WriterProperties> writerProperties(const arrow::Schema& schema, const ParquetWriteConfig& config)
  {
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
    std::function<void(const std::string&, const std::shared_ptr<arrow::DataType>&, bool)> leaves =
      [&](const std::string& path, const std::shared_ptr<arrow::DataType>& type, bool plain)
    {
      if (type->id() == arrow::Type::LIST)
      {
        leaves(path + ".list.element", std::static_pointer_cast<arrow::ListType>(type)->value_type(), plain);
      }
      else if (type->id() == arrow::Type::STRUCT)
      {
        for (const auto& child : type->fields()) { leaves(path + "." + child->name(), child->type(), plain); }
      }
      else if (plain || arrow::is_floating(type->id())) { builder.disable_dictionary(path); }
    };
    for (const auto& field : schema.fields())
    {
      const bool plain = std::find(config.plain_columns.begin(), config.plain_columns.end(), field->name()) != config.plain_columns.end();
      leaves(field->name(), field->type(), plain);
    }
    return builder.build();
  }

  /// The parallel write (writeTableToParquet() with threads > 1, ParquetTableStream): the row groups WriteTable cuts
  /// (row_group_size rows each, at most the writer's max_row_group_length, whatever the chunks of the tables are) as
  /// buffered row groups whose columns are encoded in parallel (use_threads with an executor of its own: Arrow's global
  /// CPU pool has as many threads as the machine has hardware threads, whatever -threads says). Same encodings,
  /// dictionary rules, compression and statistics as the serial write. A group continues across the tables of a
  /// stream. The pool outlives the writer (Close() joins the column tasks; the writer is destroyed before the pool is
  /// shut down).
  class BufferedGroupWriter
  {
  public:
    ~BufferedGroupWriter()
    {
      writer_.reset();
      if (pool_) { (void)pool_->Shutdown(); }
    }

    arrow::Status open(const std::shared_ptr<arrow::Schema>& schema, const std::shared_ptr<arrow::io::OutputStream>& sink,
                       const ParquetWriteConfig& config)
    {
      ARROW_ASSIGN_OR_RAISE(pool_, arrow::internal::ThreadPool::Make(std::max(1, config.threads)));
      const auto writer_properties = writerProperties(*schema, config);
      auto arrow_properties = parquet::ArrowWriterProperties::Builder().store_schema()->set_use_threads(true)
                                ->set_executor(pool_.get())->build();
      ARROW_ASSIGN_OR_RAISE(writer_, parquet::arrow::FileWriter::Open(*schema, arrow::default_memory_pool(), sink,
                                                                      writer_properties, arrow_properties));
      group_rows_ = std::max<int64_t>(1, std::min(config.row_group_size, writer_properties->max_row_group_length()));
      in_group_ = group_rows_; // (full: the first rows open a group)
      return arrow::Status::OK();
    }

    /// The pool of the column tasks has @p threads from now on (at least one)
    arrow::Status setThreads(int threads)
    {
      return pool_->SetCapacity(std::max(1, threads));
    }

    /// Appends the rows of @p table: a batch is a chunk or the rest of a group
    arrow::Status write(const arrow::Table& table)
    {
      arrow::TableBatchReader reader(table);
      reader.set_chunksize(group_rows_);
      std::shared_ptr<arrow::RecordBatch> batch;
      while (true)
      {
        ARROW_RETURN_NOT_OK(reader.ReadNext(&batch));
        if (!batch) { break; }
        for (int64_t offset = 0; offset < batch->num_rows();)
        {
          if (in_group_ == group_rows_)
          {
            ARROW_RETURN_NOT_OK(writer_->NewBufferedRowGroup());
            in_group_ = 0;
          }
          const int64_t rows = std::min(batch->num_rows() - offset, group_rows_ - in_group_);
          ARROW_RETURN_NOT_OK(writer_->WriteRecordBatch(*batch->Slice(offset, rows)));
          in_group_ += rows;
          offset += rows;
        }
      }
      return arrow::Status::OK();
    }

    /// Adds @p keys to the footer (written by close()), after the keys of the schema: in the file's key-value
    /// metadata, and in the Arrow schema stored in it (ARROW:schema, serialized at open() by store_schema; a reader
    /// takes the schema metadata from there, so it is serialized again with the keys), as if the schema had them
    arrow::Status addKeyValueMetadata(const std::shared_ptr<arrow::Schema>& schema, const std::map<std::string, std::string>& keys)
    {
      if (keys.empty()) { return arrow::Status::OK(); }
      auto metadata = schema->metadata() ? schema->metadata()->Copy() : std::make_shared<arrow::KeyValueMetadata>();
      for (const auto& [k, v] : keys) { if (metadata->Contains(k)) { ARROW_RETURN_NOT_OK(metadata->Set(k, v)); } else { metadata->Append(k, v); } }
      ARROW_ASSIGN_OR_RAISE(std::shared_ptr<arrow::Buffer> serialized,
                            arrow::ipc::SerializeSchema(*schema->WithMetadata(metadata), arrow::default_memory_pool()));
      auto added = std::make_shared<arrow::KeyValueMetadata>();
      for (const auto& [k, v] : keys) { added->Append(k, v); }
      added->Append("ARROW:schema", arrow::util::base64_encode(std::string_view(reinterpret_cast<const char*>(serialized->data()), static_cast<size_t>(serialized->size()))));
      return writer_->AddKeyValueMetadata(added);
    }

    arrow::Status close()
    {
      return writer_->Close();
    }

  private:
    std::shared_ptr<arrow::internal::ThreadPool> pool_;
    std::unique_ptr<parquet::arrow::FileWriter> writer_;
    int64_t group_rows_ = 0, in_group_ = 0;
  };

  /// Removes @p filename after a failed write (logged): FileOutputStream::Open created (and truncated) the file, and
  /// a truncated Parquet file has no footer, so a reader reports it as corrupt rather than as the smaller table it
  /// looks like. Closed before: on Windows an open handle blocks the unlink.
  bool abandonFile(const std::shared_ptr<arrow::io::FileOutputStream>& outfile, const std::string& filename, const std::string& what)
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: " << what << std::endl;
    if (outfile) { (void)outfile->Close(); }
    if (File::exists(filename) && !File::remove(filename))
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to remove incomplete output " << filename << std::endl;
    }
    return false;
  }
} // namespace

bool validateQPXTable(
  const std::shared_ptr<arrow::Table>& table,
  const std::string& filename,
  const ParquetWriteConfig& config,
  const std::shared_ptr<const arrow::KeyValueMetadata>& metadata)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: null table passed to validateQPXTable (" << filename << ")" << std::endl;
    return false;
  }
  return validateQPXValuesForWrite(metadata ? table->ReplaceSchemaMetadata(metadata) : table, filename,
                                   config.qpx_extra_columns, config.threads);
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
  if (!config.validated && !validateQPXValuesForWrite(table, filename, config.qpx_extra_columns, config.threads)) { return false; }

  auto file_result = arrow::io::FileOutputStream::Open(filename);
  if (!file_result.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to open " << filename
                     << " for writing: " << file_result.status().ToString() << std::endl;
    return false;
  }
  const auto& outfile = file_result.ValueOrDie();

  arrow::Status status;
  try
  {
    if (config.threads <= 1 || table->num_rows() == 0 || config.row_group_size <= 0) // (WriteTable rejects the size)
    {
      auto arrow_properties = parquet::ArrowWriterProperties::Builder().store_schema()->build();
      status = parquet::arrow::WriteTable(
        *table,
        arrow::default_memory_pool(),
        outfile,
        config.row_group_size,
        writerProperties(*table->schema(), config),
        arrow_properties);
    }
    else
    {
      BufferedGroupWriter writer;
      status = writer.open(table->schema(), outfile, config);
      if (status.ok()) { status = writer.write(*table); }
      if (status.ok()) { status = writer.close(); }
    }
  }
  catch (...) // (e.g. std::bad_alloc: no partial file is left behind, which a staged write would publish)
  {
    abandonFile(outfile, filename, "Failed to write " + filename + " (exception)");
    throw;
  }

  if (!status.ok())
  {
    return abandonFile(outfile, filename, "Failed to write " + filename + ": " + status.ToString());
  }

  auto close_status = outfile->Close();
  if (!close_status.ok())
  {
    return abandonFile(outfile, filename, "Failed to close " + filename + ": " + close_status.ToString());
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

  std::shared_ptr<arrow::Table> table;
  if (tables.size() == 1 && tables.front()) // (one table: as it is)
  {
    table = tables.front();
  }
  else
  {
    auto concat_result = arrow::ConcatenateTables(tables);
    if (!concat_result.ok())
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to concatenate tables for " << filename
                       << ": " << concat_result.status().ToString() << std::endl;
      return false;
    }
    table = concat_result.ValueOrDie();
  }
  // The concatenated table inherits the first input's schema metadata, which carries that
  // input's uuid/creation_date. Stamp one identity for the merged file instead.
  if (metadata) { table = table->ReplaceSchemaMetadata(metadata); }
  return writeTableToParquet(table, filename, config);
}

// ---------------------------------------------------------------------------
// ParquetTableStream
// ---------------------------------------------------------------------------

struct ParquetTableStream::Impl
{
  std::string filename;
  ParquetWriteConfig config;
  std::shared_ptr<const arrow::KeyValueMetadata> metadata;
  std::shared_ptr<arrow::Schema> schema; // of the file (the first table's, with the metadata), once open
  std::shared_ptr<arrow::io::FileOutputStream> outfile;
  std::unique_ptr<BufferedGroupWriter> writer;
  std::optional<QPXValueValidation> validator; // for a QPX view (the metadata's file_type)
  Size rows = 0;
  bool failed = false, closed = false, opened = false; // (opened: the file exists, created by append())

  /// The validator of the QPX view the metadata names, if any: as validateQPXValuesForWrite() decides it
  bool makeValidator()
  {
    const auto& meta = schema->metadata();
    if (!meta || !meta->Contains("qpx_version")) { return true; }
    const auto file_type = meta->Get("file_type");
    if (!file_type.ok())
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: QPX table for " << filename << " has no file_type metadata." << std::endl;
      return false;
    }
    QPXValueValidation::View view;
    if      (*file_type == "psm_file")     { view = QPXValueValidation::View::PSM; }
    else if (*file_type == "feature_file") { view = QPXValueValidation::View::FEATURE; }
    else if (*file_type == "pg_file")      { view = QPXValueValidation::View::PROTEIN_GROUP; }
    else { return true; }
    validator.emplace(view);
    validator->setExtraColumnsAllowed(config.qpx_extra_columns);
    return true;
  }

  /// The stream has failed (terminal): the writer is destroyed before the file is closed and removed (logged)
  bool fail(const std::string& what)
  {
    failed = true;
    writer.reset();
    if (opened) { abandonFile(outfile, filename, what); }
    else { OPENMS_LOG_ERROR << "ArrowIOHelpers: " << what << std::endl; }
    outfile.reset();
    return false;
  }
};

ParquetTableStream::ParquetTableStream(std::string filename, const ParquetWriteConfig& config,
                                       std::shared_ptr<const arrow::KeyValueMetadata> metadata) :
  impl_(std::make_unique<Impl>())
{
  impl_->filename = std::move(filename);
  impl_->config = config;
  impl_->config.threads = std::max(1, config.threads);
  impl_->metadata = std::move(metadata);
}

ParquetTableStream::~ParquetTableStream()
{
  abandon();
}

bool ParquetTableStream::append(const std::shared_ptr<arrow::Table>& table, int threads)
{
  Impl& s = *impl_;
  if (s.failed || s.closed)
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: no table can be appended to " << s.filename << (s.closed ? " (closed)" : " (failed)") << std::endl;
    return false;
  }
  if (!table) { return s.fail("null table passed to ParquetTableStream (" + s.filename + ")"); }
  if (threads > 0) { s.config.threads = threads; }
  try // (every failure ends the stream: an exception too, e.g. std::bad_alloc, so that no later table opens the file anew)
  {
    if (!s.writer)
    {
      s.schema = s.metadata ? table->schema()->WithMetadata(s.metadata) : table->schema();
      if (!s.makeValidator()) { return s.fail("Failed to write " + s.filename + ": no QPX view"); }
      auto file_result = arrow::io::FileOutputStream::Open(s.filename);
      if (!file_result.ok())
      {
        return s.fail("Failed to open " + s.filename + " for writing: " + file_result.status().ToString());
      }
      s.outfile = file_result.ValueOrDie();
      s.opened = true;
      s.writer = std::make_unique<BufferedGroupWriter>();
      const arrow::Status status = s.writer->open(s.schema, s.outfile, s.config);
      if (!status.ok()) { return s.fail("Failed to write " + s.filename + ": " + status.ToString()); }
    }
    else if (!table->schema()->Equals(*s.schema, /*check_metadata=*/false))
    {
      return s.fail("Failed to write " + s.filename + ": the schema of a table differs from the first table's");
    }
    if (s.validator)
    {
      s.validator->setThreads(s.config.threads);
      const auto validation = s.validator->validate(table);
      if (!validation.valid)
      {
        return s.fail("refusing invalid QPX " + s.schema->metadata()->Get("file_type").ValueOr("") + " table for " + s.filename + ": " + validation.toString());
      }
    }
    arrow::Status status = s.writer->setThreads(s.config.threads);
    if (status.ok()) { status = s.writer->write(*table); }
    if (!status.ok()) { return s.fail("Failed to write " + s.filename + ": " + status.ToString()); }
    s.rows += static_cast<Size>(table->num_rows());
    return true;
  }
  catch (const std::exception& e)
  {
    return s.fail("Failed to write " + s.filename + " (exception): " + e.what());
  }
}

bool ParquetTableStream::close(const std::map<std::string, std::string>& keys)
{
  Impl& s = *impl_;
  if (s.failed || s.closed) { return !s.failed; }
  if (!s.writer) { s.closed = true; return true; } // (nothing appended: no file)
  try
  {
    arrow::Status status = s.writer->addKeyValueMetadata(s.schema, keys);
    if (status.ok()) { status = s.writer->close(); }
    if (!status.ok()) { return s.fail("Failed to write " + s.filename + ": " + status.ToString()); }
    s.writer.reset();
    status = s.outfile->Close();
    if (!status.ok()) { return s.fail("Failed to close " + s.filename + ": " + status.ToString()); }
  }
  catch (const std::exception& e)
  {
    return s.fail("Failed to write " + s.filename + " (exception): " + e.what());
  }
  s.outfile.reset();
  s.closed = true;
  return true;
}

void ParquetTableStream::abandon()
{
  Impl& s = *impl_;
  if (s.closed || s.failed) { return; }
  s.failed = true;
  s.writer.reset();
  if (s.outfile) { (void)s.outfile->Close(); }
  s.outfile.reset();
  if (s.opened && File::exists(s.filename) && !File::remove(s.filename))
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to remove incomplete output " << s.filename << std::endl;
  }
}

Size ParquetTableStream::rows() const { return impl_->rows; }
bool ParquetTableStream::failed() const { return impl_->failed; }
bool ParquetTableStream::closed() const { return impl_->closed; }
const std::string& ParquetTableStream::filename() const { return impl_->filename; }

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

std::string metaValueString(const DataValue& value)
{
  // Shortest decimal that reads back as the same double. NaN as DataValue::toString() writes it, since
  // StringUtils::toDouble() reads "NaN" and "inf"/"-inf" but std::to_chars writes "nan".
  const auto append_double = [](double d, std::string& target)
  {
    if (std::isnan(d)) { target += "NaN"; return; }
    char buf[32];
    const auto result = std::to_chars(buf, buf + sizeof(buf), d);
    target.append(buf, result.ptr);
  };

  std::string s;
  if (value.valueType() == DataValue::DOUBLE_VALUE)
  {
    append_double(static_cast<double>(value), s);
  }
  else if (value.valueType() == DataValue::DOUBLE_LIST)
  {
    // the list format of DataValue::toString(), "[a, b, c]"
    s += '[';
    const DoubleList list = value.toDoubleList();
    for (Size i = 0; i < list.size(); ++i)
    {
      if (i > 0) s += ", ";
      append_double(list[i], s);
    }
    s += ']';
  }
  else
  {
    s = value.toString();
  }
  return s;
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
      try { target.setMetaValue(name, StringUtils::toDouble(value_str)); }
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
