// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include "IdentificationDataFileSupport.h"

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <set>
#include <system_error>
#include <thread>

namespace OpenMS
{
namespace
{
  namespace IO = Internal::IdentificationDataIO;
  namespace fs = std::filesystem;
  using ID = IdentificationData;
  using File = IdentificationDataFile;
  using IO::append;
  using IO::appendText;
  using IO::check;
  using IO::invalid;
  using IO::Json;
  using IO::number;
  using IO::optionalNumber;
  using IO::text;
  const std::string FORMAT = "OpenMS.IdentificationData";

  [[noreturn]] void fileError(const fs::path& path, const std::string& message)
  { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path.string(), message); }

  std::shared_ptr<arrow::Schema> matchSchema(const std::vector<std::string>& score_columns);

  bool replaceableBundle(const fs::path& path)
  {
    try
    {
      return File::isNativeFile(path.string());
    }
    catch (const Exception::BaseException&)
    {
      return false; // a malformed manifest is not ours to delete
    }
  }

  /// Renames, retrying briefly: on Windows a scanner or indexer may still hold a handle.
  std::error_code renameWithRetry(const fs::path& from, const fs::path& to)
  {
    std::error_code ec;
    for (int attempt = 0; attempt < 5; ++attempt)
    {
      fs::rename(from, to, ec);
      if (! ec) return ec;
      std::this_thread::sleep_for(std::chrono::milliseconds(50 * (attempt + 1)));
    }
    return ec;
  }

  // All outputs remain private until every table is closed and the manifest is written.
  class StagedDirectory
  {
  public:
    StagedDirectory(const std::string& target, bool replace_existing): target_(target), replace_(replace_existing)
    {
      if (! target_.empty() && ! target_.has_filename()) target_ = target_.parent_path(); // "out.idparquet/"
      if (target_.empty()) invalid("Output destination is empty");
      std::error_code ec;
      const bool exists = fs::exists(target_, ec);
      if (ec) fileError(target_, ec.message());
      if (exists && ! (replace_ && replaceableBundle(target_)))
        invalid(replace_ ? "Output destination exists and is not a native identification bundle" : "Output destination already exists");
      fs::path parent = target_.parent_path();
      if (parent.empty()) parent = ".";
      if (! fs::is_directory(parent, ec)) fileError(parent, "Output parent directory does not exist");
      // A previous interrupted process may leave a private staging directory.
      // Reserve a fresh name atomically without modifying another writer's data.
      for (Size attempt = 0; attempt < 100; ++attempt)
      {
        path = sibling(".tmp-");
        if (fs::create_directory(path, ec)) return;
        if (ec) fileError(path, ec.message());
      }
      invalid("Cannot create staging directory");
    }
    ~StagedDirectory()
    {
      if (! published_)
      {
        std::error_code ec;
        fs::remove_all(path, ec);
      }
    }
    void publish()
    {
      std::error_code ec;
      const bool exists = fs::exists(target_, ec);
      if (ec) fileError(target_, ec.message());
      if (! exists)
      {
        ec = renameWithRetry(path, target_);
        if (ec) fileError(target_, ec.message());
        published_ = true;
        return;
      }
      if (! replace_ || ! replaceableBundle(target_)) invalid("Output destination appeared while writing");
      // Move the previous bundle aside, publish, then delete it. If publication fails, restore it.
      const fs::path previous = sibling(".old-");
      ec = renameWithRetry(target_, previous);
      if (ec) fileError(target_, "Cannot replace the existing bundle: " + ec.message());
      ec = renameWithRetry(path, target_);
      if (ec)
      {
        std::error_code restore;
        fs::rename(previous, target_, restore);
        fileError(target_, ec.message());
      }
      published_ = true;
      fs::remove_all(previous, ec); // best effort; the new bundle is already published
    }
    fs::path path;

  private:
    fs::path sibling(const std::string& infix) const
    {
      fs::path parent = target_.parent_path();
      if (parent.empty()) parent = ".";
      return parent / (target_.filename().string() + infix + std::to_string(UniqueIdGenerator::getUniqueId()));
    }
    fs::path target_;
    bool replace_ = false;
    bool published_ = false;
  };
  using IO::validateJsonText;
  Json readManifest(const std::string& path)
  try
  {
    std::ifstream in(fs::path(path) / "manifest.json", std::ios::binary);
    if (! in) invalid("Cannot open manifest.json");
    Json manifest;
    try
    {
      in >> manifest;
    }
    catch (const Json::exception& e)
    {
      invalid(e.what());
    }
    validateJsonText(manifest);
    if (manifest.at("format") != FORMAT || IO::integer<unsigned>(manifest.at("schema_version")) != 1)
      invalid("Unsupported native identification format");
    if (! manifest.at("runs").is_array() || ! manifest.at("inference").is_array()) invalid("Malformed run/result descriptors");
    std::map<std::string, std::vector<std::pair<UInt64, UInt64>>> ranges;
    std::set<std::pair<std::string, std::string>> partitions;
    // A slice belongs to the run (UUID) or inference result (identifier) that declares it.
    const auto claim = [&](const Json& entry, const std::string& partition, const std::string& table) {
      const auto name = IO::tablePath({}, entry.at("path").get<std::string>()).generic_string();
      if (name != table + ".parquet") invalid("Expected a single named Parquet file per table");
      const auto start = IO::integer<UInt64>(entry.at("start"));
      const auto count = IO::integer<UInt64>(entry.at("count"));
      if (start > std::numeric_limits<UInt64>::max() - count) invalid("Overflowing table range");
      if (entry.at("partition").get<std::string>() != partition || !partitions.emplace(name, partition).second)
        invalid("Invalid or duplicate table partition");
      if (count) ranges[name].emplace_back(start, start + count);
    };
    for (const auto& run : manifest.at("runs"))
    {
      const auto& tables = run.at("tables");
      if (!tables.is_object() || !tables.contains("queries") || !tables.contains("matches")
          || tables.size() != (tables.contains("database_sequences") ? 3u : 2u)) invalid("Invalid run table declarations");
      for (auto table = tables.begin(); table != tables.end(); ++table)
        claim(table.value(), run.at("uuid").get<std::string>(), table.key());
    }
    for (const auto& result : manifest.at("inference"))
    {
      const auto& tables = result.at("tables");
      if (! tables.is_object() || tables.size() != 3) invalid("Inference requires three typed tables");
      for (const auto* name : {"inputs", "proteins", "groups"})
        claim(tables.at(name), result.at("identifier").get<std::string>(), name);
    }
    for (auto& [name, slices] : ranges)
    {
      std::sort(slices.begin(), slices.end());
      for (Size i = 1; i < slices.size(); ++i)
        if (slices[i].first < slices[i - 1].second) invalid("Overlapping table slices");
    }
    std::optional<std::vector<ID::ScoreDefinition>> expected_scores;
    std::optional<UInt32> expected_primary;
    for (const auto& run : manifest.at("runs"))
    {
      if (! run.at("scores").is_array()) invalid("Run score definitions must be an array");
      const auto& primary = run.at("primary_score");
      const auto settings = IO::readSettingsJson(run.at("settings"));
      const bool catalog
        = settings.metaValueExists("identification:catalog") && settings.getMetaValue("identification:catalog").toString() == "true";
      if (catalog && (! run.at("scores").empty() || ! primary.is_null())) invalid("A sequence catalog cannot declare PSM scores");
      if (run.at("scores").empty() && primary.is_null() && (IO::integer<UInt64>(run.at("match_count")) == 0 || catalog)) continue;
      if (primary.is_null()) invalid("Configured run must select a primary PSM score");
      const auto column = IO::integer<UInt32>(primary);
      if (column >= run.at("scores").size()) invalid("Primary PSM score is outside the declared schema");
      std::vector<ID::ScoreDefinition> scores;
      for (const auto& value : run.at("scores"))
        scores.push_back(IO::readScoreJson(value));
      if (expected_scores && (scores != *expected_scores || column != *expected_primary))
        invalid("Ordered PSM score schema or primary selection differs between run descriptors");
      expected_scores = std::move(scores);
      expected_primary = column;
    }
    // Score columns are named in the manifest, in score schema order; readers find them by name.
    const auto& columns = manifest.at("score_columns");
    if (! columns.is_array() || columns.size() != (expected_scores ? expected_scores->size() : 0))
      invalid("Score columns do not match the score schema");
    std::set<std::string> names {IO::partitionColumn("matches")};
    const auto schema = matchSchema(columns.get<std::vector<std::string>>());
    for (const auto& field : schema->fields())
      if (field->name().empty() || ! names.insert(field->name()).second) invalid("Empty or duplicate match column name: " + field->name());
    std::set<std::string> uuids, identifiers;
    for (const auto& run : manifest.at("runs"))
    {
      if (! uuids.insert(run.at("uuid").get<std::string>()).second || ! identifiers.insert(run.at("identifier").get<std::string>()).second)
        invalid("Duplicate run UUID or identifier");
    }
    return manifest;
  }
  catch (const Json::exception& error)
  {
    invalid(error.what());
  }
  std::vector<std::string> scoreColumns(const Json& manifest)
  { return manifest.at("score_columns").get<std::vector<std::string>>(); }
  void writeManifest(const fs::path& directory, const Json& manifest)
  {
    validateJsonText(manifest);
    std::ofstream out(directory / "manifest.json", std::ios::binary);
    if (! out) invalid("Cannot write manifest");
    out << manifest.dump(2) << '\n';
    out.close();
    if (! out) invalid("Failed to close manifest");
  }
  template<class Builder, class T>
  void appendOptional(arrow::ArrayBuilder& builder, const std::optional<T>& v)
  {
    if (v) append<Builder>(builder, *v);
    else
      check(builder.AppendNull());
  }
  std::shared_ptr<arrow::DataType> identityType()
  { return arrow::struct_({arrow::field("database", arrow::utf8(), false), arrow::field("accession", arrow::utf8(), false)}); }
  std::shared_ptr<arrow::Schema> querySchema()
  {
    return arrow::schema({arrow::field("query_id", arrow::uint64(), false), arrow::field("source_id", arrow::uint32(), false),
                          arrow::field("data_id", arrow::utf8(), false), arrow::field("rt", arrow::float64()), arrow::field("mz", arrow::float64()),
                          arrow::field("selected_match_id", arrow::uint64()), arrow::field("metadata", IO::metadataType(), false)});
  }
  std::shared_ptr<arrow::Schema> matchSchema(const std::vector<std::string>& score_columns)
  {
    std::vector<std::shared_ptr<arrow::Field>> fields {
      arrow::field("match_id", arrow::uint64(), false),
      arrow::field("query_id", arrow::uint64(), false),
      arrow::field("representation", arrow::utf8(), false),
      arrow::field("encoding", arrow::uint8(), false),
      arrow::field("charge", arrow::int32(), false),
      arrow::field("calculated_mz", arrow::float64()),
      arrow::field("target_decoy", arrow::uint8(), false),
      arrow::field("name", arrow::utf8(), false),
      arrow::field("formula", arrow::utf8()),
      arrow::field("identifiers", arrow::list(arrow::field("item", identityType(), false)), false),
      arrow::field("adduct", arrow::struct_({arrow::field("name", arrow::utf8(), false), arrow::field("formula", arrow::utf8(), false),
                                             arrow::field("charge", arrow::int32(), false), arrow::field("multiplier", arrow::uint32(), false)})),
      arrow::field(
        "sequence_evidence",
        arrow::list(arrow::field("item",
                                 arrow::struct_({arrow::field("database", arrow::uint32(), false), arrow::field("accession", arrow::utf8(), false),
                                                 arrow::field("start", arrow::uint64()), arrow::field("end", arrow::uint64()),
                                                 arrow::field("before", arrow::utf8(), false), arrow::field("after", arrow::utf8(), false)}),
                                 false)),
        false),
      arrow::field(
        "peak_annotations",
        arrow::list(arrow::field("item",
                                 arrow::struct_({arrow::field("annotation", arrow::utf8(), false), arrow::field("charge", arrow::int32(), false),
                                                 arrow::field("mz", arrow::float64(), false), arrow::field("intensity", arrow::float64(), false)}),
                                 false)),
        false),
      arrow::field("metadata", IO::metadataType(), false)};
    for (const auto& name : score_columns)
      fields.push_back(arrow::field(name, arrow::float64()));
    return arrow::schema(fields);
  }
  struct QueryView
  {
    UInt64 query_id;
    UInt32 source_id;
    const ID::Observation& data;
    std::optional<UInt64> selected_match_id;
  };
  struct MatchView
  {
    UInt64 match_id, query_id;
    const ID::MatchData& data;
    const std::vector<double>& scores;
  };
  std::optional<double> storedScore(double score)
  { return std::isnan(score) ? std::nullopt : std::optional<double>(score); }
  const std::optional<double>& storedScore(const std::optional<double>& score) { return score; }
  template<class Query>
  Size queryBytes(const Query& q)
  { return 64 + q.data.data_id.size() + IO::metadataBytes(q.data); }
  template<class Match>
  Size matchBytes(const Match& m)
  {
    const auto& d = m.data;
    Size bytes = 128 + d.representation.size() + d.name.size() + (d.formula ? d.formula->size() : 0) + m.scores.size() * 9 + IO::metadataBytes(d);
    for (const auto& i : d.identifiers)
      bytes += 16 + i.database.size() + i.accession.size();
    if (d.adduct) bytes += 32 + d.adduct->getName().size() + d.adduct->getEmpiricalFormula().toString().size();
    for (const auto& e : d.sequence_evidence)
      bytes += 40 + e.accession.size() + e.before.size() + e.after.size();
    for (const auto& a : d.peak_annotations)
      bytes += 32 + a.annotation.size();
    return bytes;
  }
  template<class Query>
  void writeQuery(IO::TableWriter& writer, const Query& q, const IO::Dictionary& dictionary)
  {
    append<arrow::UInt64Builder>(writer.column(0), q.query_id);
    append<arrow::UInt32Builder>(writer.column(1), q.source_id);
    appendText(writer.column(2), q.data.data_id);
    appendOptional<arrow::DoubleBuilder>(writer.column(3), q.data.rt);
    appendOptional<arrow::DoubleBuilder>(writer.column(4), q.data.mz);
    appendOptional<arrow::UInt64Builder>(writer.column(5), q.selected_match_id);
    IO::appendMetadata(writer.column(6), q.data, dictionary);
    writer.finishRow(queryBytes(q));
  }
  template<class Match>
  void writeMatch(IO::TableWriter& writer, const Match& m, const IO::Dictionary& dictionary)
  {
    const auto& d = m.data;
    append<arrow::UInt64Builder>(writer.column(0), m.match_id);
    append<arrow::UInt64Builder>(writer.column(1), m.query_id);
    appendText(writer.column(2), d.representation);
    append<arrow::UInt8Builder>(writer.column(3), static_cast<unsigned>(d.encoding));
    append<arrow::Int32Builder>(writer.column(4), d.charge);
    appendOptional<arrow::DoubleBuilder>(writer.column(5), d.calculated_mz);
    append<arrow::UInt8Builder>(writer.column(6), static_cast<unsigned>(d.target_decoy));
    appendText(writer.column(7), d.name);
    if (d.formula) appendText(writer.column(8), *d.formula);
    else
      check(writer.column(8).AppendNull());
    auto& identifiers = static_cast<arrow::ListBuilder&>(writer.column(9));
    check(identifiers.Append());
    auto& identity = *static_cast<arrow::StructBuilder*>(identifiers.value_builder());
    for (const auto& i : d.identifiers)
    {
      check(identity.Append());
      appendText(*identity.field_builder(0), i.database);
      appendText(*identity.field_builder(1), i.accession);
    }
    auto& adduct = static_cast<arrow::StructBuilder&>(writer.column(10));
    if (! d.adduct) check(adduct.AppendNull());
    else
    {
      check(adduct.Append());
      appendText(*adduct.field_builder(0), d.adduct->getName());
      appendText(*adduct.field_builder(1), d.adduct->getEmpiricalFormula().toString());
      append<arrow::Int32Builder>(*adduct.field_builder(2), d.adduct->getCharge());
      append<arrow::UInt32Builder>(*adduct.field_builder(3), d.adduct->getMolMultiplier());
    }
    auto& evidence = static_cast<arrow::ListBuilder&>(writer.column(11));
    check(evidence.Append());
    auto& entry = *static_cast<arrow::StructBuilder*>(evidence.value_builder());
    for (const auto& e : d.sequence_evidence)
    {
      check(entry.Append());
      append<arrow::UInt32Builder>(*entry.field_builder(0), e.database.value);
      appendText(*entry.field_builder(1), e.accession);
      appendOptional<arrow::UInt64Builder>(*entry.field_builder(2), e.start);
      appendOptional<arrow::UInt64Builder>(*entry.field_builder(3), e.end);
      appendText(*entry.field_builder(4), e.before);
      appendText(*entry.field_builder(5), e.after);
    }
    auto& annotations = static_cast<arrow::ListBuilder&>(writer.column(12));
    check(annotations.Append());
    auto& annotation = *static_cast<arrow::StructBuilder*>(annotations.value_builder());
    for (const auto& a : d.peak_annotations)
    {
      check(annotation.Append());
      appendText(*annotation.field_builder(0), a.annotation);
      append<arrow::Int32Builder>(*annotation.field_builder(1), a.charge);
      append<arrow::DoubleBuilder>(*annotation.field_builder(2), a.mz);
      append<arrow::DoubleBuilder>(*annotation.field_builder(3), a.intensity);
    }
    IO::appendMetadata(writer.column(13), d, dictionary);
    // The shared table has the dataset's score columns; scoreless catalog runs leave them null.
    const Size score_columns = writer.columns() - 14;
    if (m.scores.size() > score_columns) invalid("Match has more scores than the dataset score schema");
    for (Size i = 0; i < score_columns; ++i)
    {
      if (i < m.scores.size()) appendOptional<arrow::DoubleBuilder>(writer.column(14 + i), storedScore(m.scores[i]));
      else
        check(writer.column(14 + i).AppendNull());
    }
    writer.finishRow(matchBytes(m));
  }
  template<class Callback>
  void readList(const arrow::Array& array, int64_t row, const Callback& callback)
  {
    if (array.IsNull(row)) invalid("Null required list");
    const auto& list = static_cast<const arrow::ListArray&>(array);
    const auto& entries = static_cast<const arrow::StructArray&>(*list.values());
    for (int64_t i = list.value_offset(row), end = i + list.value_length(row); i < end; ++i)
    {
      if (entries.IsNull(i)) invalid("Null list item");
      callback(entries, i);
    }
  }
  File::QueryRecord readQuery(const IO::TableReader& reader, const IO::Dictionary& dictionary)
  {
    const auto row = reader.row();
    File::QueryRecord q;
    q.query_id = number<arrow::UInt64Array>(reader.column(Size {0}), row);
    q.source_id = number<arrow::UInt32Array>(reader.column(Size {1}), row);
    q.data.data_id = text(reader.column(Size {2}), row);
    q.data.rt = optionalNumber<arrow::DoubleArray>(reader.column(Size {3}), row);
    q.data.mz = optionalNumber<arrow::DoubleArray>(reader.column(Size {4}), row);
    q.selected_match_id = optionalNumber<arrow::UInt64Array>(reader.column(Size {5}), row);
    if (reader.hasColumn(Size {6})) IO::readMetadata(reader.column(Size {6}), row, q.data, dictionary);
    if ((q.data.rt && ! std::isfinite(*q.data.rt)) || (q.data.mz && ! std::isfinite(*q.data.mz))) invalid("Nonfinite observation coordinate");
    return q;
  }
  File::MatchRecord readMatch(const IO::TableReader& reader, const IO::Dictionary& dictionary, Size score_count)
  {
    const auto row = reader.row();
    File::MatchRecord m;
    m.match_id = number<arrow::UInt64Array>(reader.column(Size {0}), row);
    m.query_id = number<arrow::UInt64Array>(reader.column(Size {1}), row);
    auto& d = m.data;
    if (reader.hasColumn(Size {2}))
    {
      d.representation = text(reader.column(Size {2}), row);
      auto encoding = number<arrow::UInt8Array>(reader.column(Size {3}), row);
      if (encoding > static_cast<unsigned>(ID::Encoding::DATABASE_ID)) invalid("Unknown molecular encoding");
      d.encoding = static_cast<ID::Encoding>(encoding);
      d.charge = number<arrow::Int32Array>(reader.column(Size {4}), row);
      d.calculated_mz = optionalNumber<arrow::DoubleArray>(reader.column(Size {5}), row);
      if (d.calculated_mz && ! std::isfinite(*d.calculated_mz)) invalid("Nonfinite calculated m/z");
      auto target_decoy = number<arrow::UInt8Array>(reader.column(Size {6}), row);
      if (target_decoy > static_cast<unsigned>(ID::TargetDecoy::BOTH)) invalid("Unknown target/decoy state");
      d.target_decoy = static_cast<ID::TargetDecoy>(target_decoy);
      d.name = text(reader.column(Size {7}), row);
      if (! reader.column(Size {8}).IsNull(row)) d.formula = text(reader.column(Size {8}), row);
      readList(reader.column(Size {9}), row,
               [&](const arrow::StructArray& a, int64_t i) { d.identifiers.push_back({text(*a.field(0), i), text(*a.field(1), i)}); });
      const auto& adduct = static_cast<const arrow::StructArray&>(reader.column(Size {10}));
      if (! adduct.IsNull(row))
      {
        d.adduct.emplace(text(*adduct.field(0), row), EmpiricalFormula(text(*adduct.field(1), row)), number<arrow::Int32Array>(*adduct.field(2), row),
                         number<arrow::UInt32Array>(*adduct.field(3), row));
        if (d.adduct->getCharge() != d.charge) invalid("Adduct and molecular charge disagree");
      }
    }
    if (reader.hasColumn(Size {11}))
      readList(reader.column(Size {11}), row, [&](const arrow::StructArray& a, int64_t i) {
        ID::SequenceEvidence e;
        e.database = {number<arrow::UInt32Array>(*a.field(0), i)};
        e.accession = text(*a.field(1), i);
        e.start = optionalNumber<arrow::UInt64Array>(*a.field(2), i);
        e.end = optionalNumber<arrow::UInt64Array>(*a.field(3), i);
        if (e.start && e.end && *e.start > *e.end) invalid("Reversed sequence evidence positions");
        e.before = text(*a.field(4), i);
        e.after = text(*a.field(5), i);
        d.sequence_evidence.push_back(std::move(e));
      });
    if (reader.hasColumn(Size {12}))
      readList(reader.column(Size {12}), row, [&](const arrow::StructArray& a, int64_t i) {
        PeptideHit::PeakAnnotation item;
        item.annotation = text(*a.field(0), i);
        item.charge = number<arrow::Int32Array>(*a.field(1), i);
        item.mz = number<arrow::DoubleArray>(*a.field(2), i);
        item.intensity = number<arrow::DoubleArray>(*a.field(3), i);
        d.peak_annotations.push_back(std::move(item));
      });
    if (reader.hasColumn(Size {13})) IO::readMetadata(reader.column(Size {13}), row, d, dictionary);
    m.scores.resize(score_count);
    for (Size i = 0; i < score_count; ++i)
    {
      if (reader.hasColumn(14 + i))
      {
        m.scores[i] = optionalNumber<arrow::DoubleArray>(reader.column(14 + i), row);
        if (m.scores[i] && ! std::isfinite(*m.scores[i])) invalid("Nonfinite score");
      }
    }
    return m;
  }
  Json sourceJson(const ID::SourceFile& source)
  {
    return {{"identifier", source.identifier}, {"path", source.path}, {"metadata", IO::metadataJson(source)}};
  }
  Json databaseJson(const ID::Database& database)
  {
    return {{"path", database.path}, {"version", database.version}, {"taxonomy", database.taxonomy}, {"metadata", IO::metadataJson(database)}};
  }
  ID::Database readDatabaseJson(const Json& j)
  {
    ID::Database database;
    database.path = j.at("path").get<std::string>();
    database.version = j.at("version").get<std::string>();
    database.taxonomy = j.at("taxonomy").get<std::string>();
    IO::readMetadataJson(j.at("metadata"), database);
    return database;
  }
  ID::SourceFile readSourceJson(const Json& j)
  {
    ID::SourceFile source;
    source.identifier = j.at("identifier").get<std::string>();
    source.path = j.at("path").get<std::string>();
    IO::readMetadataJson(j.at("metadata"), source);
    return source;
  }
  Json runJson(const ID::Run& run)
  {
    Json scores = Json::array(), sources = Json::array();
    for (const auto& s : run.getScoreDefinitions())
      scores.push_back(IO::scoreJson(s));
    for (const auto& s : run.getSources())
      sources.push_back(sourceJson(s.file));
    Json databases = Json::array();
    for (const auto& database : run.getDatabases())
      databases.push_back(databaseJson(database));
    return {{"identifier", run.getIdentifier()},
            {"uuid", run.getUuid()},
            {"molecule_kind", run.getMoleculeKind()},
            {"settings", IO::settingsJson(run.getSettings())},
            {"scores", std::move(scores)},
            {"sources", std::move(sources)},
            {"databases", std::move(databases)},
            {"primary_score", run.getPrimaryScore() ? Json(run.getPrimaryScore()->value) : Json()},
            {"next_query_id", run.getNextQueryId()},
            {"next_match_id", run.getNextMatchId()},
            {"query_count", run.getNumberOfIdentifications()},
            {"match_count", run.getNumberOfMatches()}};
  }
  ID::Run runShell(const Json& j)
  {
    unsigned kind = IO::integer<unsigned>(j.at("molecule_kind"));
    if (kind > static_cast<unsigned>(ID::MoleculeKind::COMPOUND)) invalid("Unknown molecule kind");
    ID::Run run(j.at("identifier").get<std::string>(), static_cast<ID::MoleculeKind>(kind));
    run.setSettings(IO::readSettingsJson(j.at("settings")));
    if (! j.at("scores").is_array() || ! j.at("sources").is_array() || ! j.at("databases").is_array()) invalid("Run descriptors must be arrays");
    for (const auto& s : j.at("scores"))
      run.addScore(IO::readScoreJson(s));
    for (const auto& s : j.at("sources"))
      run.addSource(readSourceJson(s));
    for (const auto& database : j.at("databases"))
      if (run.addDatabase(readDatabaseJson(database)).value != run.getDatabases().size() - 1) invalid("Duplicate database of a run");
    if (! j.at("primary_score").is_null()) run.setPrimaryScore(run.getScoreId(IO::integer<UInt32>(j.at("primary_score"))));
    return run;
  }
  File::RunDescriptor descriptor(const Json& j)
  {
    auto shell = runShell(j);
    shell.restoreIdentity(j.at("uuid").get<std::string>(), IO::integer<UInt64>(j.at("next_query_id")), IO::integer<UInt64>(j.at("next_match_id")));
    File::RunDescriptor d;
    d.identifier = shell.getIdentifier();
    d.uuid = shell.getUuid();
    d.molecule_kind = shell.getMoleculeKind();
    d.scores = shell.getScoreDefinitions();
    for (const auto& source : shell.getSources())
      d.sources.push_back(source.file);
    d.databases = shell.getDatabases();
    if (shell.getPrimaryScore()) d.primary_score = shell.getPrimaryScore()->value;
    d.next_query_id = shell.getNextQueryId();
    d.next_match_id = shell.getNextMatchId();
    d.query_count = IO::integer<UInt64>(j.at("query_count"));
    d.match_count = IO::integer<UInt64>(j.at("match_count"));
    return d;
  }
  std::vector<const Json*> selectRuns(const Json& manifest, const std::vector<std::string>& selection)
  {
    std::set<std::string> selected;
    for (const auto& key : selection)
    {
      Size found = 0;
      for (const auto& j : manifest.at("runs"))
        if (j.at("identifier") == key || j.at("uuid") == key)
        {
          selected.insert(j.at("uuid").get<std::string>());
          ++found;
        }
      if (found != 1) invalid("Unknown or ambiguous run selection: " + key);
    }
    std::vector<const Json*> runs;
    for (const auto& j : manifest.at("runs"))
      if (selection.empty() || selected.contains(j.at("uuid").get<std::string>())) runs.push_back(&j);
    return runs;
  }
  void collectDictionary(const ID::Run& run, IO::Dictionary& dictionary)
  {
    for (const auto& source : run.getSources())
      for (const auto& query : source.identifications)
      {
        dictionary.collect(query);
        for (const auto& match : query.getMatches())
          dictionary.collect(match);
      }
    if (run.getDatabaseSequences())
      for (const auto& sequence : *run.getDatabaseSequences())
        dictionary.collect(sequence);
  }
  std::vector<std::string> projectedMatchColumns(const File::Projection& projection, const std::vector<std::string>& score_columns)
  {
    std::vector<std::string> fields {"match_id", "query_id"};
    if (projection.molecule)
      for (const auto* field : {"representation", "encoding", "charge", "calculated_mz", "target_decoy", "name", "formula", "identifiers", "adduct"})
        fields.emplace_back(field);
    if (projection.evidence) fields.emplace_back("sequence_evidence");
    if (projection.annotations) fields.emplace_back("peak_annotations");
    if (projection.metadata) fields.emplace_back("metadata");
    std::set<UInt32> scores(projection.score_ids.begin(), projection.score_ids.end());
    for (UInt32 score : scores)
      if (score >= score_columns.size()) invalid("Projected score ID does not exist in the run");
    for (Size i = 0; i < score_columns.size(); ++i)
      if (projection.all_scores || scores.contains(static_cast<UInt32>(i))) fields.push_back(score_columns[i]);
    return fields;
  }
  // Both tables retain scientific ordering. A single match row is sufficient lookahead.
  template<class Begin, class Match, class End>
  void walkRun(const fs::path& root, const Json& j, const File::ScanOptions& options, const Begin& begin, const Match& consume, const End& end, const IO::Options& io)
  {
    const auto d = descriptor(j);
    IO::Dictionary dictionary;
    dictionary.load(j.at("metadata_descriptors"));
    const auto& tables = j.at("tables");
    std::vector<std::string> query_columns {"query_id", "source_id", "data_id", "rt", "mz", "selected_match_id"};
    if (options.projection.metadata) query_columns.emplace_back("metadata");
    IO::TableReader queries(root, tables.at("queries"), querySchema(), io, query_columns);
    IO::TableReader matches(root, tables.at("matches"), matchSchema(io.score_columns), io, projectedMatchColumns(options.projection, io.score_columns));
    if (queries.rows() != d.query_count || matches.rows() != d.match_count) invalid("Manifest and table row counts disagree");
    std::set<UInt64> query_ids, match_ids;
    bool has_match = matches.next();
    File::MatchRecord pending;
    if (has_match) pending = readMatch(matches, dictionary, d.scores.size());
    UInt64 seen_queries = 0, seen_matches = 0;
    UInt32 prior_source = 0;
    while (queries.next())
    {
      File::QueryRecord q = readQuery(queries, dictionary);
      ++seen_queries;
      if (! q.query_id || q.query_id >= d.next_query_id || q.source_id >= d.sources.size() || q.source_id < prior_source)
        invalid("Invalid query ID, source or source ordering");
      prior_source = q.source_id;
      if (q.selected_match_id && (! *q.selected_match_id || *q.selected_match_id >= d.next_match_id)) invalid("Invalid selected match ID");
      if (queryBytes(q) > options.buffering.max_record_bytes) invalid("Query exceeds max_record_bytes");
      if (options.validate_unique_ids && ! query_ids.insert(q.query_id).second) invalid("Duplicate query ID");
      begin(q);
      bool found_selected = ! q.selected_match_id;
      while (has_match && pending.query_id == q.query_id)
      {
        if (! pending.match_id || pending.match_id >= d.next_match_id) invalid("Invalid match ID or allocation counter");
        if (options.projection.molecule)
        {
          const auto encoding = pending.data.encoding;
          const bool compatible
            = encoding == ID::Encoding::DATABASE_ID || (d.molecule_kind == ID::MoleculeKind::PEPTIDE && encoding == ID::Encoding::AA_SEQUENCE)
              || (d.molecule_kind == ID::MoleculeKind::OLIGONUCLEOTIDE && encoding == ID::Encoding::NA_SEQUENCE)
              || (d.molecule_kind == ID::MoleculeKind::COMPOUND && (encoding == ID::Encoding::SMILES || encoding == ID::Encoding::INCHI));
          if (pending.data.representation.empty() || ! compatible) invalid("Empty or incompatible molecular representation");
        }
        if (d.molecule_kind == ID::MoleculeKind::COMPOUND && ! pending.data.sequence_evidence.empty())
          invalid("Compound candidate has sequence evidence");
        for (const auto& evidence : pending.data.sequence_evidence)
        {
          if (evidence.accession.empty()) invalid("Sequence evidence needs an accession");
          if (evidence.database.value >= d.databases.size()) invalid("Sequence evidence refers to an unknown database of its run");
        }
        if (matchBytes(pending) > options.buffering.max_record_bytes) invalid("Match exceeds max_record_bytes");
        if (options.validate_unique_ids && ! match_ids.insert(pending.match_id).second) invalid("Duplicate match ID");
        if (d.primary_score)
        {
          if (matches.hasColumn(14 + *d.primary_score) && ! pending.scores[*d.primary_score]) invalid("Primary score is missing");
        }
        found_selected = found_selected || q.selected_match_id == pending.match_id;
        consume(pending);
        ++seen_matches;
        has_match = matches.next();
        if (has_match) pending = readMatch(matches, dictionary, d.scores.size());
      }
      if (! found_selected) invalid("Selected match does not belong to its query");
      end(q);
    }
    if (has_match || seen_queries != d.query_count || seen_matches != d.match_count) invalid("Unmatched or out-of-order match rows");
  }
  ID::Run readRun(const fs::path& root, const Json& j, const IO::Options& options)
  {
    auto run = runShell(j);
    File::ScanOptions scan;
    scan.buffering = options;
    walkRun(
      root, j, scan,
      [&](File::QueryRecord& q) { run.importIdentification(run.getSourceId(q.source_id), ID::QueryId {q.query_id}, std::move(q.data)); },
      [&](File::MatchRecord& m) {
        // Runs without score definitions (scoreless catalogs) leave the shared score columns null.
        const Size declared = run.getScoreDefinitions().size();
        if (m.scores.size() > declared)
        {
          if (std::any_of(m.scores.begin() + declared, m.scores.end(), [](const auto& value) { return value.has_value(); }))
            invalid("Score value outside the run's declared score definitions");
          m.scores.resize(declared);
        }
        run.importMatch(ID::QueryId {m.query_id}, ID::MatchId {m.match_id}, std::move(m.data), m.scores);
      },
      [&](const File::QueryRecord& q) {
        if (q.selected_match_id) run.setSelectedMatch(ID::QueryId {q.query_id}, ID::MatchId {*q.selected_match_id});
      },
      options);
    const auto& tables = j.at("tables");
    if (tables.contains("database_sequences"))
    {
      IO::Dictionary dictionary;
      dictionary.load(j.at("metadata_descriptors"));
      run.setDatabaseSequences(IO::readDatabaseSequences(root, tables.at("database_sequences"), dictionary, options));
    }
    run.restoreIdentity(j.at("uuid").get<std::string>(), IO::integer<UInt64>(j.at("next_query_id")), IO::integer<UInt64>(j.at("next_match_id")));
    return run;
  }
} // namespace

bool File::isNativeFile(const std::string& path)
{
  const auto manifest = fs::path(path) / "manifest.json";
  std::error_code ec;
  if (! fs::is_regular_file(manifest, ec)) return false;
  std::ifstream in(manifest, std::ios::binary);
  if (! in) invalid("Cannot open manifest.json");
  Json j;
  try
  {
    in >> j;
  }
  catch (const Json::exception& e)
  {
    invalid(e.what());
  }
  return j.contains("format") && j.at("format") == FORMAT;
}
std::vector<File::RunDescriptor> File::inspect(const std::string& path)
try
{
  const Json manifest = readManifest(path);
  std::vector<RunDescriptor> result;
  for (const auto& j : manifest.at("runs"))
  {
    result.push_back(descriptor(j));
    if (! result.back().scores.empty()) result.back().score_columns = scoreColumns(manifest);
  }
  return result;
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
void File::store(const std::string& path, const ID& data)
{ store(path, data, Options {}); }
void File::store(const std::string& path, const ID& data, const Options& options)
{
  IO::validateOptions(options);
  data.validate();
  StagedDirectory output(path, options.replace_existing);
  IO::Options io(options);
  io.output = std::make_shared<IO::WritePool>(output.path);
  io.score_columns = IO::scoreColumns(data.getScoreDefinitions());
  Json manifest {
    {"format", FORMAT}, {"schema_version", 1}, {"score_columns", io.score_columns}, {"runs", Json::array()}, {"inference", Json::array()}};
  for (const auto& run : data.getRuns())
  {
    io.partition = run.getUuid();
    IO::Dictionary dictionary;
    collectDictionary(run, dictionary);
    Json j = runJson(run);
    Json tables {{"queries", "queries.parquet"}, {"matches", "matches.parquet"}};
    IO::TableWriter queries(output.path / tables.at("queries").get<std::string>(), querySchema(), io);
    IO::TableWriter matches(output.path / tables.at("matches").get<std::string>(), matchSchema(io.score_columns), io);
    for (const auto& source : run.getSources())
      for (const auto& query : source.identifications)
      {
        QueryView q {query.getId().value, source.id.value, query.getObservation(), std::nullopt};
        if (query.getSelectedMatch()) q.selected_match_id = query.getSelectedMatch()->value;
        writeQuery(queries, q, dictionary);
        for (const auto& match : query.getMatches())
          writeMatch(matches, MatchView {match.getId().value, query.getId().value, match.getData(), match.getScoreValues()}, dictionary);
      }
    queries.close();
    matches.close();
    tables["queries"] = queries.reference();
    tables["matches"] = matches.reference();
    if (run.getDatabaseSequences())
    {
      tables["database_sequences"] = IO::writeDatabaseSequences(output.path / "database_sequences.parquet", *run.getDatabaseSequences(), dictionary, io);
    }
    j["metadata_descriptors"] = dictionary.toJson();
    j["tables"] = std::move(tables);
    manifest["runs"].push_back(std::move(j));
  }
  for (const auto& result : data.getInferenceResults())
  {
    io.partition = result.identifier;
    Json j = IO::writeInference(output.path, result, io);
    manifest["inference"].push_back(std::move(j));
  }
  io.output->close();
  writeManifest(output.path, manifest);
  output.publish();
}
void File::load(const std::string& path, ID& data)
{ load(path, data, Options {}); }
void File::load(const std::string& path, ID& data, const Options& options)
try
{
  IO::validateOptions(options);
  const auto manifest = readManifest(path);
  IO::Options io(options);
  io.input = std::make_shared<IO::ReadPool>();
  io.score_columns = scoreColumns(manifest);
  ID temporary;
  for (const auto& j : manifest.at("runs"))
    temporary.addRun(readRun(path, j, io));
  for (const auto& j : manifest.at("inference"))
  {
    auto result = IO::readInference(path, j, io);
    temporary.addInferenceResult(std::move(result));
  }
  // addRun validates each run and the dataset score/identity contract; adding
  // inference validates its provenance without changing those runs. Do not repeat
  // the full per-record validation after these checked construction boundaries.
  data.swap(temporary);
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
ID::Run File::loadRun(const std::string& path, const std::string& run)
{ return loadRun(path, run, Options {}); }
ID::Run File::loadRun(const std::string& path, const std::string& run, const Options& options)
try
{
  IO::validateOptions(options);
  const auto manifest = readManifest(path);
  IO::Options io(options);
  io.input = std::make_shared<IO::ReadPool>();
  io.score_columns = scoreColumns(manifest);
  const auto selected = selectRuns(manifest, {run});
  auto result = readRun(path, *selected.front(), io);
  result.validate();
  return result;
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
File::ScanStatistics File::scan(const std::string& path, const ScanOptions& options, const QueryCallback& on_queries, const MatchCallback& on_matches)
try
{
  IO::validateOptions(options.buffering);
  const auto manifest = readManifest(path);
  IO::Options io(options.buffering);
  io.input = std::make_shared<IO::ReadPool>();
  io.score_columns = scoreColumns(manifest);
  ScanStatistics statistics;
  std::error_code size_error;
  statistics.descriptor_bytes = fs::file_size(fs::path(path) / "manifest.json", size_error);
  if (size_error) invalid("Cannot read manifest.json: " + size_error.message());
  for (const auto* j : selectRuns(manifest, options.runs))
  {
    const std::string uuid = j->at("uuid").get<std::string>();
    std::vector<QueryRecord> queries;
    std::vector<MatchRecord> matches;
    Size query_bytes = 0, match_bytes = 0;
    auto flush_queries = [&]() {
      if (! queries.empty() && on_queries) on_queries(uuid, queries);
      queries.clear();
      query_bytes = 0;
    };
    auto flush_matches = [&]() {
      if (! matches.empty() && on_matches) on_matches(uuid, matches);
      matches.clear();
      match_bytes = 0;
    };
    walkRun(
      path, *j, options,
      [&](QueryRecord& q) {
        ++statistics.queries;
        if (on_queries)
        {
          query_bytes += queryBytes(q);
          queries.push_back(std::move(q));
          if (queries.size() >= options.buffering.batch_rows || query_bytes >= options.buffering.batch_bytes) flush_queries();
        }
      },
      [&](MatchRecord& m) {
        ++statistics.matches;
        if (on_matches)
        {
          match_bytes += matchBytes(m);
          matches.push_back(std::move(m));
          if (matches.size() >= options.buffering.batch_rows || match_bytes >= options.buffering.batch_bytes) flush_matches();
        }
      },
      [](const QueryRecord&) {}, io);
    flush_queries();
    flush_matches();
  }
  return statistics;
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
void File::filter(const std::string& input,
                  const std::string& output,
                  const MatchPredicate& keep,
                  ID::InferencePolicy policy,
                  bool keep_empty_queries)
{ filter(input, output, keep, policy, keep_empty_queries, Options {}); }
void File::filter(const std::string& input,
                  const std::string& output,
                  const MatchPredicate& keep,
                  ID::InferencePolicy policy,
                  bool keep_empty_queries,
                  const Options& options)
try
{
  IO::validateOptions(options);
  if (! keep) invalid("Filter requires a predicate");
  if (policy != ID::InferencePolicy::PRESERVE && policy != ID::InferencePolicy::DISCARD) invalid("Unknown inference policy");
  auto manifest = readManifest(input);
  StagedDirectory staged(output, options.replace_existing);
  IO::Options io(options);
  io.input = std::make_shared<IO::ReadPool>();
  io.score_columns = scoreColumns(manifest);
  io.output = std::make_shared<IO::WritePool>(staged.path);
  const auto copyTable = [&](const Json& reference) {
    const auto name = reference.at("path").get<std::string>();
    auto dest = IO::tablePath(staged.path, name);
    std::error_code ec;
    if (! fs::exists(dest, ec))
    {
      fs::create_directories(dest.parent_path(), ec);
      if (! ec) fs::copy_file(IO::tablePath(input, name), dest, ec);
      if (ec) fileError(dest, ec.message());
    }
  };
  for (auto& j : manifest["runs"])
  {
    const std::string uuid = j.at("uuid").get<std::string>();
    io.partition = uuid;
    IO::Dictionary dictionary;
    dictionary.load(j.at("metadata_descriptors"));
    IO::TableWriter queries(staged.path / "queries.parquet", querySchema(), io);
    IO::TableWriter matches(staged.path / "matches.parquet", matchSchema(io.score_columns), io);
    File::ScanOptions scan;
    scan.buffering = options;
    bool kept_any = false, kept_selection = false;
    std::optional<UInt64> selected;
    walkRun(
      input, j, scan,
      [&](const QueryRecord& q) {
        kept_any = false;
        kept_selection = false;
        selected = q.selected_match_id;
      },
      [&](const MatchRecord& m) {
        if (keep(uuid, m))
        {
          writeMatch(matches, m, dictionary);
          kept_any = true;
          if (selected == m.match_id) kept_selection = true;
        }
      },
      [&](const QueryRecord& original) {
        QueryRecord q = original;
        if (q.selected_match_id && !kept_selection) q.selected_match_id.reset();
        if (kept_any || keep_empty_queries) writeQuery(queries, q, dictionary);
      }, io);
    queries.close();
    matches.close();
    j["query_count"] = queries.rows();
    j["match_count"] = matches.rows();
    j["tables"]["queries"] = queries.reference();
    j["tables"]["matches"] = matches.reference();
    if (j.at("tables").contains("database_sequences"))
    {
      const auto& reference = j.at("tables").at("database_sequences");
      IO::validateDatabaseSequences(input, reference, j.at("databases").size(), dictionary, io);
      copyTable(reference);
    }
  }
  if (policy == ID::InferencePolicy::DISCARD) manifest["inference"] = Json::array();
  else
    for (const auto& result : manifest.at("inference"))
    {
      IO::validateInferenceTables(input, result, io);
      for (const auto& reference : result.at("tables")) copyTable(reference);
    }
  io.output->close();
  writeManifest(staged.path, manifest);
  staged.publish();
}
catch (const Json::exception& error)
{
  invalid(error.what());
}
} // namespace OpenMS
